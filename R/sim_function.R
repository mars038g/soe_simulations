library(dplyr)
library(AICcmodavg)
library(nlme)
library(Kendall)
library(zyp)

load("sim_param.rdata")

test.series <- function(ar.order,
                        rho,
                        series.length,
                        trend.strength,
                        var,
                        nsims){
  
  #Get true trend
  b <- trend.df %>% filter(Var == trend.strength) %>% pull(Value)
  true_trend <- (b* c(1:series.length))

  #AR parameters
  order <- c(ar.order,0,0)
  
  if (ar.order == 1){
    ar.strength <- get(rho)
  } else {
    rho <- NA
    ar.strength <- c(mean_ar2_coef1,mean_ar2_coef2)
    var <- iv_mean_ar2
  }
  
  #Simulate
  if (!is.na(rho) && rho == "noAR"){
    dat <- arima.sim(list(ar = list()),
                     n=series.length,
                     rand.gen=rnorm,
                     sd = sqrt(var),
                     n.start = NA) #Gives reasonable length of burn-in period
    
  } else {
    dat <- arima.sim(list(order = order,
                          ar = ar.strength),
                     n=series.length,
                     rand.gen=rnorm,
                     sd = sqrt(var),
                     n.start = NA) #Gives reasonable length of burn-in period
  }
  

  #add autocorrelated error structure to trend
  dat <- data.frame(series = true_trend + dat,
                    time = 1:length(dat))
  
  #fit gls model
  gls_sim <- fit_lm(dat = dat, ar.order = ar.order)
  gls_chosen <- gls_sim$best_lm$model
    
  if (is.na(gls_sim[1])){
    gls_mae <- NA
    gls_rmse <- NA
    pval <- NA
    slope_pred <- NA
    slope_true <- b
    
  } else {
    newtime <- seq(1, series.length, 1)
    newdata <- data.frame(time = newtime,
                          time2 = newtime^2)
    gls_pred <- AICcmodavg::predictSE(gls_sim$model,
                                      newdata = newdata,
                                      se.fit = F)
    # if(!is.na(gls_sim$best_lm$coefs.time2)){ 
    #   slope_pred <- NA
    # } else {
    slope_pred <- gls_pred[2] - gls_pred[1]
    # }
    
    ##Get error and pvalue
    gls_rmse <- sqrt(mean((gls_pred - true_trend)^2))
    gls_mae <- mean(abs(gls_pred - true_trend))
    pval <- gls_sim$best_lm$pval
  }
  
    #Results DF
    gls_df <- data.frame(test = "gls",
                         series.length = series.length,
                         n = nsims,
                         gls_chosen = gls_chosen,
                         rho1 = ar.strength[1],
                         rho2 = ifelse(ar.order == 2, ar.strength[2], NA),
                         trend = trend.strength,
                         ar = ar.order,
                         mae = gls_mae,
                         rmse = gls_rmse,
                         p = pval,
                         slope_pred = slope_pred,
                         slope_true = b,
                         var = var)
  
  
    #---------------------------------MK---------------------------------#
    mk <- MannKendall(dat$series)
    mk_p <- unlist(mk[2])
    
    mk_df <- data.frame(test = "mk",
                        series.length = series.length,
                        n = nsims,
                        gls_chosen = NA,
                        rho1 = ar.strength[1],
                        rho2 = ifelse(ar.order == 2, ar.strength[2], NA),
                        trend = trend.strength,
                        ar = ar.order,
                        mae = NA,
                        rmse = NA,
                        p = mk_p,
                        slope_pred = NA,
                        slope_true = NA,
                        var = var)
    
    #---------------------------------MK-TFPW---------------------------------#
    pw <- zyp.trend.vector(dat$series,method='yuepilon')
    
    pw_pred <- pw[[11]] + 1:series.length * pw[[2]]
    
    pw_rmse <- sqrt(mean((pw_pred - true_trend)^2))
    pw_mae <- mean(abs(pw_pred - true_trend))
    pw_p <- pw[6]
    slope_pred <- pw[[2]]
    
    pw_df <- data.frame(test = "pw",
                        series.length = series.length,
                        n = nsims,
                        gls_chosen = NA,
                        rho1 = ar.strength[1],
                        rho2 = ifelse(ar.order == 2, ar.strength[2], NA),
                        trend = trend.strength,
                        ar = ar.order,
                        mae = pw_mae,
                        rmse = pw_rmse,
                        p = pw_p,
                        slope_pred = slope_pred,
                        slope_true = b,
                        var = var)
    
    int_df <- rbind(gls_df, mk_df, pw_df)
    
    return(int_df)
}


