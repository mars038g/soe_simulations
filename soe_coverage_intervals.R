#Mann kendall with pre-whitening coverage intervals
library(boot);library(Kendall)
library(zoo);library(zyp)
library(trend);library(dplyr)
library(AICcmodavg);library(nlme)
library(MASS)
setwd('z:/shardison/soe/simulations')
ptm <- proc.time()
fit_lm <- function(dat) {
  # Remove missing values first so that all models
  # use the same number of observations (important for AIC)
  dat <- dat %>% dplyr::filter(complete.cases(.))
  
  # Constant model (null model used to calculate 
  # overall p-value)
  constant_norm <-
    nlme::gls(series ~ 1, 
              data = dat)
  
  constant_ar1 <-
    try(nlme::gls(series ~ 1,
                  data = dat,
                  correlation = nlme::corAR1(form = ~time)))
  if (class(constant_ar1) == "try-error"){
    return(best_lm <- data.frame(model = NA,
                                 aicc  = NA,
                                 coefs..Intercept = NA,
                                 coefs.time = NA,
                                 coefs.time2 = NA,
                                 pval = NA))
    
  } 
  
  
  
  # Linear model with normal error
  linear_norm <- 
    nlme::gls(series ~ time, 
              data = dat)
  
  # Linear model with AR1 error
  linear_ar1 <- 
    try(nlme::gls(series ~ time, 
                  data = dat,
                  correlation = nlme::corAR1(form = ~time)))
  if (class(linear_ar1) == "try-error"){
    return(best_lm <- data.frame(model = NA,
                                 aicc  = NA,
                                 coefs..Intercept = NA,
                                 coefs.time = NA,
                                 coefs.time2 = NA,
                                 pval = NA))
    
  }
  
  # Polynomial model with normal error
  dat$time2 <- dat$time^2
  poly_norm <- 
    nlme::gls(series ~ time + time2, 
              data = dat)
  
  # Polynomial model with AR1 error
  poly_ar1 <- 
    try(nlme::gls(series ~ time + time2, 
                  data = dat,
                  correlation = nlme::corAR1(form = ~time)))
  if (class(poly_ar1) == "try-error"){
    return(best_lm <- data.frame(model = NA,
                                 aicc  = NA,
                                 coefs..Intercept = NA,
                                 coefs.time = NA,
                                 coefs.time2 = NA,
                                 pval = NA))
    
  }
  
  # Calculate AICs for all models
  df_aicc <-
    data.frame(model = c("poly_norm",
                         "poly_ar1",
                         "linear_norm",
                         "linear_ar1"),
               aicc  = c(AICc(poly_norm),
                         AICc(poly_ar1),
                         AICc(linear_norm),
                         AICc(linear_ar1)),
               coefs = rbind(coef(poly_norm),
                             coef(poly_ar1),
                             c(coef(linear_norm), NA),
                             c(coef(linear_ar1),  NA)),
               # Calculate overall signifiance (need to use
               # ML not REML for this)
               pval = c(anova(update(constant_norm, method = "ML"),
                              update(poly_norm, method = "ML"))$`p-value`[2],
                        anova(update(constant_ar1, method = "ML"),
                              update(poly_ar1, method = "ML"))$`p-value`[2],
                        anova(update(constant_norm, method = "ML"),
                              update(linear_norm, method = "ML"))$`p-value`[2],
                        anova(update(constant_ar1, method = "ML"),
                              update(linear_ar1, method = "ML"))$`p-value`[2]))
  
  best_lm <-
    df_aicc %>%
    dplyr::filter(aicc == min(aicc))
  
  
  if (best_lm$model == "poly_norm") {
    model <- poly_norm
  } else if (best_lm$model == "poly_ar1") {
    model <- poly_ar1
  } else if (best_lm$model == "linear_norm") {
    model <- linear_norm
  } else if (best_lm$model == "linear_ar1") {
    model <- linear_ar1
  }
  
  return(list(best_lm = best_lm, 
              model = model))
}

n = 500 #number of simulations
ARsd <- .54^.5 #standard deviation


#set phi
NOAR <- list()
weakAR <- 0.1
medAR <- 0.433
strongAR <- 0.8

#place holders for generated data
mk.ts.NOAR.ltrendweak <- NULL
mk.ts.NOAR.notrend <- NULL
mk.ts.NOAR.ltrendmed <- NULL
mk.ts.NOAR.ltrendstrong <- NULL

mk.ts.medAR.notrend <- NULL
mk.ts.medAR.ltrendweak <- NULL
mk.ts.medAR.ltrendmed <- NULL
mk.ts.medAR.ltrendstrong <- NULL

mk.ts.strongAR.notrend <- NULL
mk.ts.strongAR.ltrendweak <- NULL
mk.ts.strongAR.ltrendmed <- NULL
mk.ts.strongAR.ltrendstrong <- NULL

pw.ts.NOAR.notrend <- NULL
pw.ts.NOAR.ltrendweak <- NULL
pw.ts.NOAR.ltrendmed <- NULL
pw.ts.NOAR.ltrendstrong <- NULL

pw.ts.medAR.notrend <- NULL
pw.ts.medAR.ltrendweak <- NULL
pw.ts.medAR.ltrendmed <- NULL
pw.ts.medAR.ltrendstrong <- NULL

pw.ts.strongAR.notrend <- NULL
pw.ts.strongAR.ltrendweak <- NULL
pw.ts.strongAR.ltrendmed <- NULL
pw.ts.strongAR.ltrendstrong <- NULL

gls.ts.NOAR.notrend <- NULL
gls.ts.NOAR.ltrendweak <- NULL
gls.ts.NOAR.ltrendmed <- NULL
gls.ts.NOAR.ltrendstrong <- NULL

gls.ts.medAR.notrend <- NULL
gls.ts.medAR.ltrendweak <- NULL
gls.ts.medAR.ltrendmed <- NULL
gls.ts.medAR.ltrendstrong <- NULL

gls.ts.strongAR.notrend <- NULL
gls.ts.strongAR.ltrendweak <- NULL
gls.ts.strongAR.ltrendmed <- NULL
gls.ts.strongAR.ltrendstrong <- NULL

index <- NULL
for (m in c(20)){
  notrend <- 0
  ltrendweak <- -0.262 + (0.004 * c(1:m)) 
  ltrendmed <- -0.262 + (0.051 * c(1:m)) 
  ltrendstrong <- -0.262 + (0.147 * c(1:m)) 
  print(paste("m=",m))
  for (k in c("notrend","ltrendweak","ltrendmed","ltrendstrong")){
    
    for (j in c("NOAR","medAR","strongAR")){
      
      #Generate time series confidence interval for each parameter
      og <- arima.sim(list(ar = get(j)), n=m, rand.gen=rnorm, sd = ARsd)
      og <- get(k) + og
      og <- data.frame(series = og,
                       time = 1:length(og))
      
      newtime <- seq(1, nrow(og), 1)
      newdata <- data.frame(time = newtime,
                            time2 = newtime^2)
      
      #GLS tryCatch structure - Will produce NAs for trendline if GLS throws an error
      gls_og <- tryCatch({
        lm_out_og <- fit_lm(dat = og)
        
        if (length(lm_out_og$model$coefficients) == 2){
          lm_out_og$model$coefficients[2] <- abs(lm_out_og$model$coefficients[2]) 
        } else {
          lm_out_og$model$coefficients[2] <- abs(lm_out_og$model$coefficients[2]) 
          lm_out_og$model$coefficients[3] <- abs(lm_out_og$model$coefficients[3]) 
        }
        
        gls_og <- AICcmodavg::predictSE(lm_out_og$model, 
                                        newdata = newdata,
                                        se.fit = T)
        
      } , 
      error = function(e) {
        print('error')
      }
      )
      
      if (gls_og == "error"){
        gls_og_fit <- rep(NA, m)
        gls_og_se <- rep(NA, m)
        ci_upr_gls <- rep(NA, m)
        ci_lwr_gls <- rep(NA, m)
      } else {
        gls_og_fit <- gls_og$fit
        gls_og_se <- gls_og$se.fit 
        ci_upr_gls <- gls_og_fit + 1.96*gls_og_se
        ci_lwr_gls <- gls_og_fit - 1.96*gls_og_se
      }
      
      #mann kendall
      mk_og <- zyp.sen(series~time, og)
      mk_og <- mk_og$coefficients[1]+og$time*abs(mk_og$coefficients[2])
      
      #mann kendall with pre-whitening
      pw_og <- zyp.trend.vector(og$series, method = "yuepilon")
      pw_og <- pw_og[11]+og$time*abs(pw_og[2])
      
      ##get confidence intervals from fits using a bootstrap
      #mann kendall
      tau_func <- function(z) MannKendall(z)$tau
      boot.out <- tsboot(og$series, tau_func, R=999, l=5, sim="fixed")
      ci_upr_mk <- boot.ci(boot.out, type="norm")$normal[3]
      ci_upr_mk <- mk_og + abs(ci_upr_mk)
      ci_lwr_mk <- boot.ci(boot.out, type="norm")$normal[2]
      ci_lwr_mk <- mk_og - abs(ci_lwr_mk)  
      
      #MK-PW
      pw_tau_func <- function(z) zyp.trend.vector(z, method = "yuepilon")[5]
      boot.out <- tsboot(og$series, pw_tau_func, R=999, l=5, sim="fixed")
      ci_upr_pw <- boot.ci(boot.out, type="norm")$normal[3]
      ci_upr_pw <- pw_og + abs(ci_upr_pw)
      ci_lwr_pw <- boot.ci(boot.out, type="norm")$normal[2]
      ci_lwr_pw <- pw_og - abs(ci_lwr_pw)  
      
      for (i in 1:n){
        #generate simulations
        dat <- arima.sim(list(ar = get(j)), n=m, rand.gen=rnorm, sd = ARsd)
        
        #add autocorrelated error structure to trend
        dat <- get(k) + dat
        
        #what fraction of the confidence intervals contain the original trendline?
        dat <- data.frame(series = dat,
                          time = 1:length(dat))
        
        #simulate trends
        pw_sim <- zyp.trend.vector(dat$series,method='yuepilon')
        mk_sim <- MannKendall(dat$series)
        
        #GLS tryCatch structure - Will produce NAs for trendline if GLS throws an error
        gls_sim <- tryCatch({
          gls_sim <- fit_lm(dat = dat)
          if (length(gls_sim$model$coefficients) == 2){
            gls_sim$model$coefficients[2] <- abs(gls_sim$model$coefficients[2]) 
          } else {
            gls_sim$model$coefficients[2] <- abs(gls_sim$model$coefficients[2]) 
            gls_sim$model$coefficients[3] <- abs(gls_sim$model$coefficients[3]) 
          }
          gls_sim <- AICcmodavg::predictSE(gls_sim$model, 
                                          newdata = newdata,
                                          se.fit = T)
        }, 
        error = function(e) {
          print('error')
        })
        
        if (gls_sim == "error"){
          gls_pred <- rep(NA, m)
        } else {
          gls_pred <- gls_sim$fit #GLS trend
        }
        
        #get model fits from MK and MK-PW
        pw_pred <- pw_sim[11]+dat$time*abs(pw_sim[2])
        mk_pred <- zyp.sen(series~time, dat)
        mk_pred <- mk_pred$coefficients[1]+dat$time*abs(mk_pred$coefficients[2])


        # collect the simulations inside (1) and outside (0) the coverage intervals
        for (g in 1:m){
          
          if (is.na(gls_pred[g]) | is.na(ci_lwr_gls[g]) | is.na(ci_upr_gls[g])){
            assign(paste0("gls.ts.",j,".",k),rbind(get(paste0("gls.ts.",j,".",k)),100))
          } else if ((gls_pred[g] < ci_upr_gls[g]) & (gls_pred[g] > ci_lwr_gls[g])){
            assign(paste0("gls.ts.",j,".",k),rbind(get(paste0("gls.ts.",j,".",k)),1))
          } else if ((gls_pred[g] > ci_upr_gls[g]) | (gls_pred[g] < ci_lwr_gls[g])){
            assign(paste0("gls.ts.",j,".",k),rbind(get(paste0("gls.ts.",j,".",k)),0))
          } 
          
          if ((mk_pred[g] < ci_upr_mk[g]) & (mk_pred[g] > ci_lwr_mk[g])){
            assign(paste0("mk.ts.",j,".",k),rbind(get(paste0("mk.ts.",j,".",k)),1))
          } else {
            assign(paste0("mk.ts.",j,".",k),rbind(get(paste0("mk.ts.",j,".",k)),0))
          }
          
          if ((pw_pred[g] < ci_upr_pw[g]) & (pw_pred[g] > ci_lwr_pw[g])){
            assign(paste0("pw.ts.",j,".",k),rbind(get(paste0("pw.ts.",j,".",k)),1))
          } else {
            assign(paste0("pw.ts.",j,".",k),rbind(get(paste0("pw.ts.",j,".",k)),0))
          }
          
        }
        
        
      }
      
    } 
    
  }
  sim_results = data.frame(mk.NOAR.ltrendweak = mk.ts.NOAR.ltrendweak,
                           mk.NOAR.ltrendmed = mk.ts.NOAR.ltrendmed,
                           mk.NOAR.ltrendstrong = mk.ts.NOAR.ltrendstrong,
                           mk.NOAR.notrend = mk.ts.NOAR.notrend,
                           
                           mk.medAR.ltrendweak = mk.ts.medAR.ltrendweak,
                           mk.medAR.ltrendmed = mk.ts.medAR.ltrendmed,
                           mk.medAR.ltrendstrong = mk.ts.medAR.ltrendstrong,
                           mk.medAR.notrend = mk.ts.medAR.notrend,
                           
                           mk.strongAR.ltrendweak = mk.ts.strongAR.ltrendweak,
                           mk.strongAR.ltrendmed = mk.ts.strongAR.ltrendmed,
                           mk.strongAR.ltrendstrong = mk.ts.strongAR.ltrendstrong,
                           mk.strongAR.notrend = mk.ts.strongAR.notrend,
                           
                           pw.NOAR.ltrendweak = pw.ts.NOAR.ltrendweak,
                           pw.NOAR.ltrendmed = pw.ts.NOAR.ltrendmed,
                           pw.NOAR.ltrendstrong = pw.ts.NOAR.ltrendstrong,
                           pw.NOAR.notrend = pw.ts.NOAR.notrend,
                           
                           pw.medAR.ltrendweak = pw.ts.medAR.ltrendweak,
                           pw.medAR.ltrendmed = pw.ts.medAR.ltrendmed,
                           pw.medAR.ltrendstrong = pw.ts.medAR.ltrendstrong,
                           pw.medAR.notrend = pw.ts.medAR.notrend,
                           
                           pw.strongAR.ltrendweak = pw.ts.strongAR.ltrendweak,
                           pw.strongAR.ltrendmed = pw.ts.strongAR.ltrendmed,
                           pw.strongAR.ltrendstrong = pw.ts.strongAR.ltrendstrong,
                           pw.strongAR.notrend = pw.ts.strongAR.notrend,
                           
                           gls.NOAR.ltrendweak = gls.ts.NOAR.ltrendweak,
                           gls.NOAR.ltrendmed = gls.ts.NOAR.ltrendmed,
                           gls.NOAR.ltrendstrong = gls.ts.NOAR.ltrendstrong,
                           gls.NOAR.notrend = gls.ts.NOAR.notrend,
                           
                           gls.medAR.ltrendweak = gls.ts.medAR.ltrendweak,
                           gls.medAR.ltrendmed = gls.ts.medAR.ltrendmed,
                           gls.medAR.ltrendstrong = gls.ts.medAR.ltrendstrong,
                           gls.medAR.notrend = gls.ts.medAR.notrend,
                           
                           gls.strongAR.ltrendweak = gls.ts.strongAR.ltrendweak,
                           gls.strongAR.ltrendmed = gls.ts.strongAR.ltrendmed,
                           gls.strongAR.ltrendstrong = gls.ts.strongAR.ltrendstrong,
                           gls.strongAR.notrend = gls.ts.strongAR.notrend)
  mk.ts.NOAR.notrend <- NULL
  mk.ts.NOAR.ltrendweak <- NULL
  mk.ts.NOAR.ltrendmed <- NULL
  mk.ts.NOAR.ltrendstrong <- NULL
  
  mk.ts.medAR.notrend <- NULL
  mk.ts.medAR.ltrendweak <- NULL
  mk.ts.medAR.ltrendmed <- NULL
  mk.ts.medAR.ltrendstrong <- NULL
  
  mk.ts.strongAR.notrend <- NULL
  mk.ts.strongAR.ltrendweak <- NULL
  mk.ts.strongAR.ltrendmed <- NULL
  mk.ts.strongAR.ltrendstrong <- NULL
  
  pw.ts.NOAR.notrend <- NULL
  pw.ts.NOAR.ltrendweak <- NULL
  pw.ts.NOAR.ltrendmed <- NULL
  pw.ts.NOAR.ltrendstrong <- NULL
  
  pw.ts.medAR.notrend <- NULL
  pw.ts.medAR.ltrendweak <- NULL
  pw.ts.medAR.ltrendmed <- NULL
  pw.ts.medAR.ltrendstrong <- NULL
  
  pw.ts.strongAR.notrend <- NULL
  pw.ts.strongAR.ltrendweak <- NULL
  pw.ts.strongAR.ltrendmed <- NULL
  pw.ts.strongAR.ltrendstrong <- NULL
  
  gls.ts.NOAR.notrend <- NULL
  gls.ts.NOAR.ltrendweak <- NULL
  gls.ts.NOAR.ltrendmed <- NULL
  gls.ts.NOAR.ltrendstrong <- NULL
  
  gls.ts.medAR.notrend <- NULL
  gls.ts.medAR.ltrendweak <- NULL
  gls.ts.medAR.ltrendmed <- NULL
  gls.ts.medAR.ltrendstrong <- NULL
  
  gls.ts.strongAR.notrend <- NULL
  gls.ts.strongAR.ltrendweak <- NULL
  gls.ts.strongAR.ltrendmed <- NULL
  gls.ts.strongAR.ltrendstrong <- NULL

  
  names(sim_results) <- paste0(names(sim_results),'.',m)
  write.csv(sim_results, file = paste0("sim_results_",m,"_3_11.csv"))
  
}

proc.time() - ptm


