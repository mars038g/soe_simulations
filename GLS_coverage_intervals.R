#Build GLS model fit
library(dplyr);library(AICcmodavg)
library(nlme)

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

n = 200 #number of simulations
x = 30 #length of time series
ARsd <-  .54^.5 #standard deviation
ltrendstrong <- -0.262 + (0.051 * c(1:x)) #trend

#set phi
NOAR <- list()
weakAR <- 0.1
medAR <- 0.5
strongAR <- 0.9

set.seed(100)

#place holders for generated data
series <- NULL
inside.coverage.NOAR <- NULL
inside.coverage.weakAR <- NULL
inside.coverage.medAR <- NULL
inside.coverage.strongAR <- NULL

outside.coverage.NOAR <- NULL
outside.coverage.weakAR <- NULL
outside.coverage.medAR <- NULL
outside.coverage.strongAR <- NULL

gls.fit <- NULL

par(mar = c(4,4,3,1),mfrow = c(3,2))

for (j in c("NOAR","weakAR","medAR","strongAR")){
  # generate original data to test coverage intervals
  og <- as.numeric(arima.sim(list(ar = get(j)), n=x, rand.gen=rnorm, sd =ARsd))
  og <- ltrendstrong + og
  og <- data.frame(series = og,
                   time = 1:length(og))
  lm_out_og <- fit_lm(dat = og)
  newtime <- seq(1, nrow(og), 1)
  newdata <- data.frame(time = newtime,
                        time2 = newtime^2)
  
  lm_pred_og <- AICcmodavg::predictSE(lm_out_og$model, 
                                      newdata = newdata,
                                      se.fit = T)
  
  #Conf intervals
  #ci_upr <- lm_pred_og$fit + 1.96*lm_pred_og$se.fit
  #ci_lwr <- lm_pred_og$fit - 1.96*lm_pred_og$se.fit
  #Simulate data and plot 
  for (i in 1:n){
    
    #generate simulations
    dat <- as.numeric(arima.sim(list(ar = get(j)), n=x, rand.gen=rnorm, sd = .54^.5))
    
    #add autocorrelated error structure to trend
    dat <- ltrendstrong + dat

    #assign each time series to a row of a matrix
    assign("series",rbind(series,dat))
    
    #turn vector to data.frame for use in fit_lm()
    dat <- data.frame(series = dat,
                      time = 1:length(dat))
    lm_out <- fit_lm(dat = dat)
    
    #find confidence interval of GLS fit
    newtime <- seq(1, length(dat$series), 1)
    newdata <- data.frame(time = newtime,
                          time2 = newtime^2)
    
    lm_pred <- AICcmodavg::predictSE(lm_out$model, 
                                     newdata = newdata,
                                     se.fit = T)
    #Conf intervals
    ci_upr <- lm_pred$fit + 1.96*lm_pred$se.fit
    ci_lwr <- lm_pred$fit - 1.96*lm_pred$se.fit
    
    # #initialize null plot
    # if (i == 1){
    #   plot(NULL, ylim = c(min(lm_pred$fit)-max(lm_pred$fit)*.25,max(lm_pred$fit)+max(lm_pred$fit)*.25),
    #        xlim = c(1,x), las = 1, xlab = "Time", ylab = "Series")
    # } 
    
    # collect the samples inside and outside the coverage intervals
    if (all(lm_pred_og$fit<ci_upr) & all(lm_pred_og$fit > ci_lwr)){
      assign(paste0("inside.coverage.",j),rbind(get(paste0("inside.coverage.",j)),lm_pred$fit))
    } else {
      assign(paste0("outside.coverage.",j),rbind(get(paste0("outside.coverage.", j)),lm_pred$fit))
    }
    
    #get line fits
    assign("gls.fit",rbind(gls.fit,lm_pred$fit))

  }
  
  
  #plotting loop
  #for (i in 1:n){
  #  points(gls.fit[i,], col = i, type = "l", lwd = 3)  
  #}
  
  #reset gls.fit for each figure
  #gls.fit <- NULL
  #plot confidence interval
  #lines(newtime, lm_pred_og$fit, col = 1, lwd = 3
  
  #plot confidence interval and lines
  # polygon(c(rev(newtime), newtime), 
  #         c(rev(ci_upr), ci_lwr), 
  #         col = 'grey95', border = NA)
  # lines(newtime, ci_upr, 
  #       lty = 'dashed', col = 'black', lwd = .5)
  # lines(newtime, ci_lwr, 
  #       lty = 'dashed', col = 'black', lwd = .5)
  # lines(newtime, lm_pred_og$fit, col = 1, lwd = 3)
}

