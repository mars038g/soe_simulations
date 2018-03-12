#Mann kendall with pre-whitening coverage intervals
library(boot);library(Kendall)
library(zoo);library(zyp)
library(trend);library(dplyr)
library(AICcmodavg)
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

n = 10 #number of simulations
ARsd <- .54^.5 #standard deviation


#set phi
NOAR <- list()
weakAR <- 0.1
medAR <- 0.433
strongAR <- 0.8

#place holders for generated data
series <- NULL
mk.coverage.NOAR.notrend <- NULL
mk.coverage.NOAR.ltrendweak <- NULL
mk.coverage.NOAR.ltrendmed <- NULL
mk.coverage.NOAR.ltrendstrong <- NULL

mk.coverage.medAR.notrend <- NULL
mk.coverage.medAR.ltrendweak <- NULL
mk.coverage.medAR.ltrendmed <- NULL
mk.coverage.medAR.ltrendstrong <- NULL

mk.coverage.strongAR.notrend <- NULL
mk.coverage.strongAR.ltrendweak <- NULL
mk.coverage.strongAR.ltrendmed <- NULL
mk.coverage.strongAR.ltrendstrong <- NULL

pw.coverage.NOAR.notrend <- NULL
pw.coverage.NOAR.ltrendweak <- NULL
pw.coverage.NOAR.ltrendmed <- NULL
pw.coverage.NOAR.ltrendstrong <- NULL

pw.coverage.medAR.notrend <- NULL
pw.coverage.medAR.ltrendweak <- NULL
pw.coverage.medAR.ltrendmed <- NULL
pw.coverage.medAR.ltrendstrong <- NULL

pw.coverage.strongAR.notrend <- NULL
pw.coverage.strongAR.ltrendweak <- NULL
pw.coverage.strongAR.ltrendmed <- NULL
pw.coverage.strongAR.ltrendstrong <- NULL

gls.coverage.NOAR.notrend <- NULL
gls.coverage.NOAR.ltrendweak <- NULL
gls.coverage.NOAR.ltrendmed <- NULL
gls.coverage.NOAR.ltrendstrong <- NULL

gls.coverage.medAR.notrend <- NULL
gls.coverage.medAR.ltrendweak <- NULL
gls.coverage.medAR.ltrendmed <- NULL
gls.coverage.medAR.ltrendstrong <- NULL

gls.coverage.strongAR.notrend <- NULL
gls.coverage.strongAR.ltrendweak <- NULL
gls.coverage.strongAR.ltrendmed <- NULL
gls.coverage.strongAR.ltrendstrong <- NULL

for (m in c(10,20,30)){
  notrend <- 0
  ltrendweak <- -0.262 + (0.004 * c(1:m)) 
  ltrendmed <- -0.262 + (0.051 * c(1:m)) 
  ltrendstrong <- -0.262 + (0.147 * c(1:m)) 
  
  for (k in c("notrend","ltrendweak","ltrendmed","ltrendstrong")){
    
    for (j in c("NOAR","medAR","strongAR")){
    
      # generate original data to test coverage intervals
      og <- arima.sim(list(ar = get(j)), n=m, rand.gen=rnorm, sd =ARsd)
      og <- get(k) + og
      og <- data.frame(series = og,
                       time = 1:length(og))
      lm_out_og <- fit_lm(dat = og)
      newtime <- seq(1, nrow(og), 1)
      newdata <- data.frame(time = newtime,
                            time2 = newtime^2)
      #gls
      if (length(lm_out_og$model$coefficients) == 2){
        gls_sim$model$coefficients[2] <- abs(lm_out_og$model$coefficients[2]) 
      } else {
        lm_out_og$model$coefficients[2] <- abs(lm_out_og$model$coefficients[2]) 
        lm_out_og$model$coefficients[3] <- abs(lm_out_og$model$coefficients[3]) 
      }
      gls_og <- AICcmodavg::predictSE(lm_out_og$model, 
                                          newdata = newdata,
                                          se.fit = T)
      gls_og <- gls_og$fit
  
      #mann kendall
      mk_og <- zyp.sen(series~time, og)
      mk_og <- mk_og$coefficients[1]+abs(og$time*mk_og$coefficients[2])
  
      #mann kendall with pre-whitening
      pw_og <- zyp.trend.vector(og$series, method = "yuepilon")
      pw_og <- pw_og[11]+abs(og$time*pw_og[2])
        
        for (i in 1:n){
          #generate simulations
          dat <- arima.sim(list(ar = get(j)), n=m, rand.gen=rnorm, sd =ARsd)
          
          #add autocorrelated error structure to trend
          dat <- get(k) + dat
          
          #what fraction of the confidence intervals contain the original trendline?
          dat <- data.frame(series = dat,
                            time = 1:length(dat))
          
          #simulate trends
          pw_sim <- zyp.trend.vector(dat$series,method='yuepilon')
          mk_sim <- MannKendall(dat$series)
          gls_sim <- try(fit_lm(dat), silent = T)
          
          
          if (length(gls_sim$model$coefficients) == 2){
            gls_sim$model$coefficients[2] <- abs(gls_sim$model$coefficients[2]) 
          } else {
            gls_sim$model$coefficients[2] <- abs(gls_sim$model$coefficients[2]) 
            gls_sim$model$coefficients[3] <- abs(gls_sim$model$coefficients[3]) 
          }

          
          #get model fits
          pw_pred <- pw_sim[11]+abs(dat$time*pw_sim[2])
          mk_pred <- zyp.sen(series~time, dat)
          mk_pred <- mk_pred$coefficients[1]+abs(dat$time*mk_pred$coefficients[2])
          gls_pred <- AICcmodavg::predictSE(gls_sim$model, 
                                           newdata = newdata,
                                           se.fit = T)
          gls_pred_se <- gls_pred$se.fit
          gls_pred <- gls_pred$fit
  
          ##get confidence intervals from fits
          #mann kendall
          tau_func <- function(z) MannKendall(z)$tau
          boot.out <- tsboot(dat$series, tau_func, R=500, l=5, sim="fixed")
          ci_upr_mk <- boot.ci(boot.out, type="norm")$normal[3]
          ci_upr_mk <- mk_pred + abs(ci_upr_mk)
          ci_lwr_mk <- boot.ci(boot.out, type="norm")$normal[2]
          ci_lwr_mk <- mk_pred - abs(ci_lwr_mk)  
          
          #MK-PW
          pw_tau_func <- function(z) zyp.trend.vector(z, method = "yuepilon")[5]
          boot.out <- tsboot(dat$series, pw_tau_func, R=500, l=5, sim="fixed")
          ci_upr_pw <- boot.ci(boot.out, type="norm")$normal[3]
          ci_upr_pw <- pw_pred + abs(ci_upr_pw)
          ci_lwr_pw <- boot.ci(boot.out, type="norm")$normal[2]
          ci_lwr_pw <- pw_pred - abs(ci_lwr_pw)  
          
          #gls
          ci_upr_gls <- gls_pred + 1.96*gls_pred_se
          ci_lwr_gls <- gls_pred - 1.96*gls_pred_se
  
          # collect the simulations inside (1) and outside (0) the coverage intervals
          if (all(mk_og < ci_upr_mk) & all(mk_og > ci_lwr_mk)){
            assign(paste0("mk.coverage.",j,".",k),rbind(get(paste0("mk.coverage.",j,".",k)),1))
          } else {
            assign(paste0("mk.coverage.",j,".",k),rbind(get(paste0("mk.coverage.", j,".",k)),0))
          }
          
          if (all(pw_og<ci_upr_pw) & all(pw_og > ci_lwr_pw)){
            assign(paste0("pw.coverage.",j,".",k),rbind(get(paste0("pw.coverage.",j,".",k)),1))
          } else {
            assign(paste0("pw.coverage.",j,".",k),rbind(get(paste0("pw.coverage.", j,".",k)),0))
          }
          
          if (all(gls_og<ci_upr_gls) & all(gls_og > ci_lwr_gls)){
            assign(paste0("gls.coverage.",j,".",k),rbind(get(paste0("gls.coverage.",j,".",k)),1))
          } else {
            assign(paste0("gls.coverage.",j,".",k),rbind(get(paste0("gls.coverage.", j,".",k)),0))
          }
    
          
      }
      
    } 
    
  }  
  
  sim_results = data.frame(mk.NOAR.ltrendweak = mk.coverage.NOAR.ltrendweak,
                            mk.NOAR.ltrendmed = mk.coverage.NOAR.ltrendmed,
                            mk.NOAR.ltrendstrong = mk.coverage.NOAR.ltrendstrong,
                            mk.NOAR.notrend = mk.coverage.NOAR.notrend,
                           
                            mk.medAR.ltrendweak = mk.coverage.medAR.ltrendweak,
                            mk.medAR.ltrendmed = mk.coverage.medAR.ltrendmed,
                            mk.medAR.ltrendstrong = mk.coverage.medAR.ltrendstrong,
                            mk.medAR.notrend = mk.coverage.medAR.notrend,
                           
                            mk.strongAR.ltrendweak = mk.coverage.strongAR.ltrendweak,
                            mk.strongAR.ltrendmed = mk.coverage.strongAR.ltrendmed,
                            mk.strongAR.ltrendstrong = mk.coverage.strongAR.ltrendstrong,
                           mk.strongAR.notrend = mk.coverage.strongAR.notrend,
                           
                           pw.NOAR.ltrendweak = pw.coverage.NOAR.ltrendweak,
                           pw.NOAR.ltrendmed = pw.coverage.NOAR.ltrendmed,
                           pw.NOAR.ltrendstrong = pw.coverage.NOAR.ltrendstrong,
                           pw.NOAR.notrend = pw.coverage.NOAR.notrend,
                           
                           pw.medAR.ltrendweak = pw.coverage.medAR.ltrendweak,
                           pw.medAR.ltrendmed = pw.coverage.medAR.ltrendmed,
                           pw.medAR.ltrendstrong = pw.coverage.medAR.ltrendstrong,
                           pw.medAR.notrend = pw.coverage.medAR.notrend,
                           
                           pw.strongAR.ltrendweak = pw.coverage.strongAR.ltrendweak,
                           pw.strongAR.ltrendmed = pw.coverage.strongAR.ltrendmed,
                           pw.strongAR.ltrendstrong = pw.coverage.strongAR.ltrendstrong,
                           pw.strongAR.notrend = pw.coverage.strongAR.notrend,
                           
                           gls.NOAR.ltrendweak = gls.coverage.NOAR.ltrendweak,
                           gls.NOAR.ltrendmed = gls.coverage.NOAR.ltrendmed,
                           gls.NOAR.ltrendstrong = gls.coverage.NOAR.ltrendstrong,
                           gls.NOAR.notrend = gls.coverage.NOAR.notrend,
                           
                           gls.medAR.ltrendweak = gls.coverage.medAR.ltrendweak,
                           gls.medAR.ltrendmed = gls.coverage.medAR.ltrendmed,
                           gls.medAR.ltrendstrong = gls.coverage.medAR.ltrendstrong,
                           gls.medAR.notrend = gls.coverage.medAR.notrend,
                           
                           gls.strongAR.ltrendweak = gls.coverage.strongAR.ltrendweak,
                           gls.strongAR.ltrendmed = gls.coverage.strongAR.ltrendmed,
                           gls.strongAR.ltrendstrong = gls.coverage.strongAR.ltrendstrong,
                           gls.strongAR.notrend = gls.coverage.strongAR.notrend)
  mk.coverage.NOAR.notrend <- NULL
  mk.coverage.NOAR.ltrendweak <- NULL
  mk.coverage.NOAR.ltrendmed <- NULL
  mk.coverage.NOAR.ltrendstrong <- NULL
  
  mk.coverage.medAR.notrend <- NULL
  mk.coverage.medAR.ltrendweak <- NULL
  mk.coverage.medAR.ltrendmed <- NULL
  mk.coverage.medAR.ltrendstrong <- NULL
  
  mk.coverage.strongAR.notrend <- NULL
  mk.coverage.strongAR.ltrendweak <- NULL
  mk.coverage.strongAR.ltrendmed <- NULL
  mk.coverage.strongAR.ltrendstrong <- NULL
  
  pw.coverage.NOAR.notrend <- NULL
  pw.coverage.NOAR.ltrendweak <- NULL
  pw.coverage.NOAR.ltrendmed <- NULL
  pw.coverage.NOAR.ltrendstrong <- NULL
  
  pw.coverage.medAR.notrend <- NULL
  pw.coverage.medAR.ltrendweak <- NULL
  pw.coverage.medAR.ltrendmed <- NULL
  pw.coverage.medAR.ltrendstrong <- NULL
  
  pw.coverage.strongAR.notrend <- NULL
  pw.coverage.strongAR.ltrendweak <- NULL
  pw.coverage.strongAR.ltrendmed <- NULL
  pw.coverage.strongAR.ltrendstrong <- NULL
  
  gls.coverage.NOAR.notrend <- NULL
  gls.coverage.NOAR.ltrendweak <- NULL
  gls.coverage.NOAR.ltrendmed <- NULL
  gls.coverage.NOAR.ltrendstrong <- NULL
  
  gls.coverage.medAR.notrend <- NULL
  gls.coverage.medAR.ltrendweak <- NULL
  gls.coverage.medAR.ltrendmed <- NULL
  gls.coverage.medAR.ltrendstrong <- NULL
  
  gls.coverage.strongAR.notrend <- NULL
  gls.coverage.strongAR.ltrendweak <- NULL
  gls.coverage.strongAR.ltrendmed <- NULL
  gls.coverage.strongAR.ltrendstrong <- NULL
  
  names(sim_results) <- paste0(names(sim_results),'.',m)
  write.csv(sim_results, file = paste0("sim_results_",m,"_3_11.csv"))
  rm(sim_results)

}
#coverage <- colSums(sim_results)/n





