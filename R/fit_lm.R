fit_lm <- function(dat, spec = FALSE, p = p){
  dat <- dat %>% dplyr::filter(complete.cases(.))
  
  # Constant model (null model used to calculate 
  # overall p-value)
  constant_norm <-nlme::gls(series ~ 1, data = dat)
  
  #spec parameter specifies whether arima.sim incorporated AR(1) error process. When there is no 
  #AR error in the time series, we switch the GLS models to all rely on a normal error generating process.
  if (!spec){
    constant_ar <-
      try(nlme::gls(series ~ 1,
                    data = dat,
                    correlation = corARMA(form = ~1|time, p = p, q = 0)))
    if (class(constant_ar) == "try-error"){
      return(best_lm <- data.frame(model = NA,
                                   aicc  = NA,
                                   coefs..Intercept = NA,
                                   coefs.time = NA,
                                   coefs.time2 = NA,
                                   pval = NA))
    }
    
  } else {
    constant_ar <-
      try(nlme::gls(series ~ 1,
                    data = dat))
    if (class(constant_ar) == "try-error"){
      return(best_lm <- data.frame(model = NA,
                                   aicc  = NA,
                                   coefs..Intercept = NA,
                                   coefs.time = NA,
                                   coefs.time2 = NA,
                                   pval = NA))
    }
  }
  # Linear model with normal error
  linear_norm <- nlme::gls(series ~ time, data = dat)
  
  # Linear model with AR1 error
  if (!spec){
    linear_ar <- 
      try(nlme::gls(series ~ time, 
                    data = dat,
                    correlation = corARMA(form = ~1|time, p = p, q = 0)))
    if (class(linear_ar) == "try-error"){
      return(best_lm <- data.frame(model = NA,
                                   aicc  = NA,
                                   coefs..Intercept = NA,
                                   coefs.time = NA,
                                   coefs.time2 = NA,
                                   pval = NA))
    }
  } else {
    linear_ar <- 
      try(nlme::gls(series ~ time, 
                    data = dat))
    if (class(linear_ar) == "try-error"){
      return(best_lm <- data.frame(model = NA,
                                   aicc  = NA,
                                   coefs..Intercept = NA,
                                   coefs.time = NA,
                                   coefs.time2 = NA,
                                   pval = NA))
    }
  }
  linear_phi <- linear_ar$modelStruct$corStruct
  linear_phi <-coef(linear_phi, unconstrained = FALSE)
  
  # Polynomial model with normal error
  dat$time2 <- dat$time^2
  poly_norm <- nlme::gls(series ~ time + time2, data = dat)
  
  # Polynomial model with AR1 error
  if (!spec){
    poly_ar <-
      try(nlme::gls(series ~ time + time2,
                    data = dat,
                    correlation = corARMA(form = ~1|time, p = p, q = 0)))
    if (class(poly_ar) == "try-error"){
      return(best_lm <- data.frame(model = NA,
                                   aicc  = NA,
                                   coefs..Intercept = NA,
                                   coefs.time = NA,
                                   coefs.time2 = NA,
                                   pval = NA))
      
    }
  }else {
    poly_ar <-
      try(nlme::gls(series ~ time + time2,
                    data = dat))
    if (class(poly_ar) == "try-error"){
      return(best_lm <- data.frame(model = NA,
                                   aicc  = NA,
                                   coefs..Intercept = NA,
                                   coefs.time = NA,
                                   coefs.time2 = NA,
                                   pval = NA))
      
    }
  }
  poly_phi <- poly_ar$modelStruct$corStruct
  poly_phi <- coef(poly_phi, unconstrained = FALSE)
 
  # Calculate AICs for all models
  df_aicc <-
    data.frame(model = c("poly_norm",
                         "poly_ar",
                         "linear_norm",
                         "linear_ar"),
               aicc  = c(AICc(poly_norm),
                         AICc(poly_ar),
                         AICc(linear_norm),
                         AICc(linear_ar)),
               coefs = rbind(coef(poly_norm),
                             coef(poly_ar),
                             c(coef(linear_norm), NA),
                             c(coef(linear_ar),  NA)),
               phi1 = c(NA, 
                       poly_phi[1],
                       NA,
                       linear_phi[1]),
               phi2 = c(NA, 
                        poly_phi[2],
                        NA,
                        linear_phi[2]),
               # Calculate overall signifiance (need to use
               # ML not REML for this)
               pval = c(anova(update(constant_norm, method = "ML"),
                              update(poly_norm, method = "ML"))$`p-value`[2],
                        anova(update(constant_ar, method = "ML"),
                              update(poly_ar, method = "ML"))$`p-value`[2],
                        anova(update(constant_norm, method = "ML"),
                              update(linear_norm, method = "ML"))$`p-value`[2],
                        anova(update(constant_ar, method = "ML"),
                              update(linear_ar, method = "ML"))$`p-value`[2]))
  
  best_lm <-
    df_aicc %>%
    dplyr::filter(aicc == min(aicc)) 

  phi <- best_lm %>% dplyr::select(phi1,phi2)
  
  if (best_lm$model == "poly_norm") {
    model <- poly_norm
  } else if (best_lm$model == "poly_ar") {
    model <- poly_ar
  } else if (best_lm$model == "linear_norm") {
    model <- linear_norm
  } else if (best_lm$model == "linear_ar") {
    model <- linear_ar
  }
  
  
  return(list(best_lm = best_lm, 
              model = model))
}
