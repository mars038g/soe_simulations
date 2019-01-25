fit_lm <- function(dat){
  
  dat <- dat %>% dplyr::filter(complete.cases(.))
  
  #Define null models --------------------------------------------------------------------
  
  #Null model with normal error structure
  constant_norm <-nlme::gls(series ~ 1, data = dat)
  
  #Null model with AR(1) error structure - cascades through optimizers
  constant_ar1 <- try(nlme::gls(series ~ 1,
                                data = dat,
                                correlation = corARMA(form = ~time, p = 1, q = 0)))
  if (class(constant_ar1) == "try-error"){
    
    constant_ar1 <- try(nlme::gls(series ~ 1,
                                  data = dat,
                                  control = glsControl(opt = "optim",
                                                       optimMethod = "Nelder-Mead"),
                                  correlation = corARMA(form = ~time, p = 1, q = 0)))
  } 
  
  
  #Null model with AR(2) error structure
  constant_ar2 <- try(nlme::gls(series ~ 1,
                                data = dat,
                                correlation = corARMA(form = ~time, p = 2, q = 0)))
  if (class(constant_ar2) == "try-error"){
    message("BFGS optimizer has failed, defaulting to Nelder-Mead routine (NULL AR2)")
    
    constant_ar2 <- try(nlme::gls(series ~ 1,
                                  data = dat,
                                  control = glsControl(opt = "optim",
                                                       optimMethod = "Nelder-Mead"),
                                  correlation = corARMA(form = ~time, p = 2, q = 0)))
    if (class(constant_ar2) == "try-error"){
      
      message("Nelder-Mead optimizer has failed, defaulting to SANN routine (NULL AR2)")
      
      constant_ar2 <- try(nlme::gls(series ~ 1,
                                    data = dat,
                                    control = glsControl(opt = "optim",
                                                         optimMethod = "SANN"),
                                    correlation = corARMA(form = ~time, p = 2, q = 0)))
      
      if (class(constant_ar2) == "try-error") message("NUll AR2 model has failed")
      
    }
  }
  
  #Fully parameterized models-------------------------------------------------------------
  
  # Linear model with normal error
  linear_norm <- nlme::gls(series ~ time, data = dat)
  
  # Linear model with AR1 error
  linear_ar1 <- try(nlme::gls(series ~ time,
                              data = dat,
                              correlation = corARMA(form = ~time, p = 1, q = 0)))
  if (class(linear_ar1) == "try-error"){
    linear_ar1 <- try(nlme::gls(series ~ time,
                                data = dat,
                                control = glsControl(opt = "optim",
                                                     optimMethod = "Nelder-Mead"),
                                correlation = corARMA(form = ~time, p = 1, q = 0)))
  }
  
  #Linear model with AR2 error
  linear_ar2 <- try(nlme::gls(series ~ time,
                              data = dat,
                              correlation = corARMA(form = ~time, p = 2, q = 0)))
  if (class(linear_ar2) == "try-error"){
    
    message("BFGS optimizer has failed, defaulting to Nelder-Mead routine (AR2)")
    
    linear_ar2 <- try(nlme::gls(series ~ time,
                                data = dat,
                                control = glsControl(opt = "optim",
                                                     optimMethod = "Nelder-Mead"),
                                correlation = corARMA(form = ~time, p = 2, q = 0)))
    if (class(linear_ar2) == "try-error"){
      
      message("Nelder-Mead optimizer has failed, defaulting to SANN routine (AR2)")
      
      linear_ar2 <- try(nlme::gls(series ~ time,
                                  data = dat,
                                  control = glsControl(opt = "optim",
                                                       optimMethod = "SANN"),
                                  correlation = corARMA(form = ~time, p = 2, q = 0)))
      
      if (class(linear_ar2) == "try-error") message("Linear AR2 model has failed")
    
    }
    
    
    
  }
  
  #Get phi (rho in manuscript) components of models with AR(p) error structure
  linear1_phi <- linear_ar1$modelStruct$corStruct
  linear1_phi <- coef(linear1_phi, unconstrained = FALSE)
  
  linear2_phi <- linear_ar2$modelStruct$corStruct
  linear2_phi <- coef(linear2_phi, unconstrained = FALSE)
  
  
  # Calculate AICs for all models
  df_aicc <-
    data.frame(model = c("linear_norm",
                         "linear_ar1",
                         "linear_ar2"),
               aicc  = c(AICc(linear_norm),
                         AICc(linear_ar1),
                         AICc(linear_ar2)),
               coefs = rbind(coef(linear_norm),
                             coef(linear_ar1),
                             coef(linear_ar2)),
               phi1 = c(NA,
                        linear1_phi,
                        linear2_phi[1]),
               phi2 = c(NA, 
                        NA,
                        linear2_phi[2]),
               # Calculate overall signifiance (need to use
               # ML not REML for this)
               pval = c(anova(update(constant_norm, method = "ML"),
                              update(linear_norm, method = "ML"))$`p-value`[2],
                        
                        anova(update(constant_ar1, method = "ML"),
                              update(linear_ar1, method = "ML"))$`p-value`[2],
                        
                        anova(update(constant_ar2, method = "ML"),
                              update(linear_ar2, method = "ML"))$`p-value`[2]))
  
  best_lm <-
    df_aicc %>%
    dplyr::filter(aicc == min(aicc)) 
  
  phi <- best_lm %>% dplyr::select(phi1,phi2)
  

  if (best_lm$model == "linear_norm") {
    model <- linear_norm
  } else if (best_lm$model == "linear_ar1") {
    model <- linear_ar1
  } else if (best_lm$model == "linear_ar2") {
    model <- linear_ar2
  } 
  
  
  return(list(best_lm = best_lm, 
              model = model))
}

fit_lm(dat)
