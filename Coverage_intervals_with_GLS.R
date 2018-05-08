#GLS Coverage Intervals 5/8/2018
#Sean Hardison

rm(list = ls())
#Mann kendall with pre-whitening coverage intervals
library(boot);library(Kendall)
library(zoo);library(zyp)
library(trend);library(dplyr)
library(AICcmodavg);library(nlme)
library(gtools)
ptm <- proc.time()

#Orginial GLS function written by Charles Perretti (http://github.com/perretti)
fit_lm <- function(dat, ar, m, ARsd, trend){
  success <- FALSE
  while(!success){
    dat <- dat %>% dplyr::filter(complete.cases(.))
    
    # Constant model (null model used to calculate 
    # overall p-value)
    constant_norm <-nlme::gls(series ~ 1, data = dat)
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
                                   pval = NA))}
    
    # Linear model with normal error
    linear_norm <- nlme::gls(series ~ time, data = dat)
    
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
                                   pval = NA))}
    linear_phi <- linear_ar1$modelStruct$corStruct
    linear_phi <-coef(linear_phi, unconstrained = FALSE)
    
    # Polynomial model with normal error
    dat$time2 <- dat$time^2
    poly_norm <- nlme::gls(series ~ time + time2, data = dat)
    
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
                                   pval = NA))}
    poly_phi <- poly_ar1$modelStruct$corStruct
    poly_phi <- coef(poly_phi, unconstrained = FALSE)
    
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
                 phi = c(0, 
                         poly_phi,
                         0,
                         linear_phi),
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

    phi <- best_lm$phi
    success <- (phi <= 0.8 & !invalid(phi))
    
    dat <- trend + dat
    dat <- arima.sim(list(ar = ar), n=m, rand.gen=rnorm, sd = ARsd)
    dat <- data.frame(series = dat,
                      time = 1:length(dat))}

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

set.seed(123)
n = 100 #number of simulations
ARsd <- .54^.5 #standard deviation

#set phi
NOAR <- list()
weakAR <- 0.1
medAR <- 0.433
strongAR <- 0.8

#placeholders for results
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

for (m in c(30)){
  
  notrend <- rep(0,m)
  ltrendweak <- -0.262 + (0.004 * c(1:m)) 
  ltrendmed <- -0.262 + (0.051 * c(1:m)) 
  ltrendstrong <- -0.262 + (0.147 * c(1:m)) 
  print(paste("m=",m))
  
  for (k in c("notrend","ltrendweak","ltrendmed","ltrendstrong")){
    
    for (j in c("strongAR","medAR","NOAR")){
      
      #Use true trend as base
      true_trend <- get(k)
      
      for (i in 1:n){

          #generate simulations
          dat <- arima.sim(list(ar = get(j)), n=m, rand.gen=rnorm, sd = ARsd)
          
          #add autocorrelated error structure to trend
          dat <- get(k) + dat
          
          #what fraction of the confidence intervals contain the original trendline?
          dat <- data.frame(series = dat,
                            time = 1:length(dat))
          
          #---------------------------------GLS---------------------------------#
          gls_sim <- tryCatch({
            
            newtime <- seq(1, m, 1)
            newdata <- data.frame(time = newtime,
                                  time2 = newtime^2)
            gls_sim <- fit_lm(dat = dat, ar = get(j), ARsd = ARsd, m = m, trend = get(k))
            
          }, 
          error = function(e) {
            gls_sim <- "error"
          })
          
          if (is.na(gls_sim[1]) | gls_sim == "error"){
            gls_pred <- rep(NA, m)
            ci_upr_gls <- rep(NA, m)
            ci_lwr_gls <- rep(NA, m)
          } else {
            gls_pred <- AICcmodavg::predictSE(gls_sim$model,
                                              newdata = newdata,
                                              se.fit = TRUE)
            
            ci_upr_gls <- gls_pred$fit + 1.96*gls_pred$se.fit
            ci_lwr_gls <- gls_pred$fit - 1.96*gls_pred$se.fit
         
            #plot(ci_upr_gls, ylim = c(-3,3))
            #points(ci_lwr_gls)
          }
          #---------------------------------Coverage---------------------------------#
          for (g in 1:m){
            if (is.na(ci_lwr_gls)[g] | is.na(ci_upr_gls)[g]){
              assign(paste0("gls.ts.",j,".",k),rbind(get(paste0("gls.ts.",j,".",k)),100))
            } else if ((true_trend[g] < ci_upr_gls[g]) & (true_trend[g] > ci_lwr_gls[g])){
              assign(paste0("gls.ts.",j,".",k),rbind(get(paste0("gls.ts.",j,".",k)),1))
            } else if ((true_trend[g] > ci_upr_gls[g]) | (true_trend[g] < ci_lwr_gls[g])){
              assign(paste0("gls.ts.",j,".",k),rbind(get(paste0("gls.ts.",j,".",k)),0))
            
            }
            
          }
      } 
      
    }
  }
  sim_results = data.frame(gls.NOAR.ltrendweak = gls.ts.NOAR.ltrendweak,
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
  
  names(sim_results) <- paste0(names(sim_results),'.',m)
  #write.csv(sim_results, file = paste0("sim_results_",m,"_5_1.csv"))
}


#-----------------------Visualization----------------------#
library(tidyr)
sims <- gather(sim_results, var, value, gls.NOAR.ltrendweak.30:gls.strongAR.notrend.30,
                 factor_key=TRUE)
processor <- function(ar, trend, n){
  
  assign(paste0("gls_",n),NULL)
  assign(paste0("gls_",n),rbind(get(paste0("gls_",n)),
                                sims[sims$var == paste0("gls.",ar,".",trend,".",n),]$value))
  gls <- sort(get(paste0("gls_",n)))
  out <- data.frame(gls = gls)
  return(out)
  
}


barplot_function <- function(ar = ar, trend = trend, trend.text = NULL,n = n,
                             product = "Percent Coverage", first  = T, lab = lab){
  #collect
  assign(paste0("mat_",n),NULL)
  mat <- processor(ar = ar, trend = trend,n = n)
  
  assign(paste0("mat_",n),mat) 
  
  #mean 
  gls_mean <- mean(get(paste0("mat_",n)))
  
  #vector for barplot
  vec <- c(gls_mean)
  loc <- c(1)
  bg <- 1
  
  #plot
  barplot(bg,  col = "red", ylab = "",xlab = "",yaxt = 'n',xaxt = 'n')
  barplot(vec,las = 1,add = TRUE,
          ylim = c(0,1), col = 'darkolivegreen2', border = NA,
          xaxt = "n")
  box()
  #text
  if(!is.null(trend.text)){
    text(x = par()$usr[2]+1, y = 1.5, trend.text,
         srt = 270, cex = 2.25, xpd = NA)
  }
  mtext(lab,cex = 1.5, line = 0.5)
  if (first == T){
    mtext(2, text = product, line = 2, cex = .8) 
    axis(1, at = c(1,4,7), labels = c("n = 10",20,30), tck = 0,cex.axis = 1.25)
  } else {
    axis(1, at = c(1,4,7), labels = c(10,20,30), tck = 0,cex.axis = 1.25)
  }
}



par(mfrow = c(4,3),mar = c(1.5,3.5,2.5,3), xpd = F,mgp=c(2.4,.45,0))
#par(mfrow = c(1,1))
#top row
barplot_function(ar = "NOAR",trend = "ltrendstrong", first = T, lab = "no AR")
barplot_function(ar = "medAR",trend = "ltrendstrong", first = F, lab = "med AR")
barplot_function(ar = "strongAR",trend = "ltrendstrong", first = F, lab = "strong AR",
                 trend.text = "strong trend")

barplot_function(ar = "NOAR",trend = "ltrendmed", first = T, lab = "")
barplot_function(ar = "medAR",trend = "ltrendmed", first = F, lab = "")
barplot_function(ar = "strongAR",trend = "ltrendmed", first = F, lab = "",
                 trend.text = "med trend")

barplot_function(ar = "NOAR",trend = "ltrendweak", first = T, lab = "")
barplot_function(ar = "medAR",trend = "ltrendweak", first = F, lab = "")
barplot_function(ar = "strongAR",trend = "ltrendweak", first = F, lab = "",
                 trend.text = "weak trend")

barplot_function(ar = "NOAR",trend = "ltrendweak", first = T, lab = "")
barplot_function(ar = "medAR",trend = "ltrendweak", first = F, lab = "")
barplot_function(ar = "strongAR",trend = "ltrendweak", first = F, lab = "",
                 trend.text = "no trend")


