#Simulating coverage intervals with Sen's slope
#Sean Hardison
#NEFSC - Integrated Statistics

rm(list = ls())

library(Kendall)
library(zyp)
library(trend)
library(nlme)

ptm <- proc.time()

set.seed(123)
n = 100 #number of simulations
ARsd <- .54^.5 #standard deviation

#set phi
NOAR <- list()
weakAR <- 0.1
medAR <- 0.433
strongAR <- 0.8

#bootstrap parameters
b_size <- 8  # size of moving blocks
nrep <- 599  #number of bootstrap replicates

#bootstrap functions
slope_func <- function(samp, resamp){
  zyp.trend.vector(samp[resamp])[2]
}
int_func <- function(samp, resamp){
  zyp.trend.vector(samp[resamp])[11]
}

#Plot simulation confidence intervals? 
fig = F

#place holders for generated data
mbb.ts.NOAR.ltrendweak <- NULL
mbb.ts.NOAR.notrend <- NULL
mbb.ts.NOAR.ltrendmed <- NULL
mbb.ts.NOAR.ltrendstrong <- NULL

mbb.ts.medAR.notrend <- NULL
mbb.ts.medAR.ltrendweak <- NULL
mbb.ts.medAR.ltrendmed <- NULL
mbb.ts.medAR.ltrendstrong <- NULL

mbb.ts.strongAR.notrend <- NULL
mbb.ts.strongAR.ltrendweak <- NULL
mbb.ts.strongAR.ltrendmed <- NULL
mbb.ts.strongAR.ltrendstrong <- NULL

perc.ts.NOAR.ltrendweak <- NULL
perc.ts.NOAR.notrend <- NULL
perc.ts.NOAR.ltrendmed <- NULL
perc.ts.NOAR.ltrendstrong <- NULL

perc.ts.medAR.notrend <- NULL
perc.ts.medAR.ltrendweak <- NULL
perc.ts.medAR.ltrendmed <- NULL
perc.ts.medAR.ltrendstrong <- NULL

perc.ts.strongAR.notrend <- NULL
perc.ts.strongAR.ltrendweak <- NULL
perc.ts.strongAR.ltrendmed <- NULL
perc.ts.strongAR.ltrendstrong <- NULL

for (m in c(10)){ #m = time series length
  
  notrend <- rep(0,m)
  ltrendweak <- -0.262 + (0.004 * c(1:m)) 
  ltrendmed <- -0.262 + (0.051 * c(1:m)) 
  ltrendstrong <- -0.262 + (0.147 * c(1:m)) 
  print(paste("m=",m))
  
  for (k in c("notrend","ltrendweak","ltrendmed","ltrendstrong")){
    
    for (j in c("strongAR","medAR","NOAR")){
      
      #Use true trend as base
      true_trend <- get(k)
      true_slope <- lm(true_trend ~ c(1:m))$coeff[2] 
      
      for (i in 1:n){
        #generate simulations
        dat <- arima.sim(list(ar = get(j)), n=m, rand.gen=rnorm, sd = ARsd)
        
        #add autocorrelated error structure to trend
        dat <- get(k) + dat
        
        #what fraction of the confidence intervals contain the original trendline?
        dat <- data.frame(series = dat,
                          time = 1:length(dat))
        
        #---------------------------------Sen's intervals--------------------------------#
        
        #Moving block bootstrap approach (mbb; Efron and Tibshirani)
        #Helpful code written by Andreas Buja (stat.wharton.upenn.edu/~buja/STAT-961/time-series-bootstrap.R)
        
        upper.int.bt <- rep(NA,nrep)  
        lower.int.bt <- rep(NA,nrep)         
        upper.beta.bt <- rep(NA,nrep)
        lower.beta.bt <- rep(NA,nrep)
        series <- dat$series #original simulation
        
        for(irep in 1:nrep) {                   
          series.bt <- rep(NA,m)               
          for(v in 1:ceiling(m/b_size)) {      #Resampling piece   
            endpoint <- sample(b_size:m, size=1)     
            series.bt[(v-1)*b_size+1:b_size] <- series[endpoint-(b_size:1)+1]
          }
          
          series.bt <- series.bt[1:m]           # trim overflow when b_size doesn't divide m
          blocked_dat <- data.frame(series = series.bt, #Add resampled data to data.frame
                            time = 1:m)
          slope <- zyp.sen(series~time, data = blocked_dat)
          
          #Find intervals for intercept and slope from resampled data
          upper.int.bt[irep] <- confint.zyp(slope)[3]
          lower.int.bt[irep] <- confint.zyp(slope)[1]
          upper.beta.bt[irep] <- confint.zyp(slope)[4]
          lower.beta.bt[irep] <- confint.zyp(slope)[2]

        }
        
        ci_upr_mbb <- mean(upper.beta.bt)*1:m + mean(upper.int.bt) 
        ci_lwr_mbb <- mean(lower.beta.bt)*1:m + mean(lower.int.bt) 
        
        #Percentile bootstrap (Perc; R.R. Wilcox 2010)
        
        int.boot <- boot(dat$series, int_func, R = nrep)
        lower.int.perc = boot.ci(int.boot, type = "perc")$percent[4]
        upper.int.perc = boot.ci(int.boot, type = "perc")$percent[5]
        
        slope.boot <- boot(dat$series, slope_func, R = nrep)
        lower.slope.perc = boot.ci(slope.boot, type = "perc")$percent[4]
        upper.slope.perc = boot.ci(slope.boot, type = "perc")$percent[5]
        
        ci_upr_perc <- mean(upper.slope.perc)*1:m + mean(upper.int.perc) 
        ci_lwr_perc <- mean(lower.slope.perc)*1:m + mean(lower.int.perc) 
        
        #---------------------------------Coverage---------------------------------#
        
        if (fig){
          plot(dat$series, ylim = c(-15,15))
          points(ci_upr_perc, col = 'red', type = "b")
          points(ci_lwr_perc, col = 'red', type = "b")
        
          points(ci_upr_mbb, type = "b")
          points(ci_lwr_mbb, type = "b")
        }
        
        for (g in 1:m){
          if ((true_trend[g] < ci_upr_mbb[g]) & (true_trend[g] > ci_lwr_mbb[g])){
            assign(paste0("mbb.ts.",j,".",k),rbind(get(paste0("mbb.ts.",j,".",k)),1))
          } else if ((true_trend[g] > ci_upr_mbb[g]) | (true_trend[g] < ci_lwr_mbb[g])){
            assign(paste0("mbb.ts.",j,".",k),rbind(get(paste0("mbb.ts.",j,".",k)),0))
          }
          
          if ((true_trend[g] < ci_upr_perc[g]) & (true_trend[g] > ci_lwr_perc[g])){
            assign(paste0("perc.ts.",j,".",k),rbind(get(paste0("perc.ts.",j,".",k)),1))
          } else if ((true_trend[g] > ci_upr_perc[g]) | (true_trend[g] < ci_lwr_perc[g])){
            assign(paste0("perc.ts.",j,".",k),rbind(get(paste0("perc.ts.",j,".",k)),0))
          }
          
        }
        
      } 
      
    }
  }
  sim_results = data.frame(mbb.NOAR.ltrendweak = mbb.ts.NOAR.ltrendweak,
                           mbb.NOAR.ltrendmed = mbb.ts.NOAR.ltrendmed,
                           mbb.NOAR.ltrendstrong = mbb.ts.NOAR.ltrendstrong,
                           mbb.NOAR.notrend = mbb.ts.NOAR.notrend,
                           
                           mbb.medAR.ltrendweak = mbb.ts.medAR.ltrendweak,
                           mbb.medAR.ltrendmed = mbb.ts.medAR.ltrendmed,
                           mbb.medAR.ltrendstrong = mbb.ts.medAR.ltrendstrong,
                           mbb.medAR.notrend = mbb.ts.medAR.notrend,
                           
                           mbb.strongAR.ltrendweak = mbb.ts.strongAR.ltrendweak,
                           mbb.strongAR.ltrendmed = mbb.ts.strongAR.ltrendmed,
                           mbb.strongAR.ltrendstrong = mbb.ts.strongAR.ltrendstrong,
                           mbb.strongAR.notrend = mbb.ts.strongAR.notrend,
                           
                           perc.NOAR.ltrendweak = perc.ts.NOAR.ltrendweak,
                           perc.NOAR.ltrendmed = perc.ts.NOAR.ltrendmed,
                           perc.NOAR.ltrendstrong = perc.ts.NOAR.ltrendstrong,
                           perc.NOAR.notrend = perc.ts.NOAR.notrend,
                           
                           perc.medAR.ltrendweak = perc.ts.medAR.ltrendweak,
                           perc.medAR.ltrendmed = perc.ts.medAR.ltrendmed,
                           perc.medAR.ltrendstrong = perc.ts.medAR.ltrendstrong,
                           perc.medAR.notrend = perc.ts.medAR.notrend,
                           
                           perc.strongAR.ltrendweak = perc.ts.strongAR.ltrendweak,
                           perc.strongAR.ltrendmed = perc.ts.strongAR.ltrendmed,
                           perc.strongAR.ltrendstrong = perc.ts.strongAR.ltrendstrong,
                           perc.strongAR.notrend = perc.ts.strongAR.notrend)
  
  names(sim_results) <- paste0(names(sim_results),'.',m)
  #write.csv(sim_results, file = paste0("sim_results_",m,"_5_1.csv"))
}

final_coverage <- colMeans(sim_results)

proc.time() - ptm
