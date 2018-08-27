#testing coverage intervals with linear model

#number simulations
n <- 1000

#length series
m <- 50

coverage <- NULL

for (m in c(50)){
  #create original confidence intervals from random data
  og <- rnorm(m, sd = 1)
  x <- 1:length(og)
  og_mod <- lm(og ~ x)
  newdata = data.frame(x = m)
  
  ci_upr_lm <- predict(og_mod, newdata, interval="confidence")[3]
  ci_upr_lm <- (coef(og_mod)[1] + coef(og_mod)[2]*x) + abs(ci_upr_lm)
  ci_lwr_lm <- predict(og_mod, newdata, interval="confidence")[2]
  ci_lwr_lm <- (coef(og_mod)[1] + coef(og_mod)[2]*x) - abs(ci_lwr_lm)
  
  for (s in 1:n){
    #simulate random data n times
    sim <- rnorm(m, sd = 1)
    sim_mod <- lm(sim~x)
    sim_pred <- coef(sim_mod)[1] + coef(sim_mod)[2]*x
    
    #assign 1 to observation if it falls within the original confidence intervals,
    #and 0 otherwise
    for (g in 1:m){
      if (sim_pred[g] < ci_upr_lm[g] & sim_pred[g] > ci_lwr_lm[g]){
        assign('coverage',rbind(coverage, 1))
      } else {
        assign('coverage',rbind(coverage, 0))
     }
    }
  
  }  
}

#calculate coverage percentage
sum(coverage[,1])/nrow(coverage)
