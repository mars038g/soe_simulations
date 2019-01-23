conf_mat <- function(df, test, filt){
  df <- df %>% filter(series.length == filt)
  
  #True positives
  true_pos_frac <- nrow(df[df$series.length != "no trend" &
                             df$p <= 0.05 &
                             df$test == test,])
  
  true_pos_tot <- nrow(df[df$trend != "no trend" &
                            df$test == test,])
  
  true_pos_freq <- true_pos_frac/true_pos_tot 
  
  #False positives
  
  false_pos_frac <- nrow(df[df$trend == "no trend" &
                              df$p <= 0.05 &
                              df$test == test,])
  
  false_pos_tot <- nrow(df[df$trend == "no trend" &
                             df$test == test,])
  
  false_pos_freq <- false_pos_frac/false_pos_tot
  
  #False negatives
  
  false_neg_frac <- nrow(df[df$trend != "no trend" &
                              df$p >= 0.05 &
                              df$test == test,])
  
  false_neg_tot <- nrow(df[df$trend != "no trend "&
                             df$test == test,])
  
  false_neg_freq <- false_neg_frac/false_neg_tot
  
  #true
  
  true_neg_frac <- nrow(df[df$trend == "no trend" &
                             df$p >= 0.05 &
                             df$test == test,])
  
  true_neg_tot <- nrow(df[df$trend == "no trend" &
                            df$test == test,])
  
  true_neg_freq <- true_neg_frac/true_neg_tot
  
  
  conf_mat <- data.frame(x = c("actual no","actual yes","actual no","actual yes"),
                         y = c('predicted no','predicted yes','predicted yes','predicted no'),
                         val = c(round(true_neg_freq,3), round(true_pos_freq,3),
                                 round(false_pos_freq,3),round(false_neg_freq,3)))
  
  return(conf_mat)
}
