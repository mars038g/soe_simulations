conf_mat <- function(df, test, ar = NULL){
  
  df <- df %>% filter(series.length == 30)
  if (!is.null(ar)){
    ar <- ar
  } else{
    ar <- "AR"
  }
  #True positives
  true_pos_frac <- df %>% filter(trend != "No trend", p <= 0.05, Method == test) %>% count()
  true_pos_tot <- df %>% filter(trend != "No trend", Method == test) %>% count()
  true_pos_freq <- as.numeric(true_pos_frac/true_pos_tot)
  
  #False positives
  false_pos_frac <- df %>% filter(trend == "No trend", p <= 0.05, Method == test) %>% count()
  false_pos_tot <- df %>% filter(trend == "No trend", Method == test) %>% count()
  false_pos_freq <- as.numeric(false_pos_frac/false_pos_tot)
  
  #False negatives
  false_neg_frac <- df %>% filter(trend != "No trend", p >= 0.05, Method == test) %>% count()
  false_neg_tot <- df %>% filter(trend != "No trend", Method == test) %>% count()
  false_neg_freq <- as.numeric(false_neg_frac/false_neg_tot)
  
  #True negatives
  true_neg_frac <- df %>% filter(trend == "No trend", p >= 0.05, Method == test) %>% count()
  true_neg_tot <- df %>% filter(trend == "No trend", Method == test) %>% count()
  true_neg_freq <- as.numeric(true_neg_frac/true_neg_tot)
  
  conf_mat <- data.frame(x = c("actual no","actual yes","actual no","actual yes"),
                         y = c('predicted no','predicted yes','predicted yes','predicted no'),
                         val = c(round(true_neg_freq,3), round(true_pos_freq,3),
                                 round(false_pos_freq,3),round(false_neg_freq,3)))
  return(conf_mat)
}