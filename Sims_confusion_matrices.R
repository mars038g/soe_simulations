#Building confusion matrices for summarizing simulation results
#Sean Hardison
#NEFSC - Integrated Statistics
#5/8/2018

library(dplyr)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(colorspace)

rm(list = ls())

setwd('z:/shardison/soe/simulations/sim_scripts_and_data')
load("P_results_n500.Rdata" )
label <- function(variable,value){
  return(facet_names[value])
}


summary_func <- function(df, full_proc = T, type){
  agg <- df %>% group_by(`timeseries length`, `Trend strength`, `AR strength`,
                              method) %>%  
    dplyr::summarise(median = median(Value),
       mean = mean(Value),
       sd = sd(Value),
       n = n(),
       prop = length(Value[Value < 0.05])/n())
  
  if (full_proc & type == 2){
    agg <- agg %>% group_by(`timeseries length`, method) %>%
      filter(`Trend strength` != 'no trend') %>%
      dplyr::summarise(mean = mean(prop))
    return(agg)
  } else if (full_proc & type == 1){
    agg <- agg %>% group_by(`timeseries length`, method) %>%
      filter(`Trend strength` == 'no trend') %>%
    dplyr::summarise(mean = mean(prop))
    return(agg)
  }
    return(agg)
  }
 

false_pos_mat <- summary_func(p_results, type = 1)
false_rej_mat <- summary_func(p_results, type = 2)

#false rejections matrix
ggplot(data = false_rej_mat, aes(x = `timeseries length`, y = method)) +
  geom_tile(aes(fill = mean,size = 0.25), color = "grey") +
  scale_fill_gradient(low = "white", high = "steelblue") +

    geom_text(aes(x = `timeseries length`, y = method, label = round(mean,3), size = 1.4)) +
  
  theme(legend.position = "none",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.y = element_text(margin = margin(t = 0, r = -10, b = 0, l = 0)),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(margin = margin(t = -10, r = 0, b = 0, l = 0)),
        axis.ticks.x=element_blank()) +
  labs(title = expression(paste('Average rejection rate for each N (without ',beta,' = 0)'))) +
  coord_equal() 

#false positives matrix
ggplot(data = false_pos_mat, aes(x = `timeseries length`, y = method)) +
  geom_tile(aes(fill = mean,size = 0.25), color = "grey") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  
  geom_text(aes(x = `timeseries length`, y = method, label = round(mean,3), size = 1.4)) +
  
  theme(legend.position = "none",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.y = element_text(margin = margin(t = 0, r = -10, b = 0, l = 0)),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(margin = margin(t = -10, r = 0, b = 0, l = 0)),
        axis.ticks.x=element_blank()) +
  labs(title = expression(paste('Average rejection rate for each N (when ',beta,' = 0)'))) +
  coord_equal() 


#--------Confusion matrices - One per test----------#

#We want the percentage of true positives, false positives, true negatives, and false negatives

conf_mat <- function(df, test, ar = NULL){
  
  if (!is.null(ar)){
    ar <- ar
  } else{
    ar <- "AR"
  }
  #True positives
  true_pos_frac <- nrow(df[df$`Trend strength` != "no trend" &
                          df$Value <= 0.05 &
                          df$method == test &
                            grepl(ar,df$`AR strength`),])
  true_pos_tot <- nrow(df[df$`Trend strength` != "no trend" &
                            df$method == test&
                            grepl(ar,df$`AR strength`),])
  
  true_pos_freq <- true_pos_frac/true_pos_tot 
  
  #False positives
 
  false_pos_frac <- nrow(df[df$`Trend strength` == "no trend" &
                        df$Value <= 0.05 &
                        df$method == test&
                          grepl(ar,df$`AR strength`),])
  false_pos_tot <- nrow(df[df$`Trend strength` == "no trend" &
                             df$method == test&
                             grepl(ar,df$`AR strength`),])
  
  false_pos_freq <- false_pos_frac/false_pos_tot
  
  #False negatives

  false_neg_frac <- nrow(df[df$`Trend strength` != "no trend" &
                         df$Value >= 0.05 &
                         df$method == test&
                           grepl(ar,df$`AR strength`),])
  false_neg_tot <- nrow(df[df$`Trend strength` != "no trend "&
                             df$method == test&
                             grepl(ar,df$`AR strength`),])
  
  false_neg_freq <- false_neg_frac/false_neg_tot
  
  #true
  
  true_neg_frac <- nrow(df[df$`Trend strength` == "no trend" &
                         df$Value >= 0.05 &
                         df$method == test&
                           grepl(ar,df$`AR strength`),])
  true_neg_tot <- nrow(df[df$`Trend strength` == "no trend" &
                            df$method == test&
                            grepl(ar,df$`AR strength`),])
  
  true_neg_freq <- true_neg_frac/true_neg_tot

  
  conf_mat <- data.frame(x = c("actual no","actual yes","actual no","actual yes"),
                         y = c('predicted no','predicted yes','predicted yes','predicted no'),
                         val = c(round(true_neg_freq,3), round(true_pos_freq,3),
                                 round(false_pos_freq,3),round(false_neg_freq,3)))
  
  return(conf_mat)
}

mk <- cbind(conf_mat(p_results,test = "mk"), test = rep('mk',4))
pw <- cbind(conf_mat(p_results, test = "pw"), test = rep('pw',4))
gls <- cbind(conf_mat(p_results, test = "gls"), test = rep('gls',4))
fin <- rbind(mk, pw, gls)
fin$group <- factor(paste(fin$x,fin$y))


#Make matrices for white = good and orange = bad
fin_dif <- fin %>% group_by(group) %>%
  mutate(val, best_dif = ifelse(group == "actual no predicted no"|
                                  group == "actual yes predicted yes", (mean(val) - val)*-1, #best_dif is for assigning colors
                                mean(val) - val)) 

#For custom facet titles
facet_names <- list(
  'mk'="Mann-Kendall",
  'pw'="MK-TFPW",
  'gls'="GLS"
)

pal <- colorRampPalette(c("darkorange","white")) #color palette

#plot
ggplot(data = fin_dif, aes(x,y, fill = best_dif)) +
  facet_grid(. ~ test, labeller = label)+
  geom_tile(aes(size = 0.25))  +
  scale_fill_gradientn(colors = pal(10))+
  geom_text(aes(x = x, y = y, label = round(val,3), size = 1.4)) +
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_blank(),
        axis.text.y = element_text(margin = margin(t = 0, r = -10,
                                                   b = 0, l = 0),
                                   size = 15),
        axis.text.x = element_text(margin = margin(t = -10, r = 0,
                                                   b = 0, l = 0),
                                   size = 15),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = -0.1)) 

#This piece uses color to show divergence of cells from the mean of that cell class across models
fin_dif <- fin %>% group_by(group) %>%
  mutate(mean_dif = mean(val) - val)

fin_dif$group <- factor(fin_dif$group)

ggplot(data = fin_dif, aes(x,y, fill = mean_dif)) +
  facet_grid(. ~ test, labeller = label)+
  geom_tile(aes(size = 0.25))  +
  scale_fill_gradient2(low = "darkorange", mid = "white", high = "purple")+
  geom_text(aes(x = x, y = y, label = round(val,3), size = 1.4)) +
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_blank(),
        axis.text.y = element_text(margin = margin(t = 0, r = -10,
                                                   b = 0, l = 0),
                                   size = 15),
        axis.text.x = element_text(margin = margin(t = -10, r = 0,
                                                   b = 0, l = 0),
                                   size = 15),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = -0.1))


#Do the same for varying autocorrelation strengths - first strong AR class

mk <- cbind(conf_mat(p_results,test = "mk", ar = "strong"), test = rep('mk',4))
pw <- cbind(conf_mat(p_results, test = "pw", ar = "strong"), test = rep('pw',4))
gls <- cbind(conf_mat(p_results, test = "gls", ar = "strong"), test = rep('gls',4))


ar_str <- rbind(mk, pw, gls)
ar_str$group <- paste(ar_str$x,ar_str$y)

ar_str_dif <- ar_str %>% group_by(group) %>%
  mutate(mean_dif = mean(val) - val)


ar_str_dif$group <- factor(ar_str_dif$group)

ggplot(data = ar_str_dif, aes(x,y, fill = mean_dif)) +
  facet_grid(. ~ test, labeller = label)+
  geom_tile(aes(size = 0.25))  +
  scale_fill_gradient2(low = "darkorange", mid = "white", high = "purple")+
  geom_text(aes(x = x, y = y, label = round(val,3), size = 1.4)) +
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_blank(),
        axis.text.y = element_text(margin = margin(t = 0, r = -10,
                                                   b = 0, l = 0),
                                   size = 15),
        axis.text.x = element_text(margin = margin(t = -10, r = 0,
                                                   b = 0, l = 0),
                                   size = 15),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = -0.1))

#no AR class

mk <- cbind(conf_mat(p_results,test = "mk", ar = "no AR"), test = rep('mk',4))
pw <- cbind(conf_mat(p_results, test = "pw", ar = "no AR"), test = rep('pw',4))
gls <- cbind(conf_mat(p_results, test = "gls", ar = "no AR"), test = rep('gls',4))


ar_no <- rbind(mk, pw, gls)
ar_no$group <- paste(ar_no$x,ar_no$y)

ar_no_dif <- ar_no %>% group_by(group) %>%
  mutate(mean_dif = mean(val) - val)


ar_no_dif$group <- factor(ar_no_dif$group)

ggplot(data = ar_no_dif, aes(x,y, fill = mean_dif)) +
  facet_grid(. ~ test, labeller = label)+
  geom_tile(aes(size = 0.25))  +
  scale_fill_gradient2(low = "darkorange", mid = "white", high = "purple")+
  geom_text(aes(x = x, y = y, label = round(val,3), size = 1.4)) +
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_blank(),
        axis.text.y = element_text(margin = margin(t = 0, r = -10,
                                                   b = 0, l = 0),
                                   size = 15),
        axis.text.x = element_text(margin = margin(t = -10, r = 0,
                                                   b = 0, l = 0),
                                   size = 15),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = -0.1))
