#Building confusion matrices for summarizing simulation results
#Sean Hardison
#NEFSC - Integrated Statistics
#5/8/2018

library(dplyr)
library(ggplot2)
library(scales)
rm(list = ls())

setwd('z:/shardison/soe/simulations/sim_scripts_and_data')
load("P_results_n500.Rdata" )

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



