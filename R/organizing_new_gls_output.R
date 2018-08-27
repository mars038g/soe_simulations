rm(list = ls())
library(dplyr)
library(tidyr)

d1 <- read.table('output.txt', header = T)

d1$method <- "gls.bs"

d_m <- d1 %>% rename(Trend.strength= beta, AR.strength = rho, timeseries.length  = nT, prop = pvalueBoot) %>%
  dplyr::select(-c(sigma, nT, pValueChi2)) %>% 
  mutate(Trend.strength, Trend.strength = plyr::mapvalues(Trend.strength, from = c(0,0.004, 0.051, 0.147),
                                                          to = c("no trend","weak trend", "medium trend", "strong trend"))) %>%
  mutate(AR.strength, AR.strength = plyr::mapvalues(AR.strength, from = c(0, 0.43, 0.80),
                                                          to = c("no AR", "medium AR", "strong AR")))
