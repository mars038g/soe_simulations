#coverage interval analyses
library(dplyr);library(tidyr)

rm(list = ls())
setwd('z:/shardison/soe/simulations')
d <- read.csv("sim_results_10_3_11.csv")

sim_summary <- function(csv, nsims, ts_length, VOI = NULL){
  #read in data
  m <- read.csv(csv)
  
  #add indexing variables
  m$rep <- rep(1:nsims, each = ts_length) #number of simulations
  m$index <- rep(1:ts_length,nsims) #length(time series)
  
  #convert from wide to long for summarizing 
  m_long <- gather(m, Var, Value,
                   paste0("mk.NOAR.ltrendweak.",ts_length):paste0("gls.strongAR.notrend.",ts_length),
                   factor_key=TRUE)
  
  #specify a variable of interest or get summaries from all tests
  if (is.null(VOI)){
    mk_sum <- m_long[grepl("mk",m_long$Var),]
    pw_sum <- m_long[grepl("pw",m_long$Var),]
    gls_sum <- m_long[grepl("gls",m_long$Var),]
    
    #find mean coverage of all simulations
    mk_mc <- mk_sum %>% group_by(Var,index) %>% dplyr::summarise(Total = sum(Value)) %>% 
      group_by(Var) %>% dplyr::summarise(Mean_coverage = mean(Total))
    pw_mc <- pw_sum %>% group_by(Var,index) %>% dplyr::summarise(Total = sum(Value)) %>% 
      group_by(Var) %>% dplyr::summarise(Mean_coverage = mean(Total))
    gls_mc <- gls_sum %>% group_by(Var,index) %>% dplyr::summarise(Total = sum(Value)) %>% 
      group_by(Var) %>% dplyr::summarise(Mean_coverage = mean(Total))
    out <- rbind(mk_mc,pw_mc,gls_mc)
    
    #Remove rows where GLS threw error
    out[-(out$Mean_coverage == 100),]
    out$Mean_coverage <- out$Mean_coverage/nsims
    return(out)
  } else {
    sum <- m_long[grepl(paste(VOI),m_long$Var),]
    sum1 <- sum %>% group_by(Var,index) %>% dplyr::summarise(Total = sum(Value)) %>% 
      group_by(Var) %>% dplyr::summarise(Mean_coverage = mean(Total))
    
    mk_mc <- NULL
    pw_mc <- NULL
    gls_mc <- NULL
    
    assign(paste0(VOI,"_mc"), rbind(get(paste0(VOI,"_mc")),sum1))
    
    out <- get(paste0(VOI,"_mc"))
    out[-(out$Mean_coverage == 100),]
    out$Mean_coverage <- out$Mean_coverage/nsims
    return(out)
  }
  
}

n_10 <- sim_summary(csv = "sim_results_10_3_11.csv",nsims = 500, ts_length = 10, VOI = NULL)
n_20 <- sim_summary(csv = "sim_results_20_3_11.csv",nsims = 50, ts_length = 20, VOI = NULL)
n_30 <- sim_summary(csv = "sim_results_30_3_11.csv",nsims = 500, ts_length = 30, VOI = NULL)

topleft <- NULL
topleft <- rbind(n_10[grepl("NOAR",n_10$Var) & grepl("ltrendstrong",n_10$Var),],
                 n_20[grepl("NOAR",n_20$Var) & grepl("ltrendstrong",n_20$Var),],
                 n_30[grepl("NOAR",n_30$Var) & grepl("ltrendstrong",n_30$Var),])
topleft$Var <- factor(topleft$Var)

topmid <- NULL
topmid <- rbind(n_10[grepl("medAR",n_10$Var) & grepl("ltrendstrong",n_10$Var),],
                n_20[grepl("medAR",n_20$Var) & grepl("ltrendstrong",n_20$Var),],
                n_30[grepl("medAR",n_30$Var) & grepl("ltrendstrong",n_30$Var),])
topmid$Var <- factor(topmid$Var)

topright <- NULL
topright <- rbind(n_10[grepl("strongAR",n_10$Var) & grepl("ltrendstrong",n_10$Var),],
                  n_20[grepl("strongAR",n_20$Var) & grepl("ltrendstrong",n_20$Var),],
                  n_30[grepl("strongAR",n_30$Var) & grepl("ltrendstrong",n_30$Var),])
topright$Var <- factor(topright$Var)

middle_2_left <- NULL
middle_2_left <- rbind(n_10[grepl("NOAR",n_10$Var) & grepl("ltrendmed",n_10$Var),],
                       n_20[grepl("NOAR",n_20$Var) & grepl("ltrendmed",n_20$Var),],
                       n_30[grepl("NOAR",n_30$Var) & grepl("ltrendmed",n_30$Var),])
middle_2_left$Var <- factor(middle_2_left$Var)

middle_2_mid <- NULL
middle_2_mid <- rbind(n_10[grepl("medAR",n_10$Var) & grepl("ltrendmed",n_10$Var),],
                      n_20[grepl("medAR",n_20$Var) & grepl("ltrendmed",n_20$Var),],
                      n_30[grepl("medAR",n_30$Var) & grepl("ltrendmed",n_30$Var),])
middle_2_mid$Var <- factor(middle_2_mid$Var)

middle_2_right <- NULL
middle_2_right <- rbind(n_10[grepl("strongAR",n_10$Var) & grepl("ltrendweak",n_10$Var),],
                        n_20[grepl("strongAR",n_20$Var) & grepl("ltrendweak",n_20$Var),],
                        n_30[grepl("strongAR",n_30$Var) & grepl("ltrendweak",n_30$Var),])
middle_2_right$Var <- factor(middle_2_right$Var)

middle_1_left <- NULL
middle_1_left <- rbind(n_10[grepl("NOAR",n_10$Var) & grepl("ltrendweak",n_10$Var),],
                       n_20[grepl("NOAR",n_20$Var) & grepl("ltrendweak",n_20$Var),],
                       n_30[grepl("NOAR",n_30$Var) & grepl("ltrendweak",n_30$Var),])
middle_1_left$Var <- factor(middle_1_left$Var)

middle_1_mid <- NULL
middle_1_mid <- rbind(n_10[grepl("medAR",n_10$Var) & grepl("ltrendweak",n_10$Var),],
                      n_20[grepl("medAR",n_20$Var) & grepl("ltrendweak",n_20$Var),],
                      n_30[grepl("medAR",n_30$Var) & grepl("ltrendweak",n_30$Var),])
middle_1_mid$Var <- factor(middle_1_mid$Var)


middle_1_right <- NULL
middle_1_right <- rbind(n_10[grepl("strongAR",n_10$Var) & grepl("ltrendweak",n_10$Var),],
                        n_20[grepl("strongAR",n_20$Var) & grepl("ltrendweak",n_20$Var),],
                        n_30[grepl("strongAR",n_30$Var) & grepl("ltrendweak",n_30$Var),])
middle_1_right$Var <- factor(middle_1_right$Var)


bottomleft <- NULL
bottomleft <- rbind(n_10[grepl("NOAR",n_10$Var) & grepl("notrend",n_10$Var),],
                    n_20[grepl("NOAR",n_20$Var) & grepl("notrend",n_20$Var),],
                    n_30[grepl("NOAR",n_30$Var) & grepl("notrend",n_30$Var),])
bottomleft$Var <- factor(bottomleft$Var)


bottommid <- NULL
bottommid <- rbind(n_10[grepl("medAR",n_10$Var) & grepl("notrend",n_10$Var),],
                   n_20[grepl("medAR",n_20$Var) & grepl("notrend",n_20$Var),],
                   n_30[grepl("medAR",n_30$Var) & grepl("notrend",n_30$Var),])
bottommid$Var <- factor(bottommid$Var)

bottomright <- NULL
bottomright <- rbind(n_10[grepl("strongAR",n_10$Var) & grepl("notrend",n_10$Var),],
                     n_20[grepl("strongAR",n_20$Var) & grepl("notrend",n_20$Var),],
                     n_30[grepl("strongAR",n_30$Var) & grepl("notrend",n_30$Var),])
bottomright$Var <- factor(bottomright$Var)


# plotting
axis.func <- function(lab = NULL,trend = NULL,
                      y.lab = TRUE, first = FALSE){
  if (y.lab == TRUE){
    mtext(2, text = "Percent Coverage", line = 1.5)
  }
  axis(1, at=text.loc,labels=FALSE)
  text(x=text.loc, y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
       labels=axis.lab, srt=45, adj=1, xpd=TRUE)
  axis(2, pos = 0)
  mtext(lab,cex = 1.5, line = 0.9)
  
  if(!is.null(trend)){
    text(x = par()$usr[2]+1.1, y = .5, trend,
         srt = 270, cex = 2.25)
  }
  if (first == T){
    text(x = c(1.6, 5.9, 10.1), y = -0.4, c("n = 10","20","30"))
  } else {
    text(x = c(1.6, 5.9, 10.1), y = -0.4, c("10","20","30"))
  }
  
}

par(mfrow = c(4,3), mar = c(4.3,3.5,2.5,3), xpd = T,
    mgp = c(2.4,.45,0),lty = 0)
axis.lab <- rep(c("MK","MK-PW","GLS"), 3)
text.loc <- c(0.5,1.6,2.7,4.8,5.9,7,9,10.1,11.2)
length.loc <- c(1.6,5.6,9.6)

########## top row #############
barplot(rep(1,9), col = "indianred", space = c(0.1,0.1,0.1,1,0.1,0.1,1),
        ylab = "",yaxt = "n")
barplot(topleft$Mean_coverage, space = c(0.1,0.1,0.1,1,0.1,0.1,1), col = "lightgreen",add = T,yaxt = "n")
axis.func(lab = "no AR", first = T)


barplot(rep(1,9), col = "indianred", space = c(0.1,0.1,0.1,1,0.1,0.1,1),
        ylab = "",yaxt = "n")
barplot(topmid$Mean_coverage, space = c(0.1,0.1,0.1,1,0.1,0.1,1), col = "lightgreen",add = T,yaxt = "n")
axis.func(lab = "med AR",y.lab = F)


barplot(rep(1,9), col = "indianred", space = c(0.1,0.1,0.1,1,0.1,0.1,1),
        ylab = "",yaxt = "n")
barplot(topright$Mean_coverage, space = c(0.1,0.1,0.1,1,0.1,0.1,1), col = "lightgreen",add = T,yaxt = "n")
axis.func(lab = "strong AR", trend = "strong trend",y.lab = F)


############# middle row ###############

barplot(rep(1,9), col = "indianred", space = c(0.1,0.1,0.1,1,0.1,0.1,1),
        ylab = "",yaxt = "n")
barplot(middle_2_left$Mean_coverage, space = c(0.1,0.1,0.1,1,0.1,0.1,1), col = "lightgreen",add = T,yaxt = "n")
axis.func(lab = "")

barplot(rep(1,9), col = "indianred", space = c(0.1,0.1,0.1,1,0.1,0.1,1),
        ylab = "",yaxt = "n")
barplot(middle_2_mid$Mean_coverage, space = c(0.1,0.1,0.1,1,0.1,0.1,1), col = "lightgreen",add = T,yaxt = "n")
axis.func(lab = "",y.lab = F)


barplot(rep(1,9), col = "indianred", space = c(0.1,0.1,0.1,1,0.1,0.1,1),
        ylab = "",yaxt = "n")
barplot(middle_2_right$Mean_coverage, space = c(0.1,0.1,0.1,1,0.1,0.1,1), col = "lightgreen",add = T,yaxt = "n")
axis.func(lab = "", trend = "med trend",y.lab = F)

############## middle row 1 ####################


barplot(rep(1,9), col = "indianred", space = c(0.1,0.1,0.1,1,0.1,0.1,1),
        ylab = "",yaxt = "n")
barplot(middle_1_left$Mean_coverage, space = c(0.1,0.1,0.1,1,0.1,0.1,1), col = "lightgreen",add = T,yaxt = "n")
axis.func(lab = "",y.lab = F)

barplot(rep(1,9), col = "indianred", space = c(0.1,0.1,0.1,1,0.1,0.1,1),
        ylab = "",yaxt = "n",y.lab = F)
barplot(middle_1_mid$Mean_coverage, space = c(0.1,0.1,0.1,1,0.1,0.1,1), col = "lightgreen",add = T,yaxt = "n")
axis.func(lab = "",y.lab = F)


barplot(rep(1,9), col = "indianred", space = c(0.1,0.1,0.1,1,0.1,0.1,1),
        ylab = "",yaxt = "n",y.lab = F)
barplot(middle_1_right$Mean_coverage, space = c(0.1,0.1,0.1,1,0.1,0.1,1), col = "lightgreen",add = T,yaxt = "n")
axis.func(lab = "", trend = "weak trend",y.lab = F)


############### bottom row ########################


barplot(rep(1,9), col = "indianred", space = c(0.1,0.1,0.1,1,0.1,0.1,1),
        ylab = "",yaxt = "n")
barplot(bottomleft$Mean_coverage, space = c(0.1,0.1,0.1,1,0.1,0.1,1), col = "lightgreen",add = T,yaxt = "n")
axis.func(lab = "")

barplot(rep(1,9), col = "indianred", space = c(0.1,0.1,0.1,1,0.1,0.1,1),
        ylab = "",yaxt = "n")
barplot(bottommid$Mean_coverage, space = c(0.1,0.1,0.1,1,0.1,0.1,1), col = "lightgreen",add = T,yaxt = "n")
axis.func(lab = "",y.lab = F)


barplot(rep(1,9), col = "indianred", space = c(0.1,0.1,0.1,1,0.1,0.1,1),
        ylab = "",yaxt = "n")
barplot(bottomright$Mean_coverage, space = c(0.1,0.1,0.1,1,0.1,0.1,1), col = "lightgreen",add = T,yaxt = "n")
axis.func(lab = "", trend = "no trend",y.lab = F)

