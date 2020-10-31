#########
#
# Author: Samuel Rosin
# Date: October 23, 2020
#
# This script generates plots of bias and coverage
#
#

#### Plots where group size is n_i=10 and number of groups m varies
libs <- "/nas/longleaf/home/srosin/RLibs"
user_home_directory <- "/nas/longleaf/home/srosin/Chakladar_MS_Biometrics"

library(plyr)
library(dplyr)

number_of_groups <- c(10,50,100,200,300,400,500)
n_sims <- 1000

len_df <- length(number_of_groups)*2
bias <- rep(NA,len_df)
coverage.t <- rep(NA,len_df)
a <- c(rep(0,len_df/2),rep(1,len_df/2))
plot1_df <- data.frame(number_of_groups = rep(number_of_groups,2),
                       a, bias, coverage.t)

#### for each different value of m, read in bias and coverage results for alpha=0.5,
#### and store in plot1_df
for(num in number_of_groups){
  dir_path <- paste(user_home_directory,"/M",num,"_NI10",sep="")
  print(dir_path)
  setwd(dir_path)
  truthtable <- read.csv(paste("truthtable_","M",num,"_NI10",".csv",sep=""),header=T)
  
  results_dirpath <- paste(dir_path, "/results_datasets/",sep="")
  setwd(results_dirpath)
  print(results_dirpath)
  print(num)
  
  alpha <- rep(NA,9*n_sims); mu.hat0 <- rep(NA, 9*n_sims)
  mu.hat1 <- rep(NA,9*n_sims); ASE0 <- rep(NA,9*n_sims); ASE1 <- rep(NA, 9*n_sims)
  results.table <- data.frame(alpha,mu.hat0,mu.hat1,ASE0,ASE1)
  
  i <- 1
  for(arg in 1:n_sims){
    filename <- paste("results_M",num,"_NI10","_",arg,".csv",sep="")
    row.ixs <- i:(i+8)
    results.table[row.ixs,] <- read.csv(file=filename, header=T)
    i <- i + 9
  }
  
  truthtable.rel <- truthtable[,3:4]
  results.table <- cbind(results.table[,1],truthtable.rel,results.table[,2:5])
  colnames(results.table)[1] <- "alpha"
  results.mu0 <- results.table[,c(1,2,4,6)]
  results.mu1 <- results.table[,c(1,3,5,7)]
  
  # remove NA results, where propensity model did not fit 
  results.mu0 <- results.mu0[complete.cases(results.mu0),]
  results.mu1 <- results.mu1[complete.cases(results.mu1),]
  sims.worked <- nrow(results.mu0)/9
  
  #### mu0 results table
  results.mu0$bias <-  results.mu0$mu0 - results.mu0$mu.hat0
  results.mu0$lowerci <- results.mu0$mu.hat0 - qnorm(.975)*results.mu0$ASE0
  results.mu0$upperci <- results.mu0$mu.hat0 + qnorm(.975)*results.mu0$ASE0
  results.mu0$covers <-  ((results.mu0$lowerci <= results.mu0$mu0) &  
                            (results.mu0$mu0 <= results.mu0$upperci))
  results.mu0$lowerci.t <- results.mu0$mu.hat0 - qt(.975, num - 9)*results.mu0$ASE0
  results.mu0$upperci.t <- results.mu0$mu.hat0 + qt(.975, num - 9)*results.mu0$ASE0
  results.mu0$covers.t <- ((results.mu0$lowerci.t <= results.mu0$mu0) &  
                             (results.mu0$mu0 <= results.mu0$upperci.t))
  results.mu0.table <- results.mu0 %>% 
    dplyr::group_by(alpha) %>% 
    dplyr::summarise(mu0 = mean(mu0),
                     Bias = mean(bias),
                     ESE = sd(mu.hat0),
                     ASE = median(ASE0),
                     #mean.ASE = mean(ASE0),
                     EC = mean(covers)*100,
                     EC.t = mean(covers.t)*100,
                     .groups = "drop")
  
  # store bias and coverage for this number of groups and a=0
  plot1_df[(plot1_df$number_of_groups==num) & (a == 0),]$bias <- 
    results.mu0.table[results.mu0.table$alpha==0.5,]$Bias
  plot1_df[(plot1_df$number_of_groups==num) & (a == 0),]$coverage.t <- 
    results.mu0.table[results.mu0.table$alpha==0.5,]$EC.t
  
  ##### mu1 results table
  results.mu1$bias <- results.mu1$mu1 - results.mu1$mu.hat1
  results.mu1$lowerci <- results.mu1$mu.hat1 - qnorm(.975)*results.mu1$ASE1
  results.mu1$upperci <- results.mu1$mu.hat1 + qnorm(.975)*results.mu1$ASE1
  results.mu1$covers <-  ((results.mu1$lowerci <= results.mu1$mu1) &  
                            (results.mu1$mu1 <= results.mu1$upperci))
  results.mu1$lowerci.t <- results.mu1$mu.hat1 - qt(.975, num - 9)*results.mu1$ASE1
  results.mu1$upperci.t <- results.mu1$mu.hat1 + qt(.975, num - 9)*results.mu1$ASE1
  results.mu1$covers.t <- ((results.mu1$lowerci.t <= results.mu1$mu1) &  
                             (results.mu1$mu1 <= results.mu1$upperci.t))
  results.mu1.table <- results.mu1 %>% 
    dplyr::group_by(alpha) %>% 
    dplyr::summarise(mu1 = mean(mu1),
                     Bias = mean(bias),
                     ESE = sd(mu.hat1),
                     ASE = median(ASE1),
                     # mean.ASE = mean(ASE1),
                     EC = mean(covers)*100,
                     EC.t = mean(covers.t)*100,
                     .groups = "drop")
  
  # store bias and coverage for this number of groups and a=1
  plot1_df[(plot1_df$number_of_groups==num) & (a == 1),]$bias <- 
    results.mu1.table[results.mu0.table$alpha==0.5,]$Bias
  plot1_df[(plot1_df$number_of_groups==num) & (a == 1),]$coverage.t <- 
    results.mu1.table[results.mu0.table$alpha==0.5,]$EC.t
}
#View(plot1_df)

setwd(user_home_directory)
pdf("ni10_varying_number_of_groups_plots.pdf",width=7,height=4)
par(mfrow=c(1,2),
    mar=c(4, 5, 1, 1),
    cex.axis = .85,
    cex.lab = .85)
plot(plot1_df$number_of_groups, plot1_df$bias, type="n",
     xlab="Number of Groups",ylab="Bias",
     xaxt='n',yaxt='n',ylim=c(0,0.03))
axis(side=1,at=c(0,100,200,300,400,500))
axis(side=2,at=c(0.000,.010,.020,.030))
lines(plot1_df$number_of_groups[plot1_df$a==0],plot1_df$bias[plot1_df$a==0],col=1,type="l",lty=1)
lines(plot1_df$number_of_groups[plot1_df$a==1],plot1_df$bias[plot1_df$a==1],col=2,type="l",lty=2)
legend("topright",legend=c("a = 0","a = 1"),col=c("black","red"),lty=1:2,cex=.85)
abline(a=0,b=0,lty="dotted")


plot(x=plot1_df$number_of_groups, y=plot1_df$coverage.t, type="n",
     xlab="Number of Groups",ylab="Coverage (%)",
     xaxt='n',yaxt='n',ylim=c(0,100))
axis(side=1,at=c(0,100,200,300,400,500))
axis(side=2,at=c(0,20,40,60,80,100))
lines(plot1_df$number_of_groups[plot1_df$a==0],plot1_df$coverage.t[plot1_df$a==0],col=1,type="l",lty=1)
lines(plot1_df$number_of_groups[plot1_df$a==1],plot1_df$coverage.t[plot1_df$a==1],col=2,type="l",lty=2)
#legend("topright",legend=c("a=0","a=1"),col=c("black","red"),lty=1:2)
abline(a=95,b=0,lty="dotted")
dev.off()

############ Now repeat, but for m=500 groups and varying group size n_i
#group_size <- c(10,30,50,75,100,200)
group_size <- c(10,30,50,75,100,200)

len_df2 <- length(group_size)*2
bias <- rep(NA,len_df2/2)
coverage.t <- rep(NA,len_df2/2)
a <- c(rep(0,len_df2/2),rep(1,len_df2/2))
plot2_df <- data.frame(group_size = rep(group_size,2),
                       a, bias, coverage.t)

#### for each different value of m, read in bias and coverage results for alpha=0.5,
#### and store in plot1_df
for(size in group_size){
  dir_path <- paste(user_home_directory,"/M500_NI",size,sep="")
  print(dir_path)
  setwd(dir_path)
  truthtable <- read.csv(paste("truthtable_M500_NI",size,".csv",sep=""),header=T)
  
  results_dirpath <- paste(dir_path, "/results_datasets/",sep="")
  setwd(results_dirpath)
  print(results_dirpath)
  print(size)
  
  alpha <- rep(NA,9*n_sims); mu.hat0 <- rep(NA, 9*n_sims)
  mu.hat1 <- rep(NA,9*n_sims); ASE0 <- rep(NA,9*n_sims); ASE1 <- rep(NA, 9*n_sims)
  results.table <- data.frame(alpha,mu.hat0,mu.hat1,ASE0,ASE1)
  
  i <- 1
  for(arg in 1:n_sims){
    filename <- paste("results_M500_NI",size,"_",arg,".csv",sep="")
    row.ixs <- i:(i+8)
    results.table[row.ixs,] <- read.csv(file=filename, header=T)
    i <- i + 9
  }
  
  truthtable.rel <- truthtable[,3:4]
  results.table <- cbind(results.table[,1],truthtable.rel,results.table[,2:5])
  colnames(results.table)[1] <- "alpha"
  results.mu0 <- results.table[,c(1,2,4,6)]
  results.mu1 <- results.table[,c(1,3,5,7)]
  
  #### mu0 results table
  results.mu0$bias <-  results.mu0$mu0 - results.mu0$mu.hat0
  results.mu0$lowerci <- results.mu0$mu.hat0 - qnorm(.975)*results.mu0$ASE0
  results.mu0$upperci <- results.mu0$mu.hat0 + qnorm(.975)*results.mu0$ASE0
  results.mu0$covers <-  ((results.mu0$lowerci <= results.mu0$mu0) &  
                            (results.mu0$mu0 <= results.mu0$upperci))
  results.mu0$lowerci.t <- results.mu0$mu.hat0 - qt(.975, 500 - 9)*results.mu0$ASE0
  results.mu0$upperci.t <- results.mu0$mu.hat0 + qt(.975, 500 - 9)*results.mu0$ASE0
  results.mu0$covers.t <- ((results.mu0$lowerci.t <= results.mu0$mu0) &  
                             (results.mu0$mu0 <= results.mu0$upperci.t))
  results.mu0.table <- results.mu0 %>% 
    dplyr::group_by(alpha) %>% 
    dplyr::summarise(mu0 = mean(mu0),
                     Bias = mean(bias),
                     ESE = sd(mu.hat0),
                     ASE = median(ASE0),
                     #mean.ASE = mean(ASE0),
                     EC = mean(covers)*100,
                     EC.t = mean(covers.t)*100,
                     .groups = "drop")
  
  # store bias and coverage, under alpha=0.5, for this group size and a=0
  plot2_df[(plot2_df$group_size==size) & (a == 0),]$bias <- 
    results.mu0.table[results.mu0.table$alpha==0.5,]$Bias
  plot2_df[(plot2_df$group_size==size) & (a == 0),]$coverage.t <- 
    results.mu0.table[results.mu0.table$alpha==0.5,]$EC.t
  
  ##### mu1 results table
  results.mu1$bias <- results.mu1$mu1 - results.mu1$mu.hat1
  results.mu1$lowerci <- results.mu1$mu.hat1 - qnorm(.975)*results.mu1$ASE1
  results.mu1$upperci <- results.mu1$mu.hat1 + qnorm(.975)*results.mu1$ASE1
  results.mu1$covers <-  ((results.mu1$lowerci <= results.mu1$mu1) &  
                            (results.mu1$mu1 <= results.mu1$upperci))
  results.mu1$lowerci.t <- results.mu1$mu.hat1 - qt(.975, 500 - 9)*results.mu1$ASE1
  results.mu1$upperci.t <- results.mu1$mu.hat1 + qt(.975, 500 - 9)*results.mu1$ASE1
  results.mu1$covers.t <- ((results.mu1$lowerci.t <= results.mu1$mu1) &  
                             (results.mu1$mu1 <= results.mu1$upperci.t))
  results.mu1.table <- results.mu1 %>% 
    dplyr::group_by(alpha) %>% 
    dplyr::summarise(mu1 = mean(mu1),
                     Bias = mean(bias),
                     ESE = sd(mu.hat1),
                     ASE = median(ASE1),
                     # mean.ASE = mean(ASE1),
                     EC = mean(covers)*100,
                     EC.t = mean(covers.t)*100,
                     .groups = "drop")
  
  # store bias and coverage, under alpha=0.5 for this group size and a=1
  plot2_df[(plot2_df$group_size==size) & (a == 1),]$bias <- 
    results.mu1.table[results.mu0.table$alpha==0.5,]$Bias
  plot2_df[(plot2_df$group_size==size) & (a == 1),]$coverage.t <- 
    results.mu1.table[results.mu0.table$alpha==0.5,]$EC.t
}

setwd(user_home_directory)
pdf("m500_varying_group_sizes.pdf",width=7,height=4)

par(mfrow=c(1,2),
    mar=c(4, 5, 1, 1),
    cex.axis = .85,
    cex.lab = .85)
plot(plot2_df$group_size, plot2_df$bias, type="n",
     xlab="Group Size",ylab="Bias",
     xaxt='n',yaxt='n',
     ylim = c(-0.005, 0.01))
axis(side=1,at=c(0,50,100,150,200))
axis(side=2,at=c(-0.005,0,0.005,.01))
lines(plot2_df$group_size[plot2_df$a==0],plot2_df$bias[plot2_df$a==0],col=1,type="l",lty=1)
lines(plot2_df$group_size[plot2_df$a==1],plot2_df$bias[plot2_df$a==1],col=2,type="l",lty=2)
legend("topright",legend=c("a = 0","a = 1"),col=c("black","red"),lty=1:2,cex=.85)
abline(a=0,b=0,lty="dotted")

plot(x=plot2_df$group_size, y=plot2_df$coverage.t, type="n",
     xlab="Group Size",ylab="Coverage (%)",
     xaxt='n',yaxt='n',ylim=c(0,100))
axis(side=1,at=c(0,50,100,150,200,250))
axis(side=2,at=c(0,20,40,60,80,100))
lines(plot2_df$group_size[plot2_df$a==0],plot2_df$coverage.t[plot2_df$a==0],col=1,type="l",lty=1)
lines(plot2_df$group_size[plot2_df$a==1],plot2_df$coverage.t[plot2_df$a==1],col=2,type="l",lty=2)
#legend("topright",legend=c("a=0","a=1"),col=c("black","red"),lty=1:2)
abline(a=95,b=0,lty="dotted")
dev.off()
