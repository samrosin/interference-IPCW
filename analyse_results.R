##########
#
# Author: Samuel Rosin
# Revised Date: October 5, 2020
#
# This file reads in the different results files and outputs a file 
# output_tables. This .Rout file contains two results tables, 
# one each for a=0 and a=1, which can be copied and pasted into LaTeX.
# 
# The group size n_i, number of groups m, and name of home directory 
# must be specified in the shell script
#

#### set the following directory
user_home_directory <- "/nas/longleaf/home/srosin/Chakladar_MS_Biometrics/"
libs <- "/nas/longleaf/home/srosin/RLibs"

## collect arguments passed in from SLURM job script
args <- commandArgs(trailingOnly=TRUE)
ran_seed <- as.numeric(args[1])
m <- as.numeric(args[2])
n_i <- as.numeric(args[3])
n_sims <- as.numeric(args[4])
subdir_name <- as.character(args[5])

library(dplyr)
library(xtable)

dir_path <- paste(user_home_directory,subdir_name,"/",sep="")
print(dir_path)
setwd(dir_path)

truthtable <- read.csv(paste("truthtable_",subdir_name,".csv",sep=""),header=T)

results_dirpath <- paste(dir_path, "results_datasets/",sep="")
setwd(results_dirpath)


alpha <- rep(NA,9*n_sims); mu.hat0 <- rep(NA, 9*n_sims)
mu.hat1 <- rep(NA,9*n_sims); ASE0 <- rep(NA,9*n_sims); ASE1 <- rep(NA, 9*n_sims)
results.table <- data.frame(alpha,mu.hat0,mu.hat1,ASE0,ASE1)

i <- 1
for(arg in 1:n_sims){
  filename <- paste("results_",subdir_name,"_",arg,".csv",sep="")
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
results.mu0$lowerci.t <- results.mu0$mu.hat0 - qt(.975, m - 9)*results.mu0$ASE0
results.mu0$upperci.t <- results.mu0$mu.hat0 + qt(.975, m - 9)*results.mu0$ASE0
results.mu0$covers.t <- ((results.mu0$lowerci.t <= results.mu0$mu0) &  
                           (results.mu0$mu0 <= results.mu0$upperci.t))
results.mu0.table <- results.mu0 %>% 
  dplyr::group_by(alpha) %>% 
  dplyr::summarise(mu0 = mean(mu0),
                   Bias = mean(bias),
                   ESE = sd(mu.hat0),
                   ASE = median(ASE0),
                   #mean.ASE = mean(ASE0),
                   EC = mean(covers),
                   EC.t = mean(covers.t))
results.mu0.table$EC <- paste0(round(results.mu0.table$EC*100),"%")
results.mu0.table$EC.t <- paste0(round(results.mu0.table$EC.t*100),"%")

##### mu1 results table
results.mu1$bias <- results.mu1$mu1 - results.mu1$mu.hat1
results.mu1$lowerci <- results.mu1$mu.hat1 - qnorm(.975)*results.mu1$ASE1
results.mu1$upperci <- results.mu1$mu.hat1 + qnorm(.975)*results.mu1$ASE1
results.mu1$covers <-  ((results.mu1$lowerci <= results.mu1$mu1) &  
                          (results.mu1$mu1 <= results.mu1$upperci))
results.mu1$lowerci.t <- results.mu1$mu.hat1 - qt(.975, m - 9)*results.mu1$ASE1
results.mu1$upperci.t <- results.mu1$mu.hat1 + qt(.975, m - 9)*results.mu1$ASE1
results.mu1$covers.t <- ((results.mu1$lowerci.t <= results.mu1$mu1) &  
                           (results.mu1$mu1 <= results.mu1$upperci.t))
results.mu1.table <- results.mu1 %>% 
  dplyr::group_by(alpha) %>% 
  dplyr::summarise(mu1 = mean(mu1),
                   Bias = mean(bias),
                   ESE = sd(mu.hat1),
                   ASE = median(ASE1),
                   # mean.ASE = mean(ASE1),
                   EC = mean(covers),
                   EC.t = mean(covers.t))
results.mu1.table$EC <- paste0(round(results.mu1.table$EC*100),"%")
results.mu1.table$EC.t <- paste0(round(results.mu1.table$EC.t*100),"%")

###### print out the results to a table 
setwd(dir_path)
outfile_name <- paste("output_tables_",subdir_name,".Rout",sep="")

sink(file = outfile_name, append = TRUE, type = "output")
print(paste("Number of sims that completed:", sims.worked, "out of", n_sims))
print(xtable(results.mu0.table,digits=c(1,1,2,2,2,2,0,0)), 
      include.rownames=FALSE)
print(xtable(results.mu1.table,digits=c(1,1,2,2,2,2,0,0)), 
      include.rownames=FALSE)
sink()


