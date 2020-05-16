##########
#
# Author: Samuel Rosin
# Date: May 7, 2020
# 
# This file reads in each different simulated datasets 
# and creates a single long-format dataset 
# with all true estimand values for later use. 
#

###### set these parameters 
n_i <- 50; m <- 500; n_sims <- 1000
user_home_directory <- "/nas/longleaf/home/srosin/" #home directory on a SLURM-managed computing cluster

setwd(paste(user_home_directory,"truthtables/truths_ni",n_i,"_m",m,sep=""))

alpha <- rep(NA,9*n_sims); mu0 <- rep(NA,9*n_sims)
mu1 <- rep(NA,9*n_sims); mu.marg <- rep(NA,9*n_sims)
truth.table <- data.frame(alpha,mu0,mu1,mu.marg)

i <- 1
#### read in truth files (files containing true estimand values for each dataset)
for(arg in 1:n_sims){
    filename <- paste("truth_",arg,".csv",sep="")
    row.ixs <- i:(i+8)
    truth.table[row.ixs,] <- read.csv(file=filename, header=T)
    i <- i + 9
}
truth.avgs.mu0 <- aggregate(formula = mu0 ~ alpha, data=truth.table, FUN=mean)
truth.avgs.mu1 <- aggregate(formula = mu1 ~ alpha, data=truth.table, FUN=mean)
truth.avgs <- merge(truth.avgs.mu0, truth.avgs.mu1, by="alpha")

setwd(paste(user_home_directory,"truthtables",sep=""))
out_filename <- paste("truthtable_ni",n_i,"_m",m,".csv",sep="")
write.csv(x = truth.avgs, file = out_filename)
