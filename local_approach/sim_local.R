##########
#
# Author: Samuel Rosin
# Date: May 16, 2020
#
# Simulation study for Section 3 of the manuscript 
# "Inverse Probability Weighted Estimators of Vaccine Effects
# Accommodating Partial Interference and Censoring", by Chakladar et al.,
# available at https://arxiv.org/abs/1910.03536
#
# Be sure to read the README.md file for further instructions. 
#
# Need to open and fully run the following files to have access to necessary functions:
# (1) sim_dataset.R simulates the base datasets
# (2) calc_estimands.R calculates the true values of the estimands
# (3) ase_est.R estimates the asymptotic SEs (ASEs) and is called by (4)
# (4) ipcw_eval.R estimates the IPCW point estimates and ASEs
# 
# Some of this code in files (1) and (2) is adapted from Dr. Bradley Saul's 
# interferenceSim package available at https://github.com/bsaul/interferenceSim
##########


##### here is where the group size, number of groups, 
##### and number of simulations are set 
n_i <- 10; m <- 10; n_sims <- 1

time.until <- 100 #the time until event that we are interested in. 100 by default.
set.seed(69)
alphas <- c(.1,.2,.3,.4,.5,.6,.7,.8,.9) #Bernoulli allocation levels

library(parfm) #parametric frailty models
library(lme4) #mixed models
library(dplyr) #data manipulation
library(plyr)
library(geex) #ASE 
library(xtable) #outputting tables to LaTeX

#### define inv.logit function
inv.logit <- function(x){exp(x)/(1+exp(x))}

#### simulate datasets
time1 <- Sys.time()
sim.datasets <- replicate(n=n_sims,sim_dataset(n=n_i,m=m))
base.datasets <- sim.datasets[1,]; potential.outcomes.list <- sim.datasets[2,]
base.datasets <- mapply(FUN=cbind, base.datasets, sim.number=as.list(1:n_sims),SIMPLIFY=FALSE) #add sim number
time2 <- Sys.time(); time2 - time1; print(paste(n_sims,"datasets simulated"))

##### print the truth
time1 <- Sys.time()
truth.list <- lapply(X=potential.outcomes.list, FUN=calc_estimands,
                     alphas=alphas, time.until=time.until)
truth.df.long <- do.call(what=rbind,args=truth.list)
truth.avgs.mu0 <- aggregate(formula = mu0 ~ alpha, data=truth.df.long, FUN=mean)
truth.avgs.mu1 <- aggregate(formula = mu1 ~ alpha, data=truth.df.long, FUN=mean)
truth.avgs <- merge(truth.avgs.mu0, truth.avgs.mu1, by="alpha")
time2 <- Sys.time(); time2 - time1; print("truth list below")
truth.avgs

##### evaluate the IPCW estimates
time1 <- Sys.time() ### can use the timing to estimate how long, e.g., 1000 sims might take
ipcw.list <- lapply(X=base.datasets, FUN=ipcw.eval.allalphas, alphas=alphas, time.until = time.until)
time2 <- Sys.time(); print(time2-time1); print("finished estimating all datasets")

##### output the results in tables 
results.table <- do.call(what = rbind, args = ipcw.list)
results.table <- cbind(truth.df.long[,1],results.table[,1:4],truth.df.long[,2:3])
colnames(results.table)[1] <- "alpha"
results.mu0 <- results.table[,c(1,2,4,6)]
results.mu1 <- results.table[,c(1,3,5,7)]


#### mu0 results table
results.mu0$bias <-  results.mu0$mu0 - results.mu0$`mu-hat0`
results.mu0$lowerci <- results.mu0$`mu-hat0` - qnorm(.975)*results.mu0$ASE0
results.mu0$upperci <- results.mu0$`mu-hat0` + qnorm(.975)*results.mu0$ASE0
results.mu0$covers <-  ((results.mu0$lowerci <= results.mu0$mu0) &  
                          (results.mu0$mu0 <= results.mu0$upperci))
results.mu0$lowerci.t <- results.mu0$`mu-hat0` - qt(.975, m - 9)*results.mu0$ASE0
results.mu0$upperci.t <- results.mu0$`mu-hat0` + qt(.975, m - 9)*results.mu0$ASE0
results.mu0$covers.t <- ((results.mu0$lowerci.t <= results.mu0$mu0) &  
                           (results.mu0$mu0 <= results.mu0$upperci.t))
results.mu0.table <- results.mu0 %>% 
  dplyr::group_by(alpha) %>% 
  dplyr::summarise(mu0 = mean(mu0),
                   Bias = mean(bias),
                   ESE = sd(`mu-hat0`),
                   med.ASE = median(ASE0),
                   mean.ASE = mean(ASE0),
                   EC = mean(covers),
                   EC.t = mean(covers.t))
results.mu0.table$EC <- paste0(round(results.mu0.table$EC*100),"%")
results.mu0.table$EC.t <- paste0(round(results.mu0.table$EC.t*100),"%")

###### final mu0 results table, and in latex form 
results.mu0.table
print(xtable(results.mu0.table,digits=c(1,1,2,2,2,2,2,0,0)), 
      include.rownames=FALSE)


#### mu1 results table
results.mu1$bias <-  results.mu1$mu1 - results.mu1$`mu-hat1`
results.mu1$lowerci <- results.mu1$`mu-hat1` - qnorm(.975)*results.mu1$ASE1
results.mu1$upperci <- results.mu1$`mu-hat1` + qnorm(.975)*results.mu1$ASE1
results.mu1$covers <-  ((results.mu1$lowerci <= results.mu1$mu1) &  
                          (results.mu1$mu1 <= results.mu1$upperci))
results.mu1$lowerci.t <- results.mu1$`mu-hat1` - qt(.975, m - 9)*results.mu1$ASE1
results.mu1$upperci.t <- results.mu1$`mu-hat1` + qt(.975, m - 9)*results.mu1$ASE1
results.mu1$covers.t <- ((results.mu1$lowerci.t <= results.mu1$mu1) &  
                           (results.mu1$mu1 <= results.mu1$upperci.t))
results.mu1.table <- results.mu1 %>% 
  dplyr::group_by(alpha) %>% 
  dplyr::summarise(mu1 = mean(mu1),
                   Bias = mean(bias),
                   ESE = sd(`mu-hat1`),
                   med.ASE = median(ASE1),
                   mean.ASE = mean(ASE1),
                   EC = mean(covers),
                   EC.t = mean(covers.t))
results.mu1.table$EC <- paste0(round(results.mu1.table$EC*100),"%")
results.mu1.table$EC.t <- paste0(round(results.mu1.table$EC.t*100),"%")

###### final mu1 results table, and in latex form 
results.mu1.table
print(xtable(results.mu1.table,digits=c(1,1,2,2,2,2,2,0,0)), 
      include.rownames=FALSE)
