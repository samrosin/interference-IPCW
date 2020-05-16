n_i <- 50; m <- 500; n_sims <- 1000
user_home_directory <- "/nas/longleaf/home/srosin/" #home directory on a SLURM-managed computing cluster

library(dplyr)
library(xtable)
#### 
setwd(paste(user_home_directory,"truthtables",sep=""))
truthtable <- read.csv(paste("truthtable_ni",n_i,"_m",m,".csv",sep=""),header=T)
results_wd <- paste("/nas/longleaf/home/srosin/results_datasets/results_ni",
                    n_i,"_m",m,sep="")
setwd(results_wd)

alpha <- rep(NA,9*n_sims); mu.hat0 <- rep(NA, 9*n_sims)
mu.hat1 <- rep(NA,9000); ASE0 <- rep(NA,9*n_sims); ASE1 <- rep(NA, 9*n_sims)
results.table <- data.frame(alpha,mu.hat0,mu.hat1,ASE0,ASE1)

i <- 1
for(arg in 1:n_sims){
  filename <- paste("results_ni",n_i,"_m",m,"_",arg,".csv",sep="")
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
results.mu0$lowerci.t <- results.mu0$mu.hat0 - qt(.975, m - 9)*results.mu0$ASE0
results.mu0$upperci.t <- results.mu0$mu.hat0 + qt(.975, m - 9)*results.mu0$ASE0
results.mu0$covers.t <- ((results.mu0$lowerci.t <= results.mu0$mu0) &  
                           (results.mu0$mu0 <= results.mu0$upperci.t))
results.mu0.table <- results.mu0 %>% 
  dplyr::group_by(alpha) %>% 
  dplyr::summarise(mu0 = mean(mu0),
                   Bias = mean(bias),
                   ESE = sd(mu.hat0),
                   med.ASE = median(ASE0),
                   mean.ASE = mean(ASE0),
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
                   med.ASE = median(ASE1),
                   mean.ASE = mean(ASE1),
                   EC = mean(covers),
                   EC.t = mean(covers.t))
results.mu1.table$EC <- paste0(round(results.mu1.table$EC*100),"%")
results.mu1.table$EC.t <- paste0(round(results.mu1.table$EC.t*100),"%")

###### print out the results to a table 
setwd(paste(user_home_directory,"output_tables",sep=""))
outfile_name <- paste("output_tables_ni",n_i,"_m",m,".Rout",sep="")

sink(file = outfile_name, append = TRUE, type = "output")
print(xtable(results.mu0.table,digits=c(1,1,2,2,2,2,0,0)), 
      include.rownames=FALSE)
print(xtable(results.mu1.table,digits=c(1,1,2,2,2,2,0,0)), 
      include.rownames=FALSE)
sink()


