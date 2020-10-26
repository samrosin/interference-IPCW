##########
#
# Author: Samuel Rosin
# Revised Date: October 5, 2020
#
# This file reads in a single simulated dataset and computes the 
# IPCW point estimates and ASEs. This is done using the array method in SLURM
# so that results can be generated near-simultaneously for, e.g., 
# 1,000 different simulations. 
# 
# The group size n_i, number of groups m, and name of home directory 
# must be specified
#

#### set the following directory
user_home_directory <- "/nas/longleaf/home/srosin/Chakladar_MS_Biometrics/"
libs <- "/nas/longleaf/home/srosin/RLibs"

library(optimx, lib.loc=libs)
library(msm, lib.loc=libs) #required for optimx
library(sn, lib.loc=libs) #required for optimx
library(parfm, lib.loc=libs) #parametric frailty models
library(lme4) #mixed models
library(plyr)
library(dplyr) #data manipulation
library(rootSolve, lib.loc=libs) #required for geex
library(geex, lib.loc=libs) #compute ASEs

## collect arguments passed in from SLURM job script
args <- commandArgs(trailingOnly=TRUE)
ran_seed <- as.numeric(args[1])
m <- as.numeric(args[2])
n_i <- as.numeric(args[3])
n_sims <- as.numeric(args[4])
subdir_name <- as.character(args[5])
print(args)

arg <- Sys.getenv("SLURM_ARRAY_TASK_ID") #### get the array value

### read in correct dataset, also store an output filename
dir_path <- paste(user_home_directory,subdir_name,"/",sep="")
base_dirpath <- paste(dir_path, "base_datasets","/",sep="")
print(arg)
print(base_dirpath)

setwd(base_dirpath)
basefile <- paste("base_",arg,".csv",sep="")
base <- read.csv(basefile,header = T)


results_dirpath <- paste(dir_path, "results_datasets/",sep="")
print(results_dirpath)
dir.create(results_dirpath)
output_filename <- paste("results_",subdir_name,"_",arg,".csv",sep="")
print(output_filename)
setwd(results_dirpath)


##### specifying some simulation parameters
alphas <- c(.1,.2,.3,.4,.5,.6,.7,.8,.9)

#### define inv.logit function
inv.logit <- function(x){exp(x)/(1+exp(x))}

##### these culminate in ase_est, a function to estimate Asymptotic Standard Error
##### using sandwich variance. relies on the two geex functions, 
##### (1) an estFUN sim.estFUN and (2) an m_estimate call inside ase_est01

sim.estFUN <- function(data, prop.model, time.until){
  d.i <- sum(1-data$Delta_ij) #number of censored obs in a group
  ind.pos <- ifelse(d.i>0, 1, 0) ###indicator of positive number of censored obs
  prop.scores <- grab_psiFUN(prop.model,data)
  l <- seq(from = 0, to=d.i-1) # seq for the summation
  
  #covariates in propensity model
  data.covars <- matrix(c(rep(1,nrow(data)), data$L1_ij, data$L2_ij),nrow=3,byrow=TRUE)
  
  #covariates for censoring model
  l.tilde <- as.matrix(cbind(data$L1_ij,data$L2_ij))
  
  function(theta, alphas, time.until){
    l.tilde.theta <- data$L1_ij*theta[1] + data$L2_ij*theta[2] #fn of covs
    G0 <- theta[3]*data$X_ij ###cum hazard fn
    pt.4 <- ifelse(ind.pos==1,
                   sum(1/(theta[4]+l*theta[4]^2)),0) ### account for no censored obs 
    
    score.c1 <- sum((1-data$Delta_ij)*data$L1_ij)-
      ((d.i*theta[4]+1)*(sum(G0*exp(l.tilde.theta)*data$L1_ij)))/
      (1+ theta[4]*sum(G0*exp(l.tilde.theta)))
    
    score.c2 <- sum((1-data$Delta_ij)*data$L2_ij)-
      ((d.i*theta[4]+1)*(sum(G0*exp(l.tilde.theta)*data$L2_ij)))/
      (1+ theta[4]*sum(G0*exp(l.tilde.theta)))
    
    score.c3 <- d.i/theta[3]-
      (d.i*theta[4]+1)*sum(data$X_ij*exp(l.tilde.theta))/
      (1+ theta[4]*sum(G0*exp(l.tilde.theta)))
    
    score.c4 <- d.i/theta[4]-
      ((1/theta[4]+d.i)*sum(G0*exp(l.tilde.theta)))/
      (1+ theta[4]*sum(G0*exp(l.tilde.theta)))  +
      log(1+theta[4]*sum(G0*exp(l.tilde.theta)))/(theta[4]^2) -
      ind.pos*pt.4
    
    
    #### censoring model
    data$s.c <- data$Delta_ij*(1/(theta[4]*theta[3]*data$X_ij*exp(l.tilde.theta)+1.0))^(1/theta[4])
    
    fhats.gen <- function(alpha){
      ######## code to get fhat0 and fhat1
      data$alpha.denom <- alpha^data$A_ij*(1-alpha)^(1-data$A_ij)
      data$numer1 <- (1/data$alpha.denom)*data$A_ij*data$Delta_ij*(data$X_ij<time.until)
      data$numer0 <- (1/data$alpha.denom)*(1-data$A_ij)*data$Delta_ij*(data$X_ij<time.until)
      
      #propensity fn. b is random effect
      propensity.fn <- Vectorize(function(b){
        fixefs <- theta[5:7]
        data$h_ij <- (t(inv.logit(fixefs %*% data.covars + b))) 
        if(anyNA(data$h_ij)){return(0)} #account for propensity score of 0 or 1
        if(any(data$h_ij>.9999999999999999)){return(0)} #account for propensity score of 1
        if(any(data$h_ij<1e-322)){return(0)} #account for propensity score of 0
        dn <- dnorm(x=b,mean=0,sd=theta[8])
        if(dn==0){return(0)} 
        return(exp(sum(
          data$A_ij*log(data$h_ij/alpha)+(1-data$A_ij)*log((1-data$h_ij)/(1-alpha))
        ))*dn)
      })
      
      #generate the propensity score 
      p.hat <- stats::integrate(f=propensity.fn,lower=-Inf,upper=Inf)$value
      
      #only compute for uncensored individuals   
      data$to.sum0 <- ifelse(data$s.c==0, 0, data$numer0/(p.hat*data$s.c))
      data$to.sum1 <- ifelse(data$s.c==0, 0, data$numer1/(p.hat*data$s.c))
      
      #wrap up
      fhat0 <- (1/nrow(data))*sum(data$to.sum0)
      fhat1 <- (1/nrow(data))*sum(data$to.sum1)
      return(c(fhat0,fhat1))
    }
    
    fhats <- lapply(X=alphas, FUN=fhats.gen)
    fhats <- matrix(unlist(fhats),ncol=2,byrow=TRUE)
    
    alphas.len <- length(alphas)
    theta0.ixs <- 9:(9+alphas.len-1)
    theta1.ixs <- (9+alphas.len):(9+2*alphas.len-1)
    
    return(c(score.c1,score.c2,score.c3,score.c4,prop.scores(theta[5:8]),
             fhats[,1] - theta[theta0.ixs], fhats[,2] - theta[theta1.ixs]))
  }
}

### Estimate the vcov using geex, using the theta values from the 
### censoring model (from parfm) and propensity model ()
### Note that the variable "group" specifies the cluster 

ase_est_call <- function(base, theta, alphas, time.until, p.mixed){
  geex_obj <- m_estimate(
    estFUN = sim.estFUN,
    data = base,
    units = "group",
    compute_roots = FALSE,
    roots = theta,
    outer_args = list(prop.model = p.mixed, time.until = time.until),
    inner_args = list(alphas = alphas, time.until = time.until),
    deriv_control = setup_deriv_control(method="simple")
  )
  grab.ases <- function(a){
    if(a==0){
      ixs <- 9:(8+length(alphas)) ### values 9,8 have to do with length of theta vector, and where the causal parameters are in the vcov. careful!
    } else {
      if(a==1){ ixs <- (9+length(alphas)):(8+length(alphas)*2)} ### values 9,8 have to do with length of theta vector, and where the causal parameters are in the vcov. careful!
    } 
    return(sqrt(diag(vcov(geex_obj)[ixs,ixs])))
  }
  
  ase_0s <- grab.ases(0); ase_1s <- grab.ases(1)
  return(matrix(c(ase_0s,ase_1s),nrow=length(alphas),ncol=2))
}

#### now estimate the ASE for both muhat.0 and muhat.1 
ase_est <- function(base, alphas, theta.c, theta.h, theta.r, p.mixed,
                    mu.mat, time.until){
  thetas.c <- c(theta.c, theta.h, theta.r)
  thetas.prop <- unlist(getME(p.mixed,c('beta','theta'))) 
  ### note that grab_psiFUN wants the ranef SD, not variance, 
  ### per the help(grab_psiFUN) page
  muhat.0s <- mu.mat[,1]; muhat.1s <- mu.mat[,2]
  
  theta <- as.vector(c(thetas.c,thetas.prop,muhat.0s, muhat.1s))
  
  ases <- ase_est_call(base, theta = theta, alphas, time.until, p.mixed)
  
  return(ases)
}

############
#
# Component functions are:
# (1) fhat_i_aalpha implements the group-level estimator for a given value of a, alpha, and t
# (2) ipcw.eval.alpha evaluates the IPCW point estimator for a given dataset and value of alpha,
#     for both values a=0,1. It also evaluates the two ASEs by calling ase_est from ase_est.R
# (3) ipcw.eval.allalphas evaluates the IPCW estimators for a given dataset for all 
#     values of alpha. Here is where the censoring and propensity models are fit. 


##### This function returns the group-level causal estimator
fhat_i_aalpha <- function(grpdt, alpha, time.until = 100,
                          theta.s, theta.c, p.mixed, fixefs,
                          theta.h, theta.r){
  #numerators for the estimators under a=1 and a=0
  grpdt$alpha.denom <- alpha^grpdt$A_ij*(1-alpha)^(1-grpdt$A_ij)
  grpdt$numer1 <- (1/grpdt$alpha.denom)*grpdt$A_ij*grpdt$Delta_ij*(grpdt$X_ij<time.until)
  grpdt$numer0 <- (1/grpdt$alpha.denom)*(1-grpdt$A_ij)*grpdt$Delta_ij*(grpdt$X_ij<time.until)
  
  #propensity fn. b is random effect
  propensity.fn <- Vectorize(function(b){
    grpdt.covars <- matrix(c(rep(1,nrow(grpdt)), grpdt$L1_ij, grpdt$L2_ij),nrow=3,byrow=TRUE)
    grpdt$h_ij <- (t(inv.logit(fixefs %*% grpdt.covars + b))) 
    if(anyNA(grpdt$h_ij)){return(0)} #account for propensity score of 0 or 1
    if(any(grpdt$h_ij>.9999999999999999)){return(0)} #account for propensity score of 1
    if(any(grpdt$h_ij<1e-322)){return(0)} #account for propensity score of 0
    dn <- dnorm(x=b,mean=0,sd=sqrt(theta.s))
    if(dn==0){return(0)} 
    return(exp(sum(
      grpdt$A_ij*log(grpdt$h_ij/alpha)+(1-grpdt$A_ij)*log((1-grpdt$h_ij)/(1-alpha))
    ))*dn)
  })
  
  #generate the propensity score 
  p.hat <- stats::integrate(f=propensity.fn,lower=-Inf,upper=Inf)$value
  
  #### censoring model
  l.tilde <- as.matrix(cbind(grpdt$L1_ij,grpdt$L2_ij))
  grpdt$l.tilde.theta <- l.tilde %*% theta.c
  grpdt$s.c <- grpdt$Delta_ij*(1/(theta.r*theta.h*grpdt$X_ij*exp(grpdt$l.tilde.theta)+1.0))^(1/theta.r)
  
  #only compute for uncensored individuals   
  grpdt$to.sum0 <- ifelse(grpdt$s.c==0, 0, grpdt$numer0/(p.hat*grpdt$s.c))
  grpdt$to.sum1 <- ifelse(grpdt$s.c==0, 0, grpdt$numer1/(p.hat*grpdt$s.c))
  
  #wrap up
  fhat.0 <- (1/nrow(grpdt))*sum(grpdt$to.sum0)
  fhat.1 <- (1/nrow(grpdt))*sum(grpdt$to.sum1)
  return(c(fhat.0,fhat.1))
}

##### evaluate ipcw estimator for a specific alpha, simulated dataset called "base"
ipcw.eval.alpha <- function(base, alpha, theta.c,fixefs,theta.s,
                            theta.h, theta.r, p.mixed, time.until = 100){
  base.bygroup <- split(base,base$group)
  
  grp_avg_po_est <- lapply(X=base.bygroup, FUN=fhat_i_aalpha, alpha=alpha,
                           theta.s=theta.s, theta.c=theta.c, theta.h=theta.h, theta.r=theta.r, 
                           fixefs=fixefs, time.until = time.until)
  grp_avg_po_est <- matrix(unlist(grp_avg_po_est), ncol=2, byrow = T)
  pop_avg_po_est <- matrix(apply(grp_avg_po_est, 2, mean),ncol=2)
  
  return(pop_avg_po_est)
}

#### evaluate the ipcw estimator for all levels of alpha, 1 dataset
ipcw.eval.allalphas <- function(base, alphas, time.until){
  pf <- parfm(Surv(X_ij,1-Delta_ij) ~ L1_ij + L2_ij, cluster="group",
              data=base, dist="exponential", frailty="gamma")
  theta.c <- as.matrix(coef(pf),nrow=2,ncol=1)
  p.mixed <- glmer(A_ij ~ L1_ij + L2_ij + (1 | group), data=base, family=binomial)
  if(isSingular(p.mixed)){
    all.results <- matrix(data=NA, nrow=length(alphas), ncol=4,
                        dimnames=list(alphas,
                                      c("mu-hat0","mu-hat1","ASE0","ASE1")))
    return(all.results)
  }
  fixefs <- matrix(fixef(p.mixed),nrow=1)
  theta.s <- (getME(p.mixed,"theta"))^2 #square it b/c it returns as sd 
  theta.h <- pf[,"ESTIMATE"]["lambda"] #baseline hazard
  theta.r <- pf[,"ESTIMATE"]["theta"] #gamma frailty 
  
  all.ipcws <- lapply(X=alphas,FUN=ipcw.eval.alpha, base=base,theta.c=theta.c,
                      fixefs=fixefs, theta.s=theta.s, theta.h=theta.h, theta.r=theta.r, 
                      p.mixed=p.mixed, time.until=time.until)
  all.ipcws <- matrix(unlist(all.ipcws),ncol=2,byrow=TRUE)
  ases <- ase_est(base,alphas,theta.c,theta.h,theta.r,p.mixed,
                  mu.mat = all.ipcws,time.until = 100)
  all.results <- cbind(all.ipcws, ases)
  print(all.results)
  
  rownames(all.results) <- alphas
  colnames(all.results) <- c("mu-hat0","mu-hat1","ASE0","ASE1")
  return(all.results)
}

results <- ipcw.eval.allalphas(base = base, alphas = alphas, time.until=100)
write.csv(x = results, file = output_filename)
