##########
#
# Author: Samuel Rosin
# Date: May 7, 2020
# 
# This file evaluates the IPCW estimator for a given dataset
#
# Component functions are:
# (1) fhat_i_aalpha implements the group-level estimator for a given value of a, alpha, and t
# (2) ipcw.eval.alpha evaluates the IPCW point estimator for a given dataset and value of alpha,
#     for both values a=0,1, by calling fhat_i_aalpha repeatedly. 
#     It also evaluates the two ASEs by calling ase_est from ase_est.R
# (3) ipcw.eval.allalphas evaluates the IPCW estimators for a given dataset for all 
#     values of alpha by calling ipcw.eval.alpha repeatedly. 
#     Here is where the censoring and propensity models are fit. 


##### This function returns the group-level causal estimator
fhat_i_aalpha <- function(grpdt, alpha, time.until = time.until,
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
                            theta.h, theta.r, p.mixed, time.until = time.until){
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
              data=base, dist="exponential", frailty="gamma") ##frailty model
  theta.c <- as.matrix(coef(pf),nrow=2,ncol=1) 
  p.mixed <- glmer(A_ij ~ L1_ij + L2_ij + (1 | group), data=base, family=binomial) ##mixed propensity model
  ##### below code could be used when there are issues fitting the mixed model
  # if(isSingular(p.mixed)){
  #   all.ipcws <- matrix(data=NA, nrow=9, ncol=2,
  #                       dimnames=list(c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
  #                                     c("mu-hat0","mu-hat1")))
  #   return(all.ipcws)
  # }
  fixefs <- matrix(fixef(p.mixed),nrow=1)
  theta.s <- (getME(p.mixed,"theta"))^2 #square it b/c it is given as SD, but we use the variance
  theta.h <- pf[,"ESTIMATE"]["lambda"] #baseline hazard
  theta.r <- pf[,"ESTIMATE"]["theta"] #gamma frailty 
  
  all.ipcws <- lapply(X=alphas,FUN=ipcw.eval.alpha, base=base,theta.c=theta.c,
                      fixefs=fixefs, theta.s=theta.s, theta.h=theta.h, theta.r=theta.r, 
                      p.mixed=p.mixed, time.until=time.until)
  all.ipcws <- matrix(unlist(all.ipcws),ncol=2,byrow=TRUE)
  ases <- ase_est(base,alphas,theta.c,theta.h,theta.r,p.mixed,
                  mu.mat = all.ipcws,time.until = time.until)
  all.results <- cbind(all.ipcws, ases)   
  print(all.results)
  
  rownames(all.results) <- alphas
  colnames(all.results) <- c("mu-hat0","mu-hat1","ASE0","ASE1")
  return(all.results)
}