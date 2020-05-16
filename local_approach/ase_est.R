##########
#
# Author: Samuel Rosin
# Date: May 7, 2020
# 
# This file evaluates the asymptotic standard errors of given estimators. 
#
##### These functions culminate in ase_est, a function to estimate Asymptotic Standard Error
##### using sandwich variance. The approach relies on the two geex functions, 
##### (1) an estFUN, sim.estFUN and (2) an m_estimate call inside ase_est_call

sim.estFUN <- function(data, prop.model, time.until){
  d.i <- sum(1-data$Delta_ij) #number of censored obs in a group
  ind.pos <- ifelse(d.i>0, 1, 0) ###indicator of positive number of censored obs
  prop.scores <- grab_psiFUN(prop.model,data)
  l <- seq(from = 0, to=d.i-1) # seq for the summation
  
  data.covars <- matrix(c(rep(1,nrow(data)), data$L1_ij, data$L2_ij),nrow=3,byrow=TRUE)
  
  #for censoring model
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
### censoring model (from parfm) and propensity model.
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
  
  ### fn to grab the relevant ASEs from the returned vcov matrix 
  grab.ases <- function(a){
    if(a==0){
      ixs <- 9:(8+length(alphas))
    } else {
      if(a==1){ ixs <- (9+length(alphas)):(8+length(alphas)*2)}
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
  
  ases <- ase_est_call(base, theta = theta, alphas, time.until = time.until, p.mixed)
  
  return(ases)
}