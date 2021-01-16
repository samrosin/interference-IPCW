##########
#
# Author: Samuel Rosin
# Date: October 12, 2020
#
# This file generates simulated datasets and stores them in a subdirectory 
# entitled "base_datasets". A single long-format dataset of the true estimand values
# is also produced and stored in a .csv file. 
# 
# The random effects b_i in the propensity model (i.e. which have an effect on treatment propensity)
# are generated from a Gumbel distribution 
#

#### these two directories need to be specified 
libs <- "/nas/longleaf/home/srosin/RLibs"
user_home_directory <- "/nas/longleaf/home/srosin/Chakladar_MS_Biometrics/"

library(plyr)
library(dplyr) #data manipulation

## collect arguments passed in from SLURM job script
args <- commandArgs(trailingOnly=TRUE)
print(args)
ran_seed <- as.numeric(args[1])
m <- as.numeric(args[2])
n_i <- as.numeric(args[3])
n_sims <- as.numeric(args[4])
subdir_name <- as.character(args[5])

set.seed(ran_seed)
alphas <- c(.1,.2,.3,.4,.5,.6,.7,.8,.9)

#### define inv.logit function
inv.logit <- function(x){exp(x)/(1+exp(x))}

sim_dataset <- function(n_i = n_i,
                        m = m
){
  n <- n_i*m #total sample size 
  
  
  df <- data.frame(group = rep(seq(from = 1, to = n/n_i, by=1),each=n_i), n_i = n_i,
                   V_ij=rep(NA,n), r1_i=rep(NA,n), r2_ij=rep(NA,n), 
                   L1_ij=rep(NA,n), log.L2_ij=rep(NA,n), L2_ij=rep(NA,n)) 
  
  #### Step 1: Generate two baseline covariates
  df$V_ij <- rexp(n=n, rate=1/20)
  
  r1_is <- rnorm(n=n/n_i, mean=0, sd=.1) #sd = sqrt(.1) = .316
  df$r1_i <- rep(r1_is, each=n_i)
  df$r2_ij <- rnorm(n=n, mean=0, sd=.1)
  
  df$L1_ij <- pmin(df$V_ij + df$r1_i + df$r2_ij, 100)/10 #age
  
  df$log.L2_ij <- rnorm(n=n, mean=df$r1_i + df$r2_i, sd = 0.75)
  df$L2_ij <- exp(df$log.L2_ij)
  
  
  ##### Step 2: Generate random effects b_i for treatment model. 
  #####         There are i=500 groups each with sample size n_i=10.
  #####         Each individual in a group has the same treatment and 
  #####         same random effect.
  #### Here they are generated from a Uniform(-1.8, 1.8) distribution
  #### with mean 0 and variance 1.08
  ####
  ####
  b_is <- runif(n=n/n_i, min = -1.8, max = 1.8)
  
  #b_is <- rnorm(n=n/n_i, mean=0, sd=sqrt(1.0859)) #sd=sqrt(1.0859)=1.04
  df$b_i <- rep(b_is, each=n_i) #first 10 are b_1, next 10 are b_2, ...
  
  ##### Step 3: Treatment indicators A_ij from Bernoulli(p_ij)
  p_ij <- inv.logit(0.2727 - 0.0387*df$L1_ij + 0.2179*df$L2_ij + df$b_i)
  df$A_ij <- rbinom(n=n, size=1, prob=p_ij)
  
  ##### Step 4: Potential times-to-event T_ij(a_i)
  #####         Sampled from Exponential with mean mu_ij 
  #####         <--> rate 1/mu_ij
  #####         This function is heavily based off of the generate_po() function 
  #####         found on Bradley Saul's interferenceSim Github repository:
  #####         https://github.com/bsaul/interferenceSim
  
  df.new <- data.frame(id=seq(from=1,to=n),group=df$group,L1_ij= df$L1_ij,L2_ij= df$L2_ij)
  generate_po <- function(base_data, 
                          parameters = c(200, 100, -0.98, -0.145, 50)){
    
    ## Generates a vector of potential outcomes for a given 
    ## treatment level and number of other subjects in group treated
    po_y_ij <- function(grpdt, a, k, parameters){
      
      n_i <- nrow(grpdt)
      alpha <- k / n_i
      X <- cbind(1, a, grpdt$L1_ij, grpdt$L2_ij, alpha)
      
      mu_ij <- X %*% parameters
      
      if(length(mu_ij) != n_i){stop('Problem in po_y_ij(): length(mu_ij) != n_i')}
      return(rexp(n=n_i, rate=1/mu_ij))
      
    }
    
    ## Generates a matrix of potential outcomes for all treatment levels in 
    ## and all numbers of other subjects treated (from 0 to n_i - 1)
    gen_po_ij <- function(grpdt, parameters){
      n_i <- nrow(grpdt)
      ids <- grpdt$id
      po <- array(dim = c(n_i, n_i, 2), 
                  dimnames = list(id = ids, k = 0:(n_i-1), trt = c(0,1)))
      
      for(k in 0:(n_i-1)){
        po[ , as.character(k), 1] <- po_y_ij(grpdt, 0, k, parameters)
        po[ , as.character(k), 2] <- po_y_ij(grpdt, 1, k, parameters)
      }
      
      return(po)
    }
    
    ## Generate potential outcomes for each group ##
    po <- lapply(split(base_data, base_data$group), 
                 gen_po_ij, parameters = parameters)
    return(po)
  }
  
  base.parameters <- c(200, 100, -0.98, -0.145, 50)
  this.po <- generate_po(df.new,base.parameters)
  
  
  ##### Step 5: Random effects for censoring model e_i
  #####         Sampled from Gamma with mean 1, variance 1.25
  
  ### drop some columns at this point
  drops <- c("V_ij","r1_i","r2_ij","log.L2_ij")
  df <- df[,!names(df) %in% drops]
  
  # this gamma distribution in R is specified with shape=a, scale=s
  # where mean=a*s=1 and variance = a*s^2=1.25
  # Solving algebraically, we find a=0.8 and s=1.25. 
  e_i <- rgamma(n=n/n_i, shape=0.8, scale=1.25)
  df$e_i <- rep(e_i, each=n_i)
  
  ##### Step 6: Censoring times C_ij sampled from Exponential
  #####         with mean 1/lambda_0 <--> rate lambda_0
  lambda_0 <- 0.015 * exp(0.002*df$L1_ij + 0.015*df$L2_ij)*df$e_i
  df$C_ij <- rexp(n=n, rate=lambda_0)
  
  #### Step 7: Individual censoring indicators Delta_ij
  df.bygrp <- split(df,df$group)
  delta_gen <- function(group.df){
    group.df$Delta_ij <- NA
    group.df$X_ij <- NA
    grp <- unique(group.df$group)
    group.df$k <- NA
    
    #create k
    for(j in 1:nrow(group.df)){
      all.other.rows <- group.df[-j,]
      group.df[j,]$k <- sum(all.other.rows$A_ij)
    }
    
    for(i in 1:nrow(group.df)){
      delta <- this.po[[grp]][i,group.df[i,]$k+1,(group.df[i,]$A_ij+1)] < group.df[i,]$C_ij ## make sure k+1, correct index
      group.df[i,]$Delta_ij <- ifelse(delta==TRUE,1,0)
      group.df[i,]$X_ij <- min(this.po[[grp]][i,group.df[i,]$k+1,(group.df[i,]$A_ij+1)],group.df[i,]$C_ij) ## make sure k+1, correct index
    }
    return(group.df)
  }
  df.bygrp <- lapply(X=df.bygrp, FUN=delta_gen) #get delta indicator 
  df <- do.call(rbind,df.bygrp) #'unlist' the list of by-group dataframes
  
  ####drop a few columns. these can be added back for debugging or, in the case of n_i, 
  #### to implement varying group sizes
  df[,c("n_i","b_i","e_i")] <- NULL
  
  out <- list(df,this.po)
  return(out)
}
sim.datasets <- replicate(n=n_sims,sim_dataset(n=n_i,m=m))
base.datasets <- sim.datasets[1,]; potential.outcomes.list <- sim.datasets[2,]
base.datasets <- mapply(FUN=cbind, base.datasets, sim.number=as.list(1:n_sims),SIMPLIFY=FALSE) #add sim number

### create a directory for these sims
dir_path <- paste(user_home_directory,subdir_name,"/",sep="")
dir.create(dir_path)

#create subdirectory for base datasets and then print each dataset there
base_dirpath <- paste(dir_path, "base_datasets","/",sep="")
dir.create(base_dirpath)
setwd(base_dirpath)

for(i in 1:n_sims){
  base_filename <- paste("base_",i,".csv",sep="")
  write.csv(x = base.datasets[[i]], file = base_filename, row.names = FALSE)
}



#----------------------------#
#
# The following functions culminate in calc_estimands(), which
# is used to return 
# mu(100,0,alpha) and mu(100,1,alpha)
# for several values of alpha. 
#
#--------------------------#

#-----------------------------------------------------------------------------#
#' Individual level causal estimand mu_ij(t, a, alpha)
#' 
#' Takes n X n  x 2 array as input
#' @return n x 2 matrix with values of f-bar_ij(t, a, alpha) 
#' 
#-----------------------------------------------------------------------------#
fbar_i_aalpha <- function(outcomes, alpha, time.until=100){
  
  ## Necessary bits ##
  nn <- dim(outcomes)[1]
  kk <- 0:(nn - 1)
  aa <- alpha
  
  ## Create array for output ##
  out <- array(dim = c(nn, 2), dimnames = list(id = dimnames(outcomes)$id,
                                               trt = dimnames(outcomes)$trt))
  
  fbar <- function(po.times, alpha, time.until=100){
    ind <- ifelse(po.times < time.until, 1, 0)
    return(sum(ind * choose(nn - 1, kk) * aa^(kk) * (1 - aa)^(nn - kk - 1)))
  }
  out[ , 1] <- apply(outcomes[ , , 1], 1, fbar) 
  out[ , 2] <- apply(outcomes[ , , 2], 1, fbar)
  
  return(out)
}

#-----------------------------------------------------------------------------#
#' Individual level marginal causal estimand f-bar(t, alpha)
#' Takes the n x 2 matrix returned by fbar_ij_aalpha
#' 
#' @return n x 1 vector of fbar_ij(alpha) 
#-----------------------------------------------------------------------------#
fbar_i_alpha <- function(input, alpha){
  alpha*input[, 2] + (1 - alpha) * input[, 1]
}

#-----------------------------------------------------------------------------#
#' Calculate estimands
#' 
#' Takes a set of potential outcomes generated by \code{\link{sim_dataset}} and 
#' returns a data frame containing the population level causal estimands
#' 
#' @param potential_outcomes a list of potential outcomes generated by 
#' \code{\link{generate_po}}
#' @param alphas a vector of 'allocation schemes' in (0, 1)
#' @param time.until t, where the estimand is defined in terms of the contrast 
#' in the risk of having an event by time t. Currently defined for a single value. 
#' @return a data frame with 4 columns: 1) alpha (the allocation scheme) 2)
#' a0 - the outcome estimand for untreated, 3) a1 - the outcome estimand for
#' treated, 4) marg - the marginal outcome estimand
#' @examples 
#' calc_estimands(sim_dataset()[[2]], seq(.2, .8, by = .1), time.until = 100)
#' @export
#-----------------------------------------------------------------------------#

calc_estimands <- function(potential_outcomes, alphas, time.until = 100){
  
  ## Calculate population level outcomes based on 
  ## individual --> group --> population
  estimands <- function(alpha, potential_outcomes){
    
    # --- Individual average potential outcomes --- #  
    ind_avg_po <- lapply(potential_outcomes, function(x) fbar_i_aalpha(x, alpha, time.until))
    
    # --- Group average potential outcomes --- #                
    grp_avg_po <- lapply(ind_avg_po, function(x) apply(x, 2, mean))
    grp_avg_po <- matrix(unlist(grp_avg_po), ncol=2, byrow = T)
    # --- Population average potential outcomes --- # 
    pop_avg_po <- apply(grp_avg_po, 2, mean)
    
    # --- Marginal potential outcomes --- #
    ind_marg_po <- lapply(ind_avg_po, function(x) fbar_i_alpha(x, alpha))
    grp_marg_po <- unlist(lapply(ind_marg_po, mean))
    pop_marg_po <- mean(grp_marg_po)
    
    return(data.frame(alpha = alpha, 
                      mu0    = pop_avg_po[1],
                      mu1    = pop_avg_po[2],
                      mu.marg  = pop_marg_po))
  }
  
  ## Apply the above function to each alpha level ##
  true_effects <- lapply(alphas, estimands, 
                         potential_outcomes = potential_outcomes)
  
  truth <- do.call(rbind, true_effects)
  return(truth)
}


truth.list <- lapply(X=potential.outcomes.list, FUN=calc_estimands,
                     alphas=alphas, time.until=100)
truth.df.long <- do.call(what=rbind,args=truth.list)
truth.avgs.mu0 <- aggregate(formula = mu0 ~ alpha, data=truth.df.long, FUN=mean)
truth.avgs.mu1 <- aggregate(formula = mu1 ~ alpha, data=truth.df.long, FUN=mean)
truth.avgs <- merge(truth.avgs.mu0, truth.avgs.mu1, by="alpha")

#print truthtable at sims directory level
setwd(dir_path)
write.csv(x = truth.avgs, file = paste("truthtable_",subdir_name,".csv",sep=""))

