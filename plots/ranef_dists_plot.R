#### check mean, var, skewness
# library(e1071)
# library(evd)
# n_sims <- 10^7
# 
# gumbel.vec <- rgumbel(n_sims, loc = -0.4689848, scale = 0.8124949)
# m <- rbinom(n_sims, 1, .5)
# normal.mix.vec <- (1-m)*rnorm(n_sims,-2, sd = 1.48) + m*rnorm(n_sims, 2, sd=1.48)
# exp.vec <- rexp(n_sims, rate = 1/1.0421) - 1.0421
# beta.vec <- 6*rbeta(n_sims, 2, 4) - 2
# beta.vec.bimod <- 2.8*rbeta(n_sims, .2, .4) - (2.8*.2/.6)
# unif.vec <- runif(n_sims, -1.8, 1.8)
# 
# mean(gumbel.vec); var(gumbel.vec); skewness(gumbel.vec)
# mean(normal.mix.vec); var(normal.mix.vec); e1071::skewness(normal.mix.vec)
# mean(exp.vec); var(exp.vec); e1071::skewness(exp.vec) ###skewness 2 
# mean(beta.vec); var(beta.vec); e1071::skewness(beta.vec) ##skewness .467
# mean(beta.vec.bimod); var(beta.vec.bimod); e1071::skewness(beta.vec.bimod) ##skewness .689      
# mean(unif.vec); var(unif.vec); skewness(unif.vec)

### define the pdf of X, where X=mult*B and B is Beta(shape1, shape2)
dbeta_mult <- function(x, shape1, shape2, mult){
    density <- rep(NA, length(x))
    for(i in 1:length(x)){
        if(x[i] < 0){density[i] <- 0}
        else if(x[i] > 6){density[i] <- 0}
        else density[i] <- (1/mult)*dbeta(x[i]/mult,shape1,shape2)
    }
    density
}


######## plot the distributions
x <- seq(-4, 4, length=10000)
y <- seq(-4, 4, length=250)
normal.dist <- dnorm(x, mean=0, sd=sqrt(1.0859))
logistic.dist <- dlogis(x, location=0, scale=0.574521)
t.dist <- dt(x, df=10)
gumbel.dist <- dgumbel(x, loc = -0.4689848, scale = 0.8124949)

normal.mix.dist <- .5*dnorm(x,-1, sd = 0.3) + .5*dnorm(x, 1, sd=0.3)
exp.dist <- dexp(x +1.0421, rate = 1/1.0421) 
beta.dist <- dbeta_mult(x + 2, 2, 4, 6) 
beta.dist.bimodal <- dbeta_mult(y +(2.8*.2/.6), .2, .4, 2.8) 
unif.dist <- dunif(x,-1.8,1.8)

###plot to png or pdf
setwd("/Users/samuelrosin/Dropbox/_UNC/CIWI/_Sujatro Manuscript/Revision_3_Jan2021")
# png(file = "ranef_dists.png",
#     units = "in",
#     width = 10,
#     height = 10,
#     res= 600)
pdf(file = "ranef_dists.pdf",
    width = 10,
    height = 10)

par(mfrow=c(3,3))
plot(x, normal.dist, type='l', lty=1,
     xlab=expression(b[i]),ylab="density", main = "Normal")
plot(x, logistic.dist, type='l', lty=1,
     xlab=expression(b[i]),ylab="density", main = "Logistic")
plot(x, t.dist, type='l', lty=1, 
     xlab=expression(b[i]),ylab="density", main = "Student's-t")
plot(x, gumbel.dist, type='l', lty=1, 
     xlab=expression(b[i]),ylab="density", main = "Gumbel")
plot(x, normal.mix.dist, type='l', lty=1,
     xlab=expression(b[i]),ylab="density", main = "50:50 Normal Mixture")
plot(x, exp.dist, type='l', lty=1, 
     xlab=expression(b[i]),ylab="density", main = "Exponential")
plot(x, beta.dist, type='l', lty=1, 
     xlab=expression(b[i]),ylab="density", main = "Beta")
plot(y, beta.dist.bimodal, type='l', lty=1,
     xlab=expression(b[i]),ylab="density", main = "Beta (bimodal)")
plot(x, unif.dist, type='l', lty=1,
     xlab=expression(b[i]),ylab="density", main = "Uniform(-1.8, 1.8)")

dev.off()


