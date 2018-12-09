library(sads)
library(parallel)
source("functions.R")

## Data preparation ##
dados <- read.csv2("data.csv", as.is=TRUE)
y <- dados$N.ind
Sobs <- length(y)
## Total number of trees (average density x area)
Tot.t <- 567*5.5e8
## Proportion of total trees in the sample
p1 <- sum(dados$N.ind)/Tot.t
## Total number of plots
N.plots <- 1945
## Total area hectares
Tot.A <- 5.79e8
## Sampled area ha
Samp.A <- 2.048e3
## Amazon RAD sent by Hans
load("steege_files/AmazonRAD.RData")
## Negative binomial (already generated, too slow to rerun)
nb.pred.full <- read.csv("NB_RAD.csv")$x

## Linear regression k ~ density in the plots
## estimating k parameter of a NB for each species 
dados$dens.ha <- dados$N.ind/Samp.A
dados$k <- est.kv(mu=dados$dens.ha, nzeroes=N.plots-dados$N.plots, Nplots=N.plots)
## The regression
lm.k <-lm(log(k)~log(dens.ha), data=dados, subset=k<1)
## Estimated regression standard error
lm.k.sigma <- summary(lm.k)$sigma

## Simulating samples from teh regional RADs predicted form LogSeries and Negative Binomial
nrep <- 100 #number of repetitions
nc <- detectCores() ## Number of cores available (to parallellize)

## 1. LS RAD ##
## Expected value of the aggregation parameter for clumped sample for each species in the rad
Amazon.rad.lk <- predict(lm.k, newdata=data.frame(dens.ha=Amazon.rad/Tot.A))
## Simulated samples: random (Poisson sampling) and clumped (negative binomial using the mean aggregation parameters above)
## Matrices to store samples
tmp2 <- tmp1 <- matrix(0,nrow=length(Amazon.rad), ncol=nrep)
## The simulations
for(j in 1:nrep){
    ## Samples aggregation parameters for each species according to the linear model + standard errors of this model
    k1 <- exp(rnorm(length(Amazon.rad), 
                    mean=Amazon.rad.lk, sd=lm.k.sigma))
    ## Random (Poisson) sample
    y1 <- mcmapply(rpois2, lambda = Amazon.rad/Tot.A,
                   MoreArgs=list(N = N.plots), mc.cores=nc)
    ## Clumped (negative binomial) sample
    y2 <- mcmapply(rnbinom2, mu = Amazon.rad/Tot.A, size = k1,
                   MoreArgs=list(N = N.plots), mc.cores=nc)
    ## Pick the population sizes for which abundance in the sample is larger than zero
    tmp1[1:sum(y1>0),j] <- Amazon.rad[y1>0]
    tmp2[1:sum(y2>0),j] <- Amazon.rad[y2>0]
}


## 2. NB RAD ##
nb.pred.full.lk <- predict(lm.k, newdata=data.frame(dens.ha=nb.pred.full/Tot.A))
tmp4 <- tmp3 <- matrix(0,nrow=length(nb.pred.full), ncol=nrep)
for(j in 1:nrep){
    k2 <- exp(rnorm(length(nb.pred.full), 
                    mean=nb.pred.full.lk, sd=lm.k.sigma))
    y3 <- mcmapply(rpois2, lambda = nb.pred.full/Tot.A,
                   MoreArgs=list(N = N.plots), mc.cores=nc)
    y4 <- mcmapply(rnbinom2, mu = nb.pred.full/Tot.A, size = k2,
                   MoreArgs=list(N = N.plots), mc.cores=nc)
    tmp3[1:sum(y3>0),j] <- nb.pred.full[y3>0]
    tmp4[1:sum(y4>0),j] <- nb.pred.full[y4>0]
}

tmp4 <- tmp3 <- matrix(0,nrow=length(nb.pred.full), ncol=nrep)
for(j in 1:nrep){
    k2 <- exp(rnorm(length(nb.pred.full), 
                    mean=nb.pred.full.lk, sd=lm.k.sigma))
    y3 <- mapply(sim.occ, mu = nb.pred.full/Tot.A,
                   MoreArgs=list(N = N.plots))
    y4 <- mapply(sim.occ, mu = nb.pred.full/Tot.A, size = k2,
                   MoreArgs=list(N = N.plots, pois.samp=FALSE))
    tmp3[1:sum(y3),j] <- nb.pred.full[y3>0]
    tmp4[1:sum(y4),j] <- nb.pred.full[y4>0]
}



## Mean RADs for each simulation
ls.rnd <- apply(tmp1,1,mean)
ls.clump <- apply(tmp2,1,mean)
nb.rnd <- apply(tmp3,1,mean)
nb.clump <- apply(tmp4,1,mean)
## Save objects
save(ls.rnd, ls.clump, nb.rnd, nb.clump, file="RADs_Samples.RData")
