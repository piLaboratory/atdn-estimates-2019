source("../functions.R")
## library(abc)
## library(parallel)
## library(untb)
## library(sads)
load("needed_objs2013.RData")

## Simulations
## Function to parallelize simulations
f1 <- function(x, ...){
    y <- try(sim.abc(S = x,
                     N = needed.objs$Tot.t,
                     n.plots = needed.objs$N.plots,
                     tot.area= needed.objs$Tot.A,
                     nb.fit = needed.objs$y.nb2,
                     lm.sd.fit = needed.objs$lm.sd,
                     lmk.fit = needed.objs$lm.k,
                     nrep = 1, ...))
    if(class(y)=="try-error")
        return(matrix(NA, nrow=2, ncol=8))
    else
        return(y)
}


## Number of Simulated values
nsims <- 5e3
##simulated.vals <- runif(nsims, 1e4, S.orc$S.est)
simulated.vals <- runif(nsims, 1e4, 2e4)
## Samples of LS
LS.sims <- mclapply(simulated.vals, f1, sad = "ls",  lower=1e-20, upper=1e20, mc.cores=10)
## save.image()
## Samples of TNB
NB.sims <- mclapply(simulated.vals, f1, sad = "tnb", lower=1e-20, upper=1e20, mc.cores=10)
##save.image()
## Samples of Log-normal
LN.sims <- mclapply(simulated.vals, f1, sad = "lnorm", sdlog = needed.objs$pln.cf[2], lower=1e-20, upper=1e20, mc.cores=10)
##save.image()

## Assembles all simulation results in a matrix
## Simulations from LS rads ##
## exclude failed simulations with random sampling
LS.index <- sapply(LS.sims, function(x) !any(is.na(x)))
j1 <- (1:length(LS.sims))[LS.index]
## Assembles a matrix, each line the result of a simulation
all.sims <- LS.sims[[min(j1)]]
for(i in j1[-min(j1)])
    all.sims <- rbind(all.sims, LS.sims[[i]])
## Simulations from TNB rad
## Excluding failed simulations
NB.index <- sapply(NB.sims, function(x) !any(is.na(x)))
j2 <- (1:length(NB.sims))[NB.index]
## rbind the results of the simulations to the same results matrix
for(i in j2)
    all.sims <- rbind(all.sims, NB.sims[[i]])
## Simulations ftom LN rad
## Simulations from TNB rad
## Excluding failed simulations
LN.index <- sapply(LN.sims, function(x) !any(is.na(x)))
j3 <- (1:length(LN.sims))[LN.index]
## rbind the results of the simulations to the same results matrix
for(i in j3)
    all.sims <- rbind(all.sims, LN.sims[[i]])
## Vector with labels for each simulation
sim.ids <- c(rep(c("LSrnd", "LSclump"), sum(LS.index)),
             rep(c("NBrnd", "NBclump"), sum(NB.index)),
             rep(c("LNrnd", "LNclump"), sum(LN.index)))## Simulated values for each simulation
## Vector with the parameters used in each simulation
## (which is the total species richness in the regional RADs)
sim.y <- c(rep(simulated.vals[j1],each=2),
           rep(simulated.vals[j2],each=2),
           rep(simulated.vals[j3],each=2))

save(needed.objs, all.sims, sim.ids, sim.y, file="ABC2013.RData")
