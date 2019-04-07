library(parallel)
source("functions.R")
load("lists_with_all_objects.RData")


## Number of cores to use
mc.cores <- 3

################################################################################
## Bias of logseries method
################################################################################
## Checking the relationship between Estimated and real S
S2 <- round(runif(100, 1e4, 2.5e4))

sim.ls.rad.13 <- mclapply(S2, sim.rad, N = atdn.13$Tot.t, sad = "ls", nb.fit = atdn.13$y.nb2,
                          mc.cores = mc.cores, upper = 1e12)

sim.ls.samp.13 <- mclapply(sim.ls.rad.13, sim.radsamp, tot.area  = atdn.13$Tot.A,
                           n.plots = atdn.13$N.plots, lmk.fit = atdn.13$lm.k,
                           nb.fit = atdn.13$y.nb2, mc.cores = mc.cores)
sim.ls.estS.13 <- data.frame(S = S2,
                             S.est.rnd = sapply(sim.ls.samp.13, function(x) ls.estS(x$rnd.samp, N = atdn.13$Tot.t)),
                             S.est.clump = sapply(sim.ls.samp.13, function(x) ls.estS(x$clump.samp, N = atdn.13$Tot.t))
                             )

## Bias at the confidence interval of observed estimated S richness
## Bias at the range of the confidence interval of the observed estimate ##
S2 <- round( runif(100, atdn.13$S.ls.ci[1], atdn.13$S.ls.ci[2]) )
## Log-series

sim.ls.rad.13 <- mclapply(S2, sim.rad, N = atdn.13$Tot.t, sad = "ls", nb.fit = atdn.13$y.nb2,
                          mc.cores = mc.cores, upper = 1e12)

sim.ls.samp.13 <- mclapply(sim.ls.rad.13, sim.radsamp, tot.area  = atdn.13$Tot.A,
                           n.plots = atdn.13$N.plots, lmk.fit = atdn.13$lm.k,
                           nb.fit = atdn.13$y.nb2, mc.cores = mc.cores)
sim.ls.estS.13 <- data.frame(S = S2,
                             S.est.rnd = sapply(sim.ls.samp.13, function(x) ls.estS(x$rnd.samp, N = atdn.13$Tot.t)),
                             S.est.clump = sapply(sim.ls.samp.13, function(x) ls.estS(x$clump.samp, N = atdn.13$Tot.t))
                             )



################################################################################
## Bias of negative binomial method
################################################################################
f3 <- function(x, size, loglink=TRUE){
    x <- sort(x[x>0])
    result <- try(tovo(fitnbinom2(x, start.value=c(size=size, mu=mean(x))),
                       p=sum(x)/atdn.13$Tot.t, loglink=loglink))
    if(class(result)=="try-error")
        return(NA)
    else
        return(result)
}

## Bias over a wide range species richness
S2 <- round(runif(100, 1e4, 2e4))
sim.tnb.rad.13 <- mclapply(S2, sim.rad, N = atdn.13$Tot.t, sad = "tnb",
                                       nb.fit = atdn.13$y.nb2, mc.cores = mc.cores, upper = 1e20)
sim.tnb.samp.13 <- mclapply(sim.tnb.rad.13, sim.radsamp, tot.area  = atdn.13$Tot.A,
                           n.plots = atdn.13$N.plots, lmk.fit = atdn.13$lm.k,
                           nb.fit = atdn.13$y.nb2, mc.cores = mc.cores)
sim.tnb.estS.13 <- data.frame(S = S2,
                              S.est.rnd = sapply(sim.tnb.samp.13,
                                                 function(x) f3(x$rnd.samp, size = 0.018) ),
                              S.est.clump = sapply(sim.tnb.samp.13,
                                                 function(x) f3(x$clump.samp, size = 0.018) )
                              )
# Storing in lists
bias13 <- list(
    ls = list(rads=sim.ls.rad.13, samples = sim.ls.samp.13, estimates = sim.ls.estS.13),
    tnb= list(rads=sim.tnb.rad.13, samples = sim.tnb.samp.13, estimates = sim.tnb.estS.13)
    )
    
    

## Truncated negative binomial
## 2013
S2 <- round( runif(250, atdn.13$tovo.S$CIs[4,2], atdn.13$tovo.S$CIs[4,1]) )
sim.tnb.rad.13 <- mclapply(S2, sim.rad, N = atdn.13$Tot.t, sad = "tnb",
                                       nb.fit = atdn.13$y.nb2, mc.cores = mc.cores, upper = 1e20)
sim.tnb.samp.13 <- mclapply(sim.tnb.rad.13, sim.radsamp, tot.area  = atdn.13$Tot.A,
                           n.plots = atdn.13$N.plots, lmk.fit = atdn.13$lm.k,
                           nb.fit = atdn.13$y.nb2, mc.cores = mc.cores)
sim.tnb.estS.13 <- data.frame(S = S2,
                              S.est.rnd = unlist(mclapply(sim.tnb.samp.13,
                                                 function(x) f3(x$rnd.samp, size = 0.018) , mc.cores = mc.cores)),
                              S.est.clump = unlist(mclapply(sim.tnb.samp.13,
                                                            function(x) f3(x$clump.samp, size = 0.018) , mc.cores = mc.cores))
                              )


