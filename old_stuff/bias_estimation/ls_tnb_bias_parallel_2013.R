library(parallel)
source("../functions.R")
load("../lists_with_all_objects.RData")


## Number of cores to use
mc.cores <- 12

################################################################################
## Bias of logseries method
################################################################################
## Checking the relationship between Estimated and real S
S1 <- round(runif(100, 1e4, 2.5e4))

sim.ls.rad <- mclapply(S1, sim.rad, N = atdn.13$Tot.t, sad = "ls", nb.fit = atdn.13$y.nb2,
                          mc.cores = mc.cores, upper = 1e12)
sim.ls.samp <- mclapply(sim.ls.rad, sim.radsamp, tot.area  = atdn.13$Tot.A,
                           n.plots = atdn.13$N.plots, lmk.fit = atdn.13$lm.k,
                           nb.fit = atdn.13$y.nb2, mc.cores = mc.cores)
sim.ls.estS <- data.frame(S = S1,
                             S.est.rnd = unlist( mclapply(sim.ls.samp,
                                                         function(x) ls.estS(x$rnd.samp, N = atdn.13$Tot.t),
                                                         mc.cores = mc.cores)),
                             S.est.clump = unlist( mclapply(sim.ls.samp,
                                                         function(x) ls.estS(x$clump.samp, N = atdn.13$Tot.t),
                                                         mc.cores = mc.cores))
                             )
## Bias at the range of the confidence interval of the observed estimate ##
S2 <- round( runif(100, atdn.13$S.ls.ci[1], atdn.13$S.ls.ci[2]) )

sim.ls.rad.ci <- mclapply(S2, sim.rad, N = atdn.13$Tot.t, sad = "ls", nb.fit = atdn.13$y.nb2,
                          mc.cores = mc.cores, upper = 1e12)

sim.ls.samp.ci <- mclapply(sim.ls.rad.ci, sim.radsamp, tot.area  = atdn.13$Tot.A,
                           n.plots = atdn.13$N.plots, lmk.fit = atdn.13$lm.k,
                           nb.fit = atdn.13$y.nb2, mc.cores = mc.cores)
sim.ls.estS.ci <- data.frame(S = S2,
                             S.est.rnd = unlist( mclapply(sim.ls.samp.ci,
                                                         function(x) ls.estS(x$rnd.samp, N = atdn.13$Tot.t),
                                                         mc.cores = mc.cores)),
                             S.est.clump = unlist( mclapply(sim.ls.samp.ci,
                                                         function(x) ls.estS(x$clump.samp, N = atdn.13$Tot.t),
                                                         mc.cores = mc.cores))
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
S3 <- round(runif(100, 1e4, 2e4))
sim.tnb.rad <- mclapply(S3, sim.rad, N = atdn.13$Tot.t, sad = "tnb",
                                       nb.fit = atdn.13$y.nb2, mc.cores = mc.cores, upper = 1e20)
sim.tnb.samp <- mclapply(sim.tnb.rad, sim.radsamp, tot.area  = atdn.13$Tot.A,
                           n.plots = atdn.13$N.plots, lmk.fit = atdn.13$lm.k,
                           nb.fit = atdn.13$y.nb2, mc.cores = mc.cores)
sim.tnb.estS <- data.frame(S = S3,
                              S.est.rnd = unlist(mclapply(sim.tnb.samp,
                                                          function(x) f3(x$rnd.samp, size = 0.018),
                                                          mc.cores = mc.cores)),
                              S.est.clump = sapply(sim.tnb.samp,
                                                 function(x) f3(x$clump.samp, size = 0.018) )
                              )
## Bias at the range of the confidence interval of the observed estimate ##
S4 <- round( runif(250, atdn.13$tovo.S$CIs[4,2], atdn.13$tovo.S$CIs[4,1]) )

sim.tnb.rad.ci <- mclapply(S4, sim.rad, N = atdn.13$Tot.t, sad = "tnb",
                                       nb.fit = atdn.13$y.nb2, mc.cores = mc.cores, upper = 1e20)
sim.tnb.samp.ci <- mclapply(sim.tnb.rad.ci, sim.radsamp, tot.area  = atdn.13$Tot.A,
                           n.plots = atdn.13$N.plots, lmk.fit = atdn.13$lm.k,
                           nb.fit = atdn.13$y.nb2, mc.cores = mc.cores)
sim.tnb.estS.ci <- data.frame(S = S4,
                              S.est.rnd = unlist(mclapply(sim.tnb.samp.ci,
                                                 function(x) f3(x$rnd.samp, size = 0.018) , mc.cores = mc.cores)),
                              S.est.clump = unlist(mclapply(sim.tnb.samp.ci,
                                                            function(x) f3(x$clump.samp, size = 0.018) , mc.cores = mc.cores))
                              )


# Storing in lists
bias13 <- list(
    ls = list(rads=sim.ls.rad, samples = sim.ls.samp, estimates = sim.ls.estS),
    ls.ci = list(rads=sim.ls.rad.ci, samples = sim.ls.samp.ci, estimates = sim.ls.estS.ci),
    tnb = list(rads=sim.tnb.rad, samples = sim.tnb.samp, estimates = sim.tnb.estS),
    tnb.ci = list(rads=sim.tnb.rad.ci, samples = sim.tnb.samp.ci, estimates = sim.tnb.estS.ci)
    )

save(bias13, file = "bias_ls_tnb_2013.RData")
    
