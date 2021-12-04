library(parallel)
source("../functions.R")
load("../lists_with_all_objects.RData")


## Number of cores to use
mc.cores <- 12

################################################################################
## Bias of logseries method
################################################################################
## Checking the relationship between Estimated and real S
S1 <- round(runif(250, 1e4, 2.5e4))

sim.ls.rad <- mclapply(S1, sim.rad, N = atdn.19$Tot.t, sad = "ls", 
                          mc.cores = mc.cores, upper = 1e12)
sim.ls.samp <- mclapply(sim.ls.rad, sim.radsamp, tot.area  = atdn.19$Tot.A,
                           n.plots = atdn.19$N.plots, lmk.fit = atdn.19$lm.k,
                           mc.cores = mc.cores)
sim.ls.estS <- data.frame(S = S1,
                             S.est.rnd = unlist( mclapply(sim.ls.samp,
                                                         function(x) ls.estS(x$rnd.samp, N = atdn.19$Tot.t),
                                                         mc.cores = mc.cores)),
                             S.est.clump = unlist( mclapply(sim.ls.samp,
                                                         function(x) ls.estS(x$clump.samp, N = atdn.19$Tot.t),
                                                         mc.cores = mc.cores))
                             )

################################################################################
## Bias of negative binomial method
################################################################################
f3 <- function(x, size, loglink=TRUE){
    x <- sort(x[x>0])
    result <- try(tovo(fitnbinom2(x, start.value=c(size=size, mu=mean(x))),
                       p=sum(x)/atdn.19$Tot.t, loglink=loglink))
    if(class(result)=="try-error")
        return(NA)
    else
        return(result)
}

## Bias over a wide range species richness
S3 <- round(runif(100, 5e3, 1.85e4))
sim.tnb.rad <- mclapply(S3, sim.rad, N = atdn.19$Tot.t, sad = "tnb",
                                       nb.fit = atdn.19$y.nb2, mc.cores = mc.cores, upper = 1e50)
sim.tnb.samp <- mclapply(sim.tnb.rad, sim.radsamp, tot.area  = atdn.19$Tot.A,
                           n.plots = atdn.19$N.plots, lmk.fit = atdn.19$lm.k,
                           mc.cores = mc.cores)
sim.tnb.estS <- data.frame(S = S3,
                              S.est.rnd = unlist(mclapply(sim.tnb.samp,
                                                          function(x) f3(x$rnd.samp, size = 0.018),
                                                          mc.cores = mc.cores)),
                              S.est.clump = sapply(sim.tnb.samp,
                                                 function(x) f3(x$clump.samp, size = 0.018) )
                              )

# Storing in lists
bias19 <- list(
    ls = list(rads=sim.ls.rad, samples = sim.ls.samp, estimates = sim.ls.estS),
    tnb= list(rads=sim.tnb.rad, samples = sim.tnb.samp, estimates = sim.tnb.estS)
    )

save(bias19, file = "bias_ls_tnb_2019.RData")
    
