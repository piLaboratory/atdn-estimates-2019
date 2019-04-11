library(parallel)
source("../functions.R")
load("../lists_with_all_objects.RData")
##load("../abcSummaries.RData")
load("../bias_ls_tnb_2013.RData")
load("../bias_ls_tnb_2013t.RData")
load("../bias_ls_tnb_2019.RData")


## Number of cores to use
mc.cores <- 3

## Auxiliary functions ##
f3 <- function(x, obj){
    y <-  with(obj,
               try(sim.abc(rad = x, 
                       tot.area = Tot.A,
                       n.plots = N.plots,
                       lm.sd.fit = lm.sd,
                       lmk.fit = lm.k,
                       nrep=1,
                       summary = FALSE))
               )
    if(class(y)=="try-error")
        return(list(rnd.samp = NA,
                    clump.samp = NA))
    else
        return( list(rnd.samp = y$rnd.samp$with.est.error,
               clump.samp = y$clump.samp$with.est.error) )
}

f4 <- function(x, obj, index, ...){
    y <- with(obj,
         try(ulrich(x[[index]], lm.sd.fit = lm.sd , boot = TRUE, ...)$S[1,2])
         )
    if(class(y)=="try-error")
        return(NA)
    else
        return(y)
    }

################################################################################
## 2013
################################################################################
## Logseries rad
sim.ls.rad <- bias13$ls$rads
sim.ls.samp <- mclapply(sim.ls.rad, f3, obj = atdn.13, mc.cores=mc.cores)
sim.ls.estS <- data.frame(S = sapply(sim.ls.rad, length),
                         S.est.rnd = unlist(mclapply(sim.ls.samp, f4, obj = atdn.13, index = 1, mc.cores = mc.cores)),
                         S.est.clump = unlist(mclapply(sim.ls.samp, f4, obj = atdn.13, index = 2, mc.cores = mc.cores))
                         )
## TNB rad
sim.tnb.rad <- bias13$tnb$rads
sim.tnb.samp <- mclapply(sim.tnb.rad, f3, obj = atdn.13, mc.cores=mc.cores)
sim.tnb.estS <- data.frame(S = sapply(sim.tnb.rad, length),
                         S.est.rnd = unlist(mclapply(sim.tnb.samp, f4, obj = atdn.13, index = 1, mc.cores = mc.cores)),
                         S.est.clump = unlist(mclapply(sim.tnb.samp, f4, obj = atdn.13, index = 2, mc.cores = mc.cores))
                         )

## Stores in a list
bias.lse.13 <- list(
    ls = list(rads = sim.ls.rad, samples = sim.ls.samp, estimates = sim.ls.estS),
    tnb = list(rads = sim.tnb.rad, samples = sim.tnb.samp, estimates = sim.tnb.estS)
)

################################################################################
## 2013 with updated taxonomy
################################################################################
## Logseries rad
sim.ls.rad <- bias13t$ls$rads
sim.ls.samp <- mclapply(sim.ls.rad, f3, obj = atdn.13.tax, mc.cores=mc.cores)
sim.ls.estS <- data.frame(S = sapply(sim.ls.rad, length),
                          S.est.rnd = unlist(mclapply(sim.ls.samp, f4, obj = atdn.13.tax, index = 1, mc.cores = mc.cores)),
                          S.est.clump = unlist(mclapply(sim.ls.samp, f4, obj = atdn.13.tax, index = 2, mc.cores = mc.cores))
                          )
## TNB rad
sim.tnb.rad <- bias13t$tnb$rads
sim.tnb.samp <- mclapply(sim.tnb.rad, f3, obj = atdn.13.tax, mc.cores=mc.cores)
sim.tnb.estS <- data.frame(S = sapply(sim.tnb.rad, length),
                           S.est.rnd = unlist(mclapply(sim.tnb.samp, f4, obj = atdn.13.tax, index = 1, mc.cores = mc.cores)),
                           S.est.clump = unlist(mclapply(sim.tnb.samp, f4, obj = atdn.13.tax, index = 2, mc.cores = mc.cores))
                           )

## Stores in a list
bias.lse.13t <- list(
    ls = list(rads = sim.ls.rad, samples = sim.ls.samp, estimates = sim.ls.estS),
    tnb = list(rads = sim.tnb.rad, samples = sim.tnb.samp, estimates = sim.tnb.estS)
)

################################################################################
## 2019
################################################################################
## Logseries rad
sim.ls.rad <- bias19$ls$rads
sim.ls.samp <- mclapply(sim.ls.rad, f3, obj = atdn.19, mc.cores=mc.cores)
sim.ls.estS <- data.frame(S = sapply(sim.ls.rad, length),
                         S.est.rnd = unlist(mclapply(sim.ls.samp, f4, obj = atdn.19, index = 1, mc.cores = mc.cores)),
                         S.est.clump = unlist(mclapply(sim.ls.samp, f4, obj = atdn.19, index = 2, mc.cores = mc.cores))
                         )
## TNB rad
sim.tnb.rad <- bias19$tnb$rads
sim.tnb.samp <- mclapply(sim.tnb.rad, f3, obj = atdn.19, mc.cores=mc.cores)
sim.tnb.estS <- data.frame(S = sapply(sim.tnb.rad, length),
                         S.est.rnd = unlist(mclapply(sim.tnb.samp, f4, obj = atdn.19, index = 1, mc.cores = mc.cores)),
                         S.est.clump = unlist(mclapply(sim.tnb.samp, f4, obj = atdn.19, index = 2, mc.cores = mc.cores))
                         )

## Stores in a list
bias.lse.19 <- list(
    ls = list(rads = sim.ls.rad, samples = sim.ls.samp, estimates = sim.ls.estS),
    tnb = list(rads = sim.tnb.rad, samples = sim.tnb.samp, estimates = sim.tnb.estS)
)

## Saves all objects
save(bias.lse.13, bias.lse.13t, bias.lse.19, file = "bias_LSE.RData")
