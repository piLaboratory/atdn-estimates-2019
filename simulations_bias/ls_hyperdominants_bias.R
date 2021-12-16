library(parallel)
library(tidyverse)
source("../functions.R")
load("../lists_with_all_objects.RData")


## Number of cores to use
mc.cores <- 6

################################################################################
## Bias of logseries method
################################################################################
## Checking the relationship between Estimated and real S
## a set of 250 values of species richness that encompass the range of estimated species richness
S1 <- round(runif(250, 1e4, 2.5e4))
## Simulates a regional log-series RAD for each species richness,
## assuming the observed total regional abundance (taht is, estimated
## total number of trees in the Amazon basin)
sim.ls.rad <- mclapply(S1, sim.rad, N = atdn.13$Tot.t, sad = "ls", 
                       mc.cores = mc.cores, upper = 1e12)
## Simulates samples of trees from the RADs created above. Each sample is
## composed of the observed number of plots, and is performed with and
## without conspecific aggregation over the plots.
sim.ls.samp <- mclapply(sim.ls.rad, sim.radsamp, tot.area  = atdn.13$Tot.A,
                        n.plots = atdn.13$N.plots, lmk.fit = atdn.13$lm.k,
                        mc.cores = mc.cores)
## Estimates species richness form each sample using the logseries
## model.  The resulting data.frame shows the true richness in the
## community and the richness estimated from samples with
## (S.est.clump) and without (S.est.rnd) conspecific aggregation.  As
## conspecif aggregation cause a contant bia on the richness
## estimates, this data can the be used to estimate this bias and
## apply it to the estimated richness from the empirical data with
## logseries.
sim.ls.estS <- data.frame(S = S1,
                             S.est.rnd = unlist( mclapply(sim.ls.samp,
                                                         function(x) ls.estS(x$rnd.samp, N = atdn.13$Tot.t),
                                                         mc.cores = mc.cores)),
                             S.est.clump = unlist( mclapply(sim.ls.samp,
                                                         function(x) ls.estS(x$clump.samp, N = atdn.13$Tot.t),
                                                         mc.cores = mc.cores))
                             )
################################################################################
## Hyperdominants
################################################################################
## Counts the number of hyperdominants in the sampled LS rad
## As these RADS have already been generated, it is simple as that:
hyper.true <- sapply(sim.ls.rad, nDom, p=0.5)

## Estimates the number of hyperdominants from the Richness estimated form the sample
## First we simulate RADs assuming the richness estimated by the sample and
## the total number of trees in the metacommunity
sim.ls.rad2 <- mclapply(sim.ls.estS$S.est.clump, sim.rad, N = atdn.13$Tot.t, sad = "ls", 
                        mc.cores = mc.cores, upper = 1e12)
## And the count the number of hyperdominants in the RADs simulated above
hyper.est.clump <- sapply(sim.ls.rad2, nDom, p=0.5)

## And here's the plot of proportion of dominants in the true RAD and in the sample with clumpoing
## A data.frame to store results
sim.ls.estHyp <-
    cbind(sim.ls.estS, H = hyper.true, H.clump = hyper.est.clump) %>%
    mutate(H.prop = H/S, H.clump.prop = H.clump / S.est.clump)
## A plot of the true proportion od hyperdominants x propostions got from the estimated S in a sampling with clumping
plot(H.prop ~ H.clump.prop, data = sim.ls.estHyp)
abline(0,1)
## A nice linear relantionship. Proportion od hyperdominants estimated
## from samples with conspecific aggregation underestimates the true proportion.
