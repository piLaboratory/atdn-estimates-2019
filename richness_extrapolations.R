## Expected number of species from simulated samples of
## regional SADs with the estimated species richness and upper and lower bounds
source("functions.R")
library(dplyr)
library(tidyr)
## Basic quantities and original (biased estimates), from script 'dataprep.R'
load("lists_with_all_objects.RData")
## Auxiliary function
f.sp <- function(n.plots, ...){
    df <- sp.samp(n.plots=n.plots, ...)
    c(rnd.mean = mean(df$rnd.samp),clump.mean = mean(df$clump.samp))
}

## Estimated values: mean of estimatimates by differente methods, weighted by the inverse of thei standard errors
## Lower and upper range: minimum and maximum of the range estimates of the different methods
## Methods used: LS, LSE, CHAO (all of these bias-corrected) and ABC (no bias correction necessary)
S1  <- 
    read.csv("figs_and_tables/estimates_S_table.csv") %>%
    filter(bias.corrected==TRUE&type!="LSE TNB"&sampling=="clump"&type!="TNB"&dataset=="2019") %>%
    mutate(se = (IC.up-IC.low)/4)%>%
    summarise(w.mean = weighted.mean(mean, w=1/se), min=min(IC.low), max=max(IC.up)) %>%    
    unlist()
## Simulated regional SAD
ls.m <- rad.ls(S = S1[1], N = atdn.19$Tot.t)
ls.low <- rad.ls(S = S1[2], N = atdn.19$Tot.t)
ls.up <- rad.ls(S = S1[3], N = atdn.19$Tot.t)
## Simulated samples
n.plots <- ceiling( seq(round(atdn.13.tax$N.plots*.9), 2*atdn.19$N.plots, length=100) )
S.est.mean <- mclapply(n.plots, f.sp,  rad = ls.m$y,
     tot.area=atdn.19$Tot.A,
     lmk.fit = atdn.19$lm.k,
     nb.fit = atdn.19$y.nb2,
     mc.cores=5)
S.est.low <- mclapply(n.plots, f.sp,  rad = ls.low$y,
     tot.area=atdn.19$Tot.A,
     lmk.fit = atdn.19$lm.k,
     nb.fit = atdn.19$y.nb2,
     mc.cores=5)
S.est.upp <- mclapply(n.plots, f.sp,  rad = ls.up$y,
     tot.area=atdn.19$Tot.A,
     lmk.fit = atdn.19$lm.k,
     nb.fit = atdn.19$y.nb2,
     mc.cores=5)
S.proj.19 <- list(n.plots=n.plots, S.est.mean=S.est.mean, S.est.low=S.est.low, S.est.upp=S.est.upp)
save(S.proj.19, file="richness_extrapolation_19.RData")
