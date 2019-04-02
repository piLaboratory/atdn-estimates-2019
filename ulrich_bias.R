library(parallel)
source("functions.R")
load("lists_with_all_objects.RData")
load("abcSummaries.RData")

f1 <- function(x, obj,...){
    with(obj,
         sim.abc(x, 
                 N = Tot.t,
                 sad = "ls",
                 tot.area = Tot.A,
                 n.plots = N.plots,
                 lm.sd.fit = lm.sd,
                 lmk.fit = lm.k,
                 nb.fit =y.nb2,
                 summary = FALSE, ...)
         )
}

f2 <- function(x, i=2) ulrich(x$clump.samp[,i])$S[1,1]

## Species richness to run the simulation
S1 <- runif(100, 12000, 20000)
### 2013
ls.samp.13 <- mclapply(S1, f1, obj = atdn.13, upper=1e16, mc.cores = 3)
ls.samp.13.u <- sapply(ls.samp.13, f2)
plot(ls.samp.13.u, S1)
abline(0,1)
ulr13.b.lm <- lm(S ~ S.est, data = data.frame(S = S1, S.est=ls.samp.13.u))
predict(ulr13.b.lm, data.frame(S.est = atdn.13$S.ulrich$S[1,1]), interval = "confidence")

## 2013 updated taxonomy
### 2013
ls.samp.13t <- mclapply(S1, f1, obj = atdn.13.tax, upper=1e16, mc.cores = 3)
ls.samp.13t.u <- sapply(ls.samp.13t, f2)
plot(ls.samp.13t.u, S1)
abline(0,1)
ulr13t.b.lm <- lm(S ~ S.est, data = data.frame(S = S1, S.est=ls.samp.13t.u))
predict(ulr13t.b.lm, data.frame(S.est = atdn.13.tax$S.ulrich$S[1,1]), interval = "confidence")

### 2019
ls.samp.19 <- mclapply(S1, f1, obj = atdn.19, upper=1e16, mc.cores = 3)
ls.samp.19.u <- sapply(ls.samp.19, f2)
plot(ls.samp.19.u, S1)
abline(0,1)
ulr19.b.lm <- lm(S ~ S.est, data = data.frame(S = S1, S.est=ls.samp.19.u))
predict(ulr19.b.lm, data.frame(S.est = atdn.19$S.ulrich$S[1,1]), interval = "confidence")
