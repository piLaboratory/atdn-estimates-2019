library(parallel)
source("../functions.R")
load("../lists_with_all_objects.RData")
load("bias_ls_tnb_2019.RData")

## Number of cores to use
mc.cores <- 3

## Auxiliary functions
f10 <- function(x, size, index, ...){
    y <- x[[index]]
    x <- sort(y[y>0])
    ## TNB
    tnb <- try(fitnbinom2(x, start.value=c(size=size, mu=mean(x)), ...))
    if(class(tnb)=="try-error")
        AIC.tnb <- NA
    else
        AIC.tnb <- AIC(tnb)
    ## LS
    ls <- try(fitls(x))
    if(class(ls)=="try-error")
        AIC.ls <- NA
    else
        AIC.ls <- AIC(ls)
    return(data.frame(AIC.ls=AIC.ls, AIC.tnb=AIC.tnb))
}

aic.w2 <- function(x){
    x <- as.vector(x)
    min.aic <- min(x)
    daics <- sapply(x, function(x) x- min.aic)
    ws <- exp(-0.5*daics)
    ws/sum(ws)
    }
################################################################################
## Model selection AICs for samples of a logseries community
################################################################################
## Simulated samples (already done by bias analyses of S estimates from LS and TNB)
sim.ls.samp <- bias19$ls$samples
## AICs of each model for each sample, with and without clumping
aics.rnd <- mclapply(sim.ls.samp, f10, size = 0.18, index = 1, mc.cores = mc.cores)
aics.clump <- mclapply(sim.ls.samp, f10, size = 0.18, index = 2, mc.cores = mc.cores)
sim.ls.msel <- data.frame(S = bias19$ls$estimates$S,
                          ls.rnd = sapply(aics.rnd, function(x) unlist(x$AIC.ls)),
                          tnb.rnd = sapply(aics.rnd, function(x) unlist(x$AIC.tnb)),
                          ls.clump = sapply(aics.clump, function(x) unlist(x$AIC.ls)),
                          tnb.clump = sapply(aics.clump, function(x) unlist(x$AIC.tnb))
                          )
sim.ls.msel$ls.rnd.w <- apply(sim.ls.msel[,c("ls.rnd","tnb.rnd")], 1, function(x) aic.w2(x)[1])
sim.ls.msel$ls.clump.w <- apply(sim.ls.msel[,c("ls.clump","tnb.clump")], 1, function(x) aic.w2(x)[1])


################################################################################
## Model selection AICs for samples of a TNB community
################################################################################
## Simulated samples (already done by bias analyses of S estimates from LS and TNB)
sim.tnb.samp <- bias19$tnb$samples
## AICs of each model for each sample, with and without clumping
aics.rnd <- mclapply(sim.tnb.samp, f10, size = 0.18, index = 1, mc.cores = mc.cores, upper=1e40)
aics.clump <- mclapply(sim.tnb.samp, f10, size = 0.18, index = 2, mc.cores = mc.cores, upper = 1e40)
sim.tnb.msel <- data.frame(S = bias19$tnb$estimates$S,
                          ls.rnd = sapply(aics.rnd, function(x) unlist(x$AIC.ls)),
                          tnb.rnd = sapply(aics.rnd, function(x) unlist(x$AIC.tnb)),
                          ls.clump = sapply(aics.clump, function(x) unlist(x$AIC.ls)),
                          tnb.clump = sapply(aics.clump, function(x) unlist(x$AIC.tnb))
                          )
sim.tnb.msel$ls.rnd.w <- apply(sim.tnb.msel[,c("ls.rnd","tnb.rnd")], 1, function(x) aic.w2(x)[1])
sim.tnb.msel$ls.clump.w <- apply(sim.tnb.msel[,c("ls.clump","tnb.clump")], 1, function(x) aic.w2(x)[1])


# Storing in a list
bias.msel.19 <- list(
    ls = sim.ls.msel,
    tnb= sim.tnb.msel
    )

save(bias.msel.19, file = "bias_msel.RData")
    
