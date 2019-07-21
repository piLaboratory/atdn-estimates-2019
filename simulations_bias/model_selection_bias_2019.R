library(parallel)
source("../functions.R")
load("../lists_with_all_objects.RData")
load("../bias_ls_tnb_2019.RData")

## Number of cores to use
mc.cores <- 5

## Auxiliary functions
f10 <- function(x, size, index, ...){
    y <- x[[index]]
    x <- sort(y[y>0])
    ## TNB
    tnb <- try(fitnbinom2(x, start.value=c(size=size, mu=mean(x)), ...))
    if(class(tnb)=="try-error"){
        AIC.tnb <- NA
        BIC.tnb <- NA
    }
    else{
        AIC.tnb <- AIC(tnb)
        BIC.tnb <- BIC(tnb)
         } 
    ## LS
    ls <- try(fitls(x))
    if(class(ls)=="try-error"){
        AIC.ls <- NA
        BIC.ls <- NA
    }
    else{
        AIC.ls <- AIC(ls)
        BIC.ls <- BIC(ls)
    }
    return(data.frame(AIC.ls=AIC.ls, AIC.tnb=AIC.tnb, BIC.ls=BIC.ls, BIC.tnb=BIC.tnb))
}

ic.w2 <- function(x){
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
ics.rnd <- mclapply(sim.ls.samp, f10, size = 0.18, index = 1, mc.cores = mc.cores)
ics.clump <- mclapply(sim.ls.samp, f10, size = 0.18, index = 2, mc.cores = mc.cores)
sim.ls.msel <- data.frame(S = bias19$ls$estimates$S,
                          ls.rnd.aic = sapply(ics.rnd, function(x) unlist(x$AIC.ls)),
                          tnb.rnd.aic = sapply(ics.rnd, function(x) unlist(x$AIC.tnb)),
                          ls.clump.aic = sapply(ics.clump, function(x) unlist(x$AIC.ls)),
                          tnb.clump.aic = sapply(ics.clump, function(x) unlist(x$AIC.tnb)),
                          ls.rnd.bic = sapply(ics.rnd, function(x) unlist(x$BIC.ls)),
                          tnb.rnd.bic = sapply(ics.rnd, function(x) unlist(x$BIC.tnb)),
                          ls.clump.bic = sapply(ics.clump, function(x) unlist(x$BIC.ls)),
                          tnb.clump.bic = sapply(ics.clump, function(x) unlist(x$BIC.tnb))
                          )
sim.ls.msel$ls.rnd.w.aic <- apply(sim.ls.msel[,c("ls.rnd.aic","tnb.rnd.aic")], 1, function(x) ic.w2(x)[1])
sim.ls.msel$ls.clump.w.aic <- apply(sim.ls.msel[,c("ls.clump.aic","tnb.clump.aic")], 1, function(x) ic.w2(x)[1])
sim.ls.msel$ls.rnd.w.bic <- apply(sim.ls.msel[,c("ls.rnd.bic","tnb.rnd.bic")], 1, function(x) ic.w2(x)[1])
sim.ls.msel$ls.clump.w.bic <- apply(sim.ls.msel[,c("ls.clump.bic","tnb.clump.bic")], 1, function(x) ic.w2(x)[1])


################################################################################
## Model selection AICs for samples of a TNB community
################################################################################
## Simulated samples (already done by bias analyses of S estimates from LS and TNB)
sim.tnb.samp <- bias19$tnb$samples
## AICs of each model for each sample, with and without clumping
ics.rnd <- mclapply(sim.tnb.samp, f10, size = 0.18, index = 1, mc.cores = mc.cores, upper=1e40)
ics.clump <- mclapply(sim.tnb.samp, f10, size = 0.18, index = 2, mc.cores = mc.cores, upper = 1e40)
sim.tnb.msel <- data.frame(S = bias19$tnb$estimates$S,
                          ls.rnd.aic = sapply(ics.rnd, function(x) unlist(x$AIC.ls)),
                          tnb.rnd.aic = sapply(ics.rnd, function(x) unlist(x$AIC.tnb)),
                          ls.clump.aic = sapply(ics.clump, function(x) unlist(x$AIC.ls)),
                          tnb.clump.aic = sapply(ics.clump, function(x) unlist(x$AIC.tnb)),
                          ls.rnd.bic = sapply(ics.rnd, function(x) unlist(x$BIC.ls)),
                          tnb.rnd.bic = sapply(ics.rnd, function(x) unlist(x$BIC.tnb)),
                          ls.clump.bic = sapply(ics.clump, function(x) unlist(x$BIC.ls)),
                          tnb.clump.bic = sapply(ics.clump, function(x) unlist(x$BIC.tnb))
                          )
sim.tnb.msel$ls.rnd.w.aic <- apply(sim.tnb.msel[,c("ls.rnd.aic","tnb.rnd.aic")], 1, function(x) ic.w2(x)[1])
sim.tnb.msel$ls.clump.w.aic <- apply(sim.tnb.msel[,c("ls.clump.aic","tnb.clump.aic")], 1, function(x) ic.w2(x)[1])
sim.tnb.msel$ls.rnd.w.bic <- apply(sim.tnb.msel[,c("ls.rnd.bic","tnb.rnd.bic")], 1, function(x) ic.w2(x)[1])
sim.tnb.msel$ls.clump.w.bic <- apply(sim.tnb.msel[,c("ls.clump.bic","tnb.clump.bic")], 1, function(x) ic.w2(x)[1])


# Storing in a list
bias.msel.19 <- list(
    ls = sim.ls.msel,
    tnb= sim.tnb.msel
    )

save(bias.msel.19, file = "../bias_msel.RData")
    
