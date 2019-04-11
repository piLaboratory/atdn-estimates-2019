library(parallel)
library(vegan)
source("../functions.R")
load("../lists_with_all_objects.RData")
load("../bias_ls_tnb_2013.RData")
load("../bias_ls_tnb_2013t.RData")
load("../bias_ls_tnb_2019.RData")

## Auxiliary function
bias.chao <- function(obj){
    ## Chao S estimates for communities sampled under random and clumped distribution of conspecifics
    ## for a Logseries coomunity ##
    ls.rnd <- data.frame(t(sapply(obj$ls$samples, function(x) estimateR(x$rnd.samp[x$rnd.samp>0]))))
    ls.clump <- data.frame(t(sapply(obj$ls$samples, function(x) estimateR(x$clump.samp[x$clump.samp>0]))))
    ls.rnd$S <- obj$ls$estimates$S
    ls.clump$S <- obj$ls$estimates$S
    ls.estimates <- data.frame(S = obj$ls$estimates$S, S.est.rnd = ls.rnd$S.chao1, S.est.clump=ls.clump$S.chao1)
    ## for a TNB community ##
    tnb.rnd <- data.frame(t(sapply(obj$tnb$samples, function(x) estimateR(x$rnd.samp[x$rnd.samp>0]))))
    tnb.clump <- data.frame(t(sapply(obj$tnb$samples, function(x) estimateR(x$clump.samp[x$clump.samp>0]))))
    tnb.rnd$S <- obj$tnb$estimates$S
    tnb.clump$S <- obj$tnb$estimates$S
    tnb.estimates <- data.frame(S = obj$tnb$estimates$S, S.est.rnd = tnb.rnd$S.chao1, S.est.clump=tnb.clump$S.chao1)
    ## Storing in a list
    list(
        ls = list(rnd = ls.rnd, clump = ls.clump, samples = obj$ls$samples, estimates = ls.estimates),
        tnb = list(rnd = tnb.rnd, clump = tnb.clump, samples = obj$tnb$samples, estimates = tnb.estimates) 
    )
}

################################################################################
## 2013 ##
bias.chao.13 <- bias.chao(bias13)

## 2013 updated ##
bias.chao.13t <- bias.chao(bias13t)

## 2019 ##
bias.chao.19 <- bias.chao(bias19)

################################################################################
## Save objects in a binary file
################################################################################
save(bias.chao.13, bias.chao.13t, bias.chao.19, file = "../bias_chao.RData")
    
