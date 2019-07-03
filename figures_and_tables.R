source("functions.R")
library(dplyr)
library(tidyr)
## Basic quantities and original (biased estimates), from script 'dataprep.R'
load("lists_with_all_objects.RData")
## Summaries and posterior distributions of ABC analyses, from script 'abc_summaries.R'
load("abcSummaries.RData")
## Bias analyses for estimates from logseries and tnb fitted to the sample for 2013 dataset , from script
## 'simulation_bias/ls_tnb_bias_parallel_2013.R'
load("bias_ls_tnb_2013.RData")
## same for 2013 updated
load("bias_ls_tnb_2013t.RData")
## same for 2019
load("bias_ls_tnb_2019.RData")
## Bias analyses for model selection of ls, tnb, pln fitted to the samples, from script 'simulation_bias/model_selection_bias.R'
load("bias_msel.RData")
## Bias analyses for estimates from CHAO method from script 'simulation_bias/chao_bias.R'
load("bias_chao.RData")
## Bias analyses for estimates from LSE method on the estimated population sizes,
## from script 'simulation_bias/LSE_bias.R'
load("bias_LSE.RData")
## Extrapolated species richness from simulated samples with clumping from a logseries
## from script 'richness_extrapolations.R'
load("richness_extrapolation_19.RData")

################################################################################
## Model selection table for ls, tnb, poilog models for each dataset
################################################################################
ms13 <- (with(atdn.13, AICtab(y.nb2, y.ls, pln, weights = TRUE, sort=FALSE)))
ms13t <- (with(atdn.13.tax, AICtab(y.nb2, y.ls, pln, weights = TRUE, sort=FALSE)))
ms19 <- (with(atdn.19, AICtab(y.nb2, y.ls, pln, weights = TRUE, sort=FALSE)))
tab.ms.all <- data.frame(
    model = c("TNB", "LS", "PLN"),
    df = ms13$df,
    y.13 = ms13$dAIC,
    y.13t = ms13t$dAIC,
    y.19 = ms19$dAIC
)
write.csv(tab.ms.all, row.names=FALSE, file = "figs_and_tables/model_selection_table.csv")

################################################################################
## Table with estimated coefficients and SEs and likelihoods for ls, tnb, poilog models for each dataset
################################################################################
tab.cf <- rbind(
    cbind(with(atdn.13, summary(y.nb2)@coef[,1:2]),
          with(atdn.13.tax, summary(y.nb2)@coef[,1:2]),
          with(atdn.19, summary(y.nb2)@coef[,1:2])),
    cbind(with(atdn.13, matrix(summary(y.ls)@coef[,1:2], nrow=1)),
          with(atdn.13.tax, matrix(summary(y.ls)@coef[,1:2], nrow=1)),
          with(atdn.19, matrix(summary(y.ls)@coef[,1:2], nrow=1))),
    cbind(with(atdn.13, summary(pln)@coef[,1:2]),
          with(atdn.13.tax, summary(pln)@coef[,1:2]),
          with(atdn.19, summary(pln)@coef[,1:2]))
    )
colnames(tab.cf) <- paste(colnames(tab.cf), rep(c("2013","2013t","2019"),each=2))
rownames(tab.cf)[rownames(tab.cf)==""] <- "alpha"
write.csv(tab.cf, file = "figs_and_tables/model_coefficients.csv")

################################################################################
## ABC posterior probabilities for each model for each dataset
################################################################################
tab.abc <- rbind(
    summary(abc2013.summ$model.sel)$Prob,
    summary(abc2013t.summ$model.sel)$Prob,
    summary(abc2019.summ$model.sel)$Prob
    )
rownames(tab.abc) <- c("2013","2013t","2019")
write.csv(tab.abc, row.names=FALSE, file = "figs_and_tables/abc_post_table.csv")


################################################################################
## Tables with all estimates of species richness and CI's 95%
################################################################################
## Not exactly an elegant code, but alas, it works

## Table with original estimates and CI's (without bias correction) ##
S.estimates <- expand.grid(
    dataset = c("2013", "2013 updated", "2019"),
    type = c("CHAO", "TNB", "LS", "LSE"),
    mean = NA,
    IC.low = NA,
    IC.up = NA
)
## Add a missing column of sampling in simulations (to compatibility with the next table)
S.estimates$sampling <- NA
S.estimates <- S.estimates[,c(1:2,6,3:5)]
## Including the values
## Truncated negative binomial
S.estimates[S.estimates$type=="TNB"&S.estimates$dataset=="2013",4:6] <-
    with(atdn.13, c(tovo.S$S.est, tovo.S$CIs[4,2:1]))
S.estimates[S.estimates$type=="TNB"&S.estimates$dataset=="2013 updated",4:6] <-
    with(atdn.13.tax, c(tovo.S$S.est, tovo.S$CIs[4,2:1]))
S.estimates[S.estimates$type=="TNB"&S.estimates$dataset=="2019",4:6] <-
    with(atdn.19, c(tovo.S$S.est, tovo.S$CIs[4,2:1]))
## Log-series
S.estimates[S.estimates$type=="LS"&S.estimates$dataset=="2013",4:6] <-
    with(atdn.13, c(S.ls, S.ls.ci))
S.estimates[S.estimates$type=="LS"&S.estimates$dataset=="2013 updated",4:6] <-
    with(atdn.13.tax, c(S.ls, S.ls.ci))
S.estimates[S.estimates$type=="LS"&S.estimates$dataset=="2019",4:6] <-
    with(atdn.19, c(S.ls, S.ls.ci))
## Linear extension of RAD of estimated pop sizes (LSE)
S.estimates[S.estimates$type=="LSE"&S.estimates$dataset=="2013",4:6] <-
    with(atdn.13, S.ulrich$S[1,2:4])
S.estimates[S.estimates$type=="LSE"&S.estimates$dataset=="2013 updated",4:6] <-
    with(atdn.13.tax, S.ulrich$S[1,2:4])
S.estimates[S.estimates$type=="LSE"&S.estimates$dataset=="2019",4:6] <-
    with(atdn.19, S.ulrich$S[1,2:4])
## Chao 1
S.estimates[S.estimates$type=="CHAO"&S.estimates$dataset=="2013",4:6] <-
    with(atdn.13, Chao[c(2,4:5)])
S.estimates[S.estimates$type=="CHAO"&S.estimates$dataset=="2013 updated",4:6] <-
    with(atdn.13.tax, Chao[c(2,4:5)])
S.estimates[S.estimates$type=="CHAO"&S.estimates$dataset=="2019",4:6] <-
    with(atdn.19, Chao[c(2,4:5)])

## Table of Bias-corrected estimates ##
S.estimates.bc <- expand.grid(
    sampling = c("rnd", "clump"),
    type = c("TNB", "LS", "LSE LS", "LSE TNB", "CHAO", "ABC"),
    dataset = c("2013", "2013 updated", "2019"),    
    mean = NA,
    IC.low = NA,
    IC.up = NA,
    stringsAsFactors=FALSE
)[c(3:1,4:6)]
## Including the values
## Truncated negative binomial
type <- "TNB"
## 2013
dataset <- "2013"
obj <- bias13$tnb
S.estimates.bc[S.estimates.bc$type==type&S.estimates.bc$dataset==dataset,4:6] <-
    bias.ci(obj, ci.vector = S.estimates[S.estimates$type==type&S.estimates$dataset==dataset,4:6])
##2013 updated
dataset <- "2013 updated"
obj <- bias13t$tnb
S.estimates.bc[S.estimates.bc$type==type&S.estimates.bc$dataset==dataset,4:6] <-
    bias.ci(obj, ci.vector = S.estimates[S.estimates$type==type&S.estimates$dataset==dataset,4:6])
## 2019
dataset <- "2019"
obj <- bias19$tnb
S.estimates.bc[S.estimates.bc$type==type&S.estimates.bc$dataset==dataset,4:6] <-
    bias.ci(obj, ci.vector = S.estimates[S.estimates$type==type&S.estimates$dataset==dataset,4:6])
## Logseries
type <- "LS"
## 2013
dataset <- "2013"
obj <- bias13$ls
S.estimates.bc[S.estimates.bc$type==type&S.estimates.bc$dataset==dataset,4:6] <-
    bias.ci(obj, ci.vector = S.estimates[S.estimates$type==type&S.estimates$dataset==dataset,4:6], method="lm")
##2013 updated
dataset <- "2013 updated"
obj <- bias13t$ls
S.estimates.bc[S.estimates.bc$type==type&S.estimates.bc$dataset==dataset,4:6] <-
    bias.ci(obj, ci.vector = S.estimates[S.estimates$type==type&S.estimates$dataset==dataset,4:6], method="lm")
## 2019
dataset <- "2019"
obj <- bias19$ls
S.estimates.bc[S.estimates.bc$type==type&S.estimates.bc$dataset==dataset,4:6] <-
    bias.ci(obj, ci.vector = S.estimates[S.estimates$type==type&S.estimates$dataset==dataset,4:6], method="lm")
## LSE (linear extension of the RAD of estimated population sizes, assuming a log-series RAD)
type <- "LSE LS"
## 2013
dataset <- "2013"
obj <- bias.lse.13$ls
S.estimates.bc[S.estimates.bc$type==type&S.estimates.bc$dataset==dataset,4:6] <-
    bias.ci(obj, ci.vector = S.estimates[S.estimates$type=="LSE"&S.estimates$dataset==dataset,4:6], method="lm")
##2013 updated
dataset <- "2013 updated"
obj <- bias.lse.13t$ls
S.estimates.bc[S.estimates.bc$type==type&S.estimates.bc$dataset==dataset,4:6] <-
    bias.ci(obj, ci.vector = S.estimates[S.estimates$type=="LSE"&S.estimates$dataset==dataset,4:6], method="lm")
## 2019
dataset <- "2019"
obj <- bias.lse.19$ls
S.estimates.bc[S.estimates.bc$type==type&S.estimates.bc$dataset==dataset,4:6] <-
    bias.ci(obj, ci.vector = S.estimates[S.estimates$type=="LSE"&S.estimates$dataset==dataset,4:6], method="lm")
## LSE (linear extension of the RAD of estimated population sizes, assuming a TNB RAD)
type <- "LSE TNB"
## 2013
dataset <- "2013"
obj <- bias.lse.13$tnb
S.estimates.bc[S.estimates.bc$type==type&S.estimates.bc$dataset==dataset,4:6] <-
    bias.ci(obj, ci.vector = S.estimates[S.estimates$type=="LSE"&S.estimates$dataset==dataset,4:6], method="lm")
##2013 updated
dataset <- "2013 updated"
obj <- bias.lse.13t$tnb
S.estimates.bc[S.estimates.bc$type==type&S.estimates.bc$dataset==dataset,4:6] <-
    bias.ci(obj, ci.vector = S.estimates[S.estimates$type=="LSE"&S.estimates$dataset==dataset,4:6], method="lm")
## 2019
dataset <- "2019"
obj <- bias.lse.19$tnb
S.estimates.bc[S.estimates.bc$type==type&S.estimates.bc$dataset==dataset,4:6] <-
    bias.ci(obj, ci.vector = S.estimates[S.estimates$type=="LSE"&S.estimates$dataset==dataset,4:6], method="lm")
## CHAO
type <- "CHAO"
## 2013
dataset <- "2013"
obj <- bias.chao.13$ls
S.estimates.bc[S.estimates.bc$type==type&S.estimates.bc$dataset==dataset,4:6] <-
    bias.ci(obj, ci.vector = S.estimates[S.estimates$type=="CHAO"&S.estimates$dataset==dataset,4:6], method="lm")
##2013 updated
dataset <- "2013 updated"
obj <- bias.chao.13t$ls
S.estimates.bc[S.estimates.bc$type==type&S.estimates.bc$dataset==dataset,4:6] <-
    bias.ci(obj, ci.vector = S.estimates[S.estimates$type=="CHAO"&S.estimates$dataset==dataset,4:6], method="lm")
## 2019
dataset <- "2019"
obj <- bias.chao.19$ls
S.estimates.bc[S.estimates.bc$type==type&S.estimates.bc$dataset==dataset,4:6] <-
    bias.ci(obj, ci.vector = S.estimates[S.estimates$type=="CHAO"&S.estimates$dataset==dataset,4:6], method="lm")
## ABC (only the selected model)
S.estimates.bc[S.estimates.bc$type=="ABC"&S.estimates.bc$dataset=="2013"&S.estimates.bc$sampling=="clump",4:6] <-
    summary(abc2013.summ$S.post1)[c(4,2,6),]
S.estimates.bc[S.estimates.bc$type=="ABC"&S.estimates.bc$dataset=="2013 updated"&S.estimates.bc$sampling=="clump",4:6] <-
    summary(abc2013t.summ$S.post1)[c(4,2,6),]
S.estimates.bc[S.estimates.bc$type=="ABC"&S.estimates.bc$dataset=="2019"&S.estimates.bc$sampling=="clump",4:6] <-
    summary(abc2019.summ$S.post1)[c(4,2,6),]
## A single table with all values (for supplemmentary material) ##
S.estimates.all <- rbind(S.estimates.bc, S.estimates)
S.estimates.all$bias.corrected <- !is.na(S.estimates.all$sampling)

write.csv(S.estimates.all, row.names=FALSE, file = "figs_and_tables/estimates_S_table.csv")

################################################################################
## Figures with relationship estimated x true values of Species Richness 
################################################################################
## auxiliary ploting functions ##
## True x estimated values
bias.p1 <- function(bias.est, point, lower, upper, reg.line = TRUE, ...){
    plot(S ~ S.est.rnd, data = bias.est, type= "n", ...)
    rect(lower,0,upper,max(bias.est$S)*1.1,col="lightgrey", border=NA)
    box()
    points(S ~ S.est.rnd, data = bias.est, col="blue", cex=0.5)
    points(S ~ S.est.clump, data = bias.est, col="red", cex=0.5)
    abline(0,1)
    if(reg.line){
        lm1 <- lm(S ~ S.est.rnd, data = bias.est, subset=S.est.rnd>=lower&S.est.rnd<=upper)
        plm1 <- predict(lm1,data.frame(S.est.rnd=c(point,lower,upper)))
        lm2 <- lm(S ~ S.est.clump, data = bias.est, subset=S.est.clump>=lower&S.est.clump<=upper)
        plm2 <- predict(lm2,data.frame(S.est.clump=c(point,lower,upper)))
        lines(c(point,lower,upper), plm1, col="blue", lwd=3)
        lines(c(point,lower,upper), plm2, col="red", lwd=3)
    }
}
## Point and range bias-corrected estimates
bias.p2 <- function(dataset, type, ...){
    df <- S.estimates.all[S.estimates.all$dataset==dataset & S.estimates.all$type==type & S.estimates.all$bias.corrected,]
    df <- df[order(df$sampling),]
    plot(df$mean ~ c(1,2), type = "n",
         axes=FALSE,
         xlab="", ylab="", xlim=c(0.5,2.5), ...)
    segments(x0=1:2, x1=1:2, y0=df$IC.low, y1=df$IC.up, col=c("red", "blue"), lwd=2)
    points(1:2, df$mean, col=c("red", "blue"), cex=2)
}

### Plots for each data set ###

## 2013 ##
pdf("figs_and_tables/estSxS_2013.pdf", width = 12, height = 11)
nf <- layout(
    matrix(
        c( c(1,2:5),c(1,6:9), c(1,rep(10,4)) ),
        nrow = 3,
        ncol = 5,
        byrow = TRUE
    ),
    widths = c(0.5, rep(c(0.5,3),2)),
    heights= c(3,3,0.5))
## Y axis label
par(las = 0, mar=c(0,4,0,0), cex.axis=1.25)
plot(0,0, xlim=c(0,1), ylim =c(0,1), axes=FALSE, xlab="", ylab="", type="n")
mtext(text = "Number of species in simulated communities", side=2, cex = 1.5)
## Logseries #
par(mar = c(5, 0, 4, 0))
bias.p2("2013", "LS", ylim=atdn.13$S.ls.ci*c(.9,1.1))
par(mar = c(5, 2, 4, 4), cex.main=1.5)
p.bias(bias13$ls$estimates, atdn.13$S.ls, atdn.13$S.ls.ci[1], atdn.13$S.ls.ci[2],
       ylim = atdn.13$S.ls.ci*c(.9,1.1),
       xlim = atdn.13$S.ls.ci*c(.9,1.1),
       ylab = "", xlab = "", main = "LS")
# LSE
par(mar = c(5, 0, 4, 0))
bias.p2("2013", "LSE LS", ylim=c(1.5e4,1.95e4))
par(mar = c(5, 2, 4, 4))
p.bias(bias.lse.13$ls$estimates, atdn.13$S.ulrich$S$boot.mean[1],
       atdn.13$S.ulrich$S$boot.CI.low[1], atdn.13$S.ulrich$S$boot.CI.up[1],
       ylim = c(1.5e4,1.95e4),
       xlim = c(1.4e4,1.6e4),
       ylab = "", xlab = "", main = "LSE")
##CHAO
par(mar = c(5, 0, 4, 0))
bias.p2("2013", "CHAO", ylim=c(5e3,1.6e4))
par(mar = c(5, 2, 4, 4))
p.bias(bias.chao.13$ls$estimates, atdn.13$Chao[2], 
       atdn.13$Chao[4], atdn.13$Chao[5],
       ylim = c(5e3,1.6e4),
       xlim = atdn.13$Chao[4:5]*c(.95,1.05),
       ylab = "", xlab = "", main = "CHAO")
## TNB
par(mar = c(5, 0, 5, 0))
bias.p2("2013", "TNB", ylim=range(bias13$tnb$estimates$S))
par(mar = c(5, 2, 4, 4))
p.bias(bias13$tnb$estimates, atdn.13$tovo.S$S.est, 
       atdn.13$tovo.S$CIs[4,2], atdn.13$tovo.S$CIs[4,1],
       reg.line=FALSE,
       ylim = range(bias13$tnb$estimates$S),
       ylab = "", xlab = "", main = "TNB")
## X axis
par(mar=c(2,0,1,0))
plot(0,0, xlim=c(0,1), ylim =c(0,1), axes=FALSE, xlab="", ylab="", type="n")
mtext(text = "Estimated number of species from simulated samples", side=1, cex = 1.5)
dev.off()

## 2013 updated ##
pdf("figs_and_tables/estSxS_2013t.pdf", width = 12, height = 11)
nf <- layout(
    matrix(
        c( c(1,2:5),c(1,6:9), c(1,rep(10,4)) ),
        nrow = 3,
        ncol = 5,
        byrow = TRUE
    ),
    widths = c(0.5, rep(c(0.5,3),2)),
    heights= c(3,3,0.5))
## Y axis label
par(las = 0, mar=c(0,4,0,0), cex.axis=1.25)
plot(0,0, xlim=c(0,1), ylim =c(0,1), axes=FALSE, xlab="", ylab="", type="n")
mtext(text = "Number of species in simulated communities", side=2, cex = 1.5)
## Logseries #
par(mar = c(5, 0, 4, 0))
bias.p2("2013 updated", "LS", ylim=atdn.13.tax$S.ls.ci*c(.9,1.1))
par(mar = c(5, 2, 4, 4), cex.main=1.5)
p.bias(bias13t$ls$estimates, atdn.13.tax$S.ls, atdn.13.tax$S.ls.ci[1], atdn.13.tax$S.ls.ci[2],
       ylim = atdn.13.tax$S.ls.ci*c(.9,1.1),
       xlim = atdn.13.tax$S.ls.ci*c(.9,1.1),
       ylab = "", xlab = "", main = "LS")
# LSE
par(mar = c(5, 0, 4, 0))
bias.p2("2013 updated", "LSE LS", ylim=c(1.25e4,1.65e4))
par(mar = c(5, 2, 4, 4))
p.bias(bias.lse.13t$ls$estimates, atdn.13.tax$S.ulrich$S$boot.mean[1],
       atdn.13.tax$S.ulrich$S$boot.CI.low[1], atdn.13.tax$S.ulrich$S$boot.CI.up[1],
       ylim = c(1.25e4,1.65e4),
       xlim = c(1.3e4,1.45e4),
       ylab = "", xlab = "", main = "LSE")
##CHAO
par(mar = c(5, 0, 4, 0))
bias.p2("2013 updated", "CHAO", ylim=c(5e3,1.6e4))
par(mar = c(5, 2, 4, 4))
p.bias(bias.chao.13t$ls$estimates, atdn.13.tax$Chao[2], 
       atdn.13.tax$Chao[4], atdn.13.tax$Chao[5],
       ylim = c(5e3,1.6e4),
       xlim = atdn.13.tax$Chao[4:5]*c(.95,1.05),
       ylab = "", xlab = "", main = "CHAO")
## TNB
par(mar = c(5, 0, 5, 0))
bias.p2("2013 updated", "TNB", ylim=range(bias13t$tnb$estimates$S))
par(mar = c(5, 2, 4, 4))
p.bias(bias13t$tnb$estimates, atdn.13.tax$tovo.S$S.est, 
       atdn.13.tax$tovo.S$CIs[4,2], atdn.13.tax$tovo.S$CIs[4,1],
       reg.line=FALSE,
       ylim = range(bias13t$tnb$estimates$S),
       xlim = range(bias13$tnb$estimates$S.est.rnd),
       ylab = "", xlab = "", main = "TNB")
## X axis
par(mar=c(2,0,1,0))
plot(0,0, xlim=c(0,1), ylim =c(0,1), axes=FALSE, xlab="", ylab="", type="n")
mtext(text = "Estimated number of species from simulated samples", side=1, cex = 1.5)
dev.off()


## 2019 ##
pdf("figs_and_tables/estSxS_2019.pdf", width = 12, height = 11)
nf <- layout(
    matrix(
        c( c(1,2:5),c(1,6:9), c(1,rep(10,4)) ),
        nrow = 3,
        ncol = 5,
        byrow = TRUE
    ),
    widths = c(0.5, rep(c(0.5,3),2)),
    heights= c(3,3,0.5))
## Y axis label
par(las = 0, mar=c(0,4,0,0), cex.axis=1.25)
plot(0,0, xlim=c(0,1), ylim =c(0,1), axes=FALSE, xlab="", ylab="", type="n")
mtext(text = "Number of species in simulated communities", side=2, cex = 1.5)
## Logseries #
par(mar = c(5, 0, 4, 0))
bias.p2("2019", "LS", ylim=atdn.19$S.ls.ci*c(.9,1.1))
par(mar = c(5, 2, 4, 4), cex.main=1.5)
p.bias(bias19$ls$estimates, atdn.19$S.ls, atdn.19$S.ls.ci[1], atdn.19$S.ls.ci[2],
       ylim = atdn.19$S.ls.ci*c(.9,1.1),
       xlim = atdn.19$S.ls.ci*c(.9,1.1),
       ylab = "", xlab = "", main = "LS")
# LSE
par(mar = c(5, 0, 4, 0))
bias.p2("2019", "LSE LS", ylim=c(1.4e4,1.72e4))
par(mar = c(5, 2, 4, 4))
p.bias(bias.lse.19$ls$estimates, atdn.19$S.ulrich$S$boot.mean[1],
       atdn.19$S.ulrich$S$boot.CI.low[1], atdn.19$S.ulrich$S$boot.CI.up[1],
       ylim = c(1.4e4,1.72e4),
       xlim = c(1.4e4,1.55e4),
       ylab = "", xlab = "", main = "LSE")
##CHAO
par(mar = c(5, 0, 4, 0))
bias.p2("2019", "CHAO", ylim=c(5e3,1.6e4))
par(mar = c(5, 2, 4, 4))
p.bias(bias.chao.19$ls$estimates, atdn.19$Chao[2], 
       atdn.19$Chao[4], atdn.19$Chao[5],
       ylim = c(5e3,1.6e4),
       xlim = atdn.19$Chao[4:5]*c(.95,1.05),
       ylab = "", xlab = "", main = "CHAO")
## TNB
par(mar = c(5, 0, 5, 0))
bias.p2("2019", "TNB", ylim=range(bias19$tnb$estimates$S))
par(mar = c(5, 2, 4, 4))
p.bias(bias19$tnb$estimates, atdn.19$tovo.S$S.est, 
       atdn.19$tovo.S$CIs[4,2], atdn.19$tovo.S$CIs[4,1],
       reg.line=FALSE,
       ylim = range(bias19$tnb$estimates$S),
       xlim = range(bias13$tnb$estimates$S.est.rnd),
       ylab = "", xlab = "", main = "TNB")
## X axis
par(mar=c(2,0,1,0))
plot(0,0, xlim=c(0,1), ylim =c(0,1), axes=FALSE, xlab="", ylab="", type="n")
mtext(text = "Estimated number of species from simulated samples", side=1, cex = 1.5)
dev.off()

################################################################################
## Old anc complex code for a complex figure with all data sets and methods##
## pdf("figs_and_tables/estSxS_ls_tnb.pdf", width = 12, height = 12)
## nf <- layout(
##     matrix(
##         c( c(1,2:10),c(1,11:19), c(1,20:28), c(1, 29:37), c(1,rep(38,9)) ),
##         nrow = 5,
##         ncol = 10,
##         byrow = TRUE
##     ),
##     widths = c(1, rep(c(0.3,0.3,3),4)),
##     heights= c(3,3,3,3,1))
## ## Logseries #
## ## 2013
## par(las = 0, mar=c(0,4,0,0), cex.axis=1.25)
## plot(0,0, xlim=c(0,1), ylim =c(0,1), axes=FALSE, xlab="", ylab="", type="n")
## mtext(text = "Number of species in simulated communities", side=2, cex = 1.5)
## par(mar = c(3, 0, 2, 0))
## with(subset(bias13$ls$estimates, S.est.clump>atdn.13$S.ls.ci[1] & S.est.clump<atdn.13$S.ls.ci[2]),
##      boxplot(S, axes=FALSE, ylim = atdn.13$S.ls.ci*c(.9,1.1), border="red"))
## with(subset(bias13$ls$estimates, S.est.rnd>atdn.13$S.ls.ci[1] & S.est.rnd<atdn.13$S.ls.ci[2]),
##      boxplot(S, axes=FALSE, ylim = atdn.13$S.ls.ci*c(.9,1.1), border="blue"))
## par(mar = c(3, 2, 2, 4))
## plot(S ~ S.est.rnd, data = bias13$ls$estimates, ylim = atdn.13$S.ls.ci*c(.9,1.1),
##      xlim = atdn.13$S.ls.ci*c(.9,1.1), col="blue", ylab = "", xlab = "", type="n", main = "2013")
## rect(atdn.13$S.ls.ci[1],0,atdn.13$S.ls.ci[2],max(bias13$ls$estimates$S),col="lightgrey", border=NA)
## box()
## axis(1)
## points(S ~ S.est.rnd, data = bias13$ls$estimates, col="blue")
## points(S ~ S.est.clump, data = bias13$ls$estimates, col="red")
## abline(0,1)
## segments(x0 = c(atdn.13$S.ls.ci),
##        y0 = c(atdn.13$S.ls.ci),
##        x1 = c(0,0),
##        y1 = c(atdn.13$S.ls.ci),
##        lty=2)
## ## 2013 updated taxonomy
## par(mar = c(3, 0, 2, 0))
## with(subset(bias13t$ls$estimates, S.est.clump>atdn.13.tax$S.ls.ci[1] & S.est.clump<atdn.13.tax$S.ls.ci[2]),
##      boxplot(S, axes=FALSE, ylim = atdn.13.tax$S.ls.ci*c(.9,1.1), border="red"))
## with(subset(bias13t$ls$estimates, S.est.rnd>atdn.13.tax$S.ls.ci[1] & S.est.rnd<atdn.13.tax$S.ls.ci[2]),
##      boxplot(S, axes=FALSE, ylim = atdn.13.tax$S.ls.ci*c(.9,1.1), border="blue"))
## par(mar = c(3, 2, 2, 4))
## plot(S ~ S.est.rnd, data = bias13t$ls$estimates, ylim = atdn.13.tax$S.ls.ci*c(.9,1.1),
##      xlim = atdn.13.tax$S.ls.ci*c(.9,1.1), col="blue", ylab = "", xlab = "", type="n", main = "2013 Updated")
## rect(atdn.13.tax$S.ls.ci[1],0,atdn.13.tax$S.ls.ci[2],max(bias13t$ls$estimates$S),col="lightgrey", border=NA)
## box()
## axis(1)
## points(S ~ S.est.rnd, data = bias13t$ls$estimates, col="blue")
## points(S ~ S.est.clump, data = bias13t$ls$estimates, col="red")
## abline(0,1)
## segments(x0 = c(atdn.13.tax$S.ls.ci),
##        y0 = c(atdn.13.tax$S.ls.ci),
##        x1 = c(0,0),
##        y1 = c(atdn.13.tax$S.ls.ci),
##        lty=2)
## ## 2019
## par(mar = c(3, 0, 2, 0))
## with(subset(bias19$ls$estimates, S.est.clump>atdn.19$S.ls.ci[1] & S.est.clump<atdn.19$S.ls.ci[2]),
##      boxplot(S, axes=FALSE, ylim = atdn.19$S.ls.ci*c(.9,1.1), border="red"))
## with(subset(bias19$ls$estimates, S.est.rnd>atdn.19$S.ls.ci[1] & S.est.rnd<atdn.19$S.ls.ci[2]),
##      boxplot(S, axes=FALSE, ylim = atdn.19$S.ls.ci*c(.9,1.1), border="blue"))
## par(mar = c(3, 2, 2, 4))
## plot(S ~ S.est.rnd, data = bias19$ls$estimates, ylim = atdn.19$S.ls.ci*c(.9,1.1),
##      xlim = atdn.19$S.ls.ci*c(.9,1.1), col="blue", ylab = "", xlab = "", type="n", main = "2019")
## rect(atdn.19$S.ls.ci[1],0,atdn.19$S.ls.ci[2],max(bias19$ls$estimates$S),col="lightgrey", border=NA)
## box()
## axis(1)
## points(S ~ S.est.rnd, data = bias19$ls$estimates, col="blue")
## points(S ~ S.est.clump, data = bias19$ls$estimates, col="red")
## abline(0,1)
## segments(x0 = c(atdn.19$S.ls.ci),
##        y0 = c(atdn.19$S.ls.ci),
##        x1 = c(0,0),
##        y1 = c(atdn.19$S.ls.ci),
##        lty=2)
## par(las = 1)
## mtext("LS", side = 4, line = 1.1, cex = 1.1)
## par(las = 0)
## ## TNB
## par(mar = c(3, 0, 2, 0))
## with(subset(bias13$tnb$estimates, S.est.clump>atdn.13$tovo.S$CIs[4,2] & S.est.clump<atdn.13$tovo.S$CIs[4,1]),
##      boxplot(S, axes=FALSE, ylim = range(bias13$tnb$estimates$S), border="red"))
## with(subset(bias13$tnb$estimates, S.est.rnd>atdn.13$tovo.S$CIs[4,2] & S.est.rnd<atdn.13$tovo.S$CIs[4,1]),
##      boxplot(S, axes=FALSE, ylim = range(bias13$tnb$estimates$S), border="blue"))
## par(mar = c(3, 2, 2, 4))
## plot(S ~ S.est.rnd, data = bias13$tnb$estimates, ylim = range(bias13$tnb$estimates$S),
##      col="blue", ylab = "", xlab = "", type="n", main = "")
## rect(atdn.13$tovo.S$CIs[4,2],0,atdn.13$tovo.S$CIs[4,1],max(bias13$tnb$estimates$S)*1.1,col="lightgrey", border=NA)
## box()
## axis(1)
## points(S ~ S.est.rnd, data = bias13$tnb$estimates, col="blue")
## points(S ~ S.est.clump, data = bias13$tnb$estimates, col="red")
## abline(0,1)
## segments(x0 = c(atdn.13$tovo.S$CIs[4,]),
##        y0 = c(atdn.13$tovo.S$CIs[4,]),
##        x1 = c(0,0),
##        y1 = c(atdn.13$tovo.S$CIs[4,]),
##        lty=2)
## ## 2013 updated taxonomy
## par(mar = c(3, 0, 2, 0))
## with(subset(bias13t$tnb$estimates, S.est.clump>atdn.13.tax$tovo.S$CIs[4,2] & S.est.clump<atdn.13.tax$tovo.S$CIs[4,1]),
##      boxplot(S, axes=FALSE, ylim = range(bias13t$tnb$estimates$S), border="red"))
## with(subset(bias13t$tnb$estimates, S.est.rnd>atdn.13.tax$tovo.S$CIs[4,2] & S.est.rnd<atdn.13.tax$tovo.S$CIs[4,1]),
##      boxplot(S, axes=FALSE, ylim = range(bias13t$tnb$estimates$S), border="blue"))
## par(mar = c(3, 2, 2, 4))
## plot(S ~ S.est.rnd, data = bias13t$tnb$estimates, ylim = range(bias13t$tnb$estimates$S),
##      xlim = range(bias13$tnb$estimates$S.est.rnd), ## To avoid some extreme values 
##      col="blue", ylab = "", xlab = "", type="n", main = "")
## rect(atdn.13.tax$tovo.S$CIs[4,2],0,atdn.13.tax$tovo.S$CIs[4,1],max(bias13t$tnb$estimates$S),col="lightgrey", border=NA)
## box()
## axis(1)
## points(S ~ S.est.rnd, data = bias13t$tnb$estimates, col="blue")
## points(S ~ S.est.clump, data = bias13t$tnb$estimates, col="red")
## abline(0,1)
## segments(x0 = c(atdn.13.tax$tovo.S$CIs[4,]),
##        y0 = c(atdn.13.tax$tovo.S$CIs[4,]),
##        x1 = c(0,0),
##        y1 = c(atdn.13.tax$tovo.S$CIs[4,]),
##        lty=2)
## ## 2019
## par(mar = c(3, 0, 2, 0))
## with(subset(bias19$tnb$estimates, S.est.clump>atdn.19$tovo.S$CIs[4,2] & S.est.clump<atdn.19$tovo.S$CIs[4,1]),
##      boxplot(S, axes=FALSE, ylim = range(bias19$tnb$estimates$S), border="red"))
## with(subset(bias19$tnb$estimates, S.est.rnd>atdn.19$tovo.S$CIs[4,2] & S.est.rnd<atdn.19$tovo.S$CIs[4,1]),
##      boxplot(S, axes=FALSE, ylim = range(bias19$tnb$estimates$S), border="blue"))
## par(mar = c(3, 2, 2, 4))
## plot(S ~ S.est.rnd, data = bias19$tnb$estimates, ylim = range(bias19$tnb$estimates$S),
##      xlim = range(bias13$tnb$estimates$S.est.rnd),
##      col="blue", ylab = "", xlab = "", type="n", main = "")
## rect(atdn.19$tovo.S$CIs[4,2],0,atdn.19$tovo.S$CIs[4,1],max(bias19$tnb$estimates$S),col="lightgrey", border=NA)
## box()
## axis(1)
## points(S ~ S.est.rnd, data = bias19$tnb$estimates, col="blue")
## points(S ~ S.est.clump, data = bias19$tnb$estimates, col="red")
## abline(0,1)
## segments(x0 = c(atdn.19$tovo.S$CIs[4,]),
##        y0 = c(atdn.19$tovo.S$CIs[4,]),
##        x1 = c(0,0),
##        y1 = c(atdn.19$tovo.S$CIs[4,]),
##        lty=2)
## par(las = 1)
## mtext("TNB", side = 4, line = 0.5, cex = 1.1)
## par(las = 0)
## ## LSE
## par(mar = c(3, 0, 2, 0))
## with(subset(bias.lse.13$ls$estimates, S.est.clump>atdn.13$S.ulrich$S[1,3] & S.est.clump<atdn.13$S.ulrich$S[1,4]),
##      boxplot(S, axes=FALSE,  ylim = c(1.5e4,1.95e4), border="red"))
## with(subset(bias.lse.13$ls$estimates, S.est.rnd>atdn.13$S.ulrich$S[1,3] & S.est.rnd<atdn.13$S.ulrich$S[1,4]),
##      boxplot(S, axes=FALSE, ylim = c(1.5e4,1.95e4), border="blue"))
## par(mar = c(3, 2, 2, 4))
## plot(S ~ S.est.rnd, data = bias.lse.13$ls$estimates,
##      ylim = c(1.5e4,1.95e4), xlim = c(1.4e4,1.6e4),
##      col="blue", ylab = "", xlab = "", type="n", main = "")
## rect(atdn.13$S.ulrich$S[1,3],0,atdn.13$S.ulrich$S[1,4],max(bias.lse.13$ls$estimates$S)*1.1,col="lightgrey", border=NA)
## box()
## axis(1)
## points(S ~ S.est.rnd, data = bias.lse.13$ls$estimates, col="blue")
## points(S ~ S.est.clump, data = bias.lse.13$ls$estimates, col="red")
## abline(0,1)
## segments(x0 = unlist(atdn.13$S.ulrich$S[1,3:4]),
##        y0 = unlist(atdn.13$S.ulrich$S[1,3:4]),
##        x1 = c(0,0),
##        y1 = unlist(atdn.13$S.ulrich$S[1,3:4]),
##        lty=2)
## ## 2013 updated taxonomy
## par(mar = c(3, 0, 2, 0))
## with(subset(bias.lse.13t$ls$estimates, S.est.clump>atdn.13.tax$S.ulrich$S[1,3] & S.est.clump<atdn.13.tax$S.ulrich$S[1,4]),
##      boxplot(S, axes=FALSE, ylim = c(1.25e4,1.65e4), border="red"))
## with(subset(bias.lse.13t$ls$estimates, S.est.rnd>atdn.13.tax$S.ulrich$S[1,3] & S.est.rnd<atdn.13.tax$S.ulrich$S[1,4]),
##      boxplot(S, axes=FALSE, ylim = c(1.25e4,1.65e4), border="blue"))
## par(mar = c(3, 2, 2, 4))
## plot(S ~ S.est.rnd, data = bias.lse.13t$ls$estimates,
##      ylim = c(1.25e4,1.65e4), xlim = c(1.3e4,1.45e4),
##      col="blue", ylab = "", xlab = "", type="n", main = "")
## rect(atdn.13.tax$S.ulrich$S[1,3],0,atdn.13.tax$S.ulrich$S[1,4],max(bias.lse.13t$ls$estimates$S),col="lightgrey", border=NA)
## box()
## axis(1)
## points(S ~ S.est.rnd, data = bias.lse.13t$ls$estimates, col="blue")
## points(S ~ S.est.clump, data = bias.lse.13t$ls$estimates, col="red")
## abline(0,1)
## segments(x0 = unlist(atdn.13.tax$S.ulrich$S[1,3:4]),
##        y0 = unlist(atdn.13.tax$S.ulrich$S[1,3:4]),
##        x1 = c(0,0),
##        y1 = unlist(atdn.13.tax$S.ulrich$S[1,3:4]),
##        lty=2)
## ## 2019
## par(mar = c(3, 0, 2, 0))
## with(subset(bias.lse.19$ls$estimates, S.est.clump>atdn.19$S.ulrich$S[1,3] & S.est.clump<atdn.19$S.ulrich$S[1,4]),
##      boxplot(S, axes=FALSE, ylim = c(1.4e4,1.72e4), border="red"))
## with(subset(bias.lse.19$ls$estimates, S.est.rnd>atdn.19$S.ulrich$S[1,3] & S.est.rnd<atdn.19$S.ulrich$S[1,4]),
##      boxplot(S, axes=FALSE, ylim = c(1.4e4,1.72e4),border="blue"))
## par(mar = c(3, 2, 2, 4))
## plot(S ~ S.est.rnd, data = bias.lse.19$ls$estimates,
##      ylim = c(1.4e4,1.72e4), xlim = c(1.4e4,1.55e4),
##      col="blue", ylab = "", xlab = "", type="n", main = "")
## rect(atdn.19$S.ulrich$S[1,3],0,atdn.19$S.ulrich$S[1,4],max(bias.lse.19$ls$estimates$S),col="lightgrey", border=NA)
## box()
## axis(1)
## points(S ~ S.est.rnd, data = bias.lse.19$ls$estimates, col="blue")
## points(S ~ S.est.clump, data = bias.lse.19$ls$estimates, col="red")
## abline(0,1)
## segments(x0 = unlist(atdn.19$S.ulrich$S[1,3:4]),
##        y0 = unlist(atdn.19$S.ulrich$S[1,3:4]),
##        x1 = c(0,0),
##        y1 = unlist(atdn.19$S.ulrich$S[1,3:4]),
##        lty=2)
## par(las = 1)
## mtext("LSE", side = 4, line = 0.5, cex = 1.1)
## par(las = 0)
## ## CHAO #
## ## 2013
## par(mar = c(3, 0, 2, 0))
## with(subset(bias.chao.13$ls$estimates, S.est.clump>atdn.13$Chao[4] & S.est.clump<atdn.13$Chao[5]),
##      boxplot(S, axes=FALSE, ylim = c(5e3,1.6e4), border="red"))
## with(subset(bias.chao.13$ls$estimates, S.est.rnd>atdn.13$Chao[4] & S.est.rnd<atdn.13$Chao[5]),
##      boxplot(S, axes=FALSE, ylim = c(5e3,1.6e4), border="blue"))
## par(mar = c(3, 2, 2, 4))
## plot(S ~ S.est.rnd, data = bias.chao.13$ls$estimates, ylim = c(5e3,1.6e4),
##      xlim = atdn.13$Chao[4:5]*c(.95,1.05), col="blue", ylab = "", xlab = "", type="n", main = "")
## rect(atdn.13$Chao[4],0,atdn.13$Chao[5],max(bias.chao.13$ls$estimates$S),col="lightgrey", border=NA)
## box()
## axis(1)
## points(S ~ S.est.rnd, data = bias.chao.13$ls$estimates, col="blue")
## points(S ~ S.est.clump, data = bias.chao.13$ls$estimates, col="red")
## abline(0,1)
## segments(x0 = c(atdn.13$Chao[4:5]),
##        y0 = c(atdn.13$Chao[4:5]),
##        x1 = c(0,0),
##        y1 = c(atdn.13$Chao[4:5]),
##        lty=2)
## ## 2013 updated taxonomy
## par(mar = c(3, 0, 2, 0))
## with(subset(bias.chao.13t$ls$estimates, S.est.clump>atdn.13.tax$Chao[4] & S.est.clump<atdn.13.tax$Chao[5]),
##      boxplot(S, axes=FALSE, ylim = c(5e3,1.6e4), border="red"))
## with(subset(bias.chao.13t$ls$estimates, S.est.rnd>atdn.13.tax$Chao[4] & S.est.rnd<atdn.13.tax$Chao[5]),
##      boxplot(S, axes=FALSE, ylim = c(5e3,1.6e4), border="blue"))
## par(mar = c(3, 2, 2, 4))
## plot(S ~ S.est.rnd, data = bias.chao.13t$ls$estimates, ylim = c(5e3,1.6e4),
##      xlim = atdn.13.tax$Chao[4:5]*c(.95,1.05), col="blue", ylab = "", xlab = "", type="n", main = "")
## rect(atdn.13.tax$Chao[4],0,atdn.13.tax$Chao[5],max(bias.chao.13t$ls$estimates$S),col="lightgrey", border=NA)
## box()
## axis(1)
## points(S ~ S.est.rnd, data = bias.chao.13t$ls$estimates, col="blue")
## points(S ~ S.est.clump, data = bias.chao.13t$ls$estimates, col="red")
## abline(0,1)
## segments(x0 = c(atdn.13.tax$Chao[4:5]),
##        y0 = c(atdn.13.tax$Chao[4:5]),
##        x1 = c(0,0),
##        y1 = c(atdn.13.tax$Chao[4:5]),
##        lty=2)
## ## 2019
## par(mar = c(3, 0, 2, 0))
## with(subset(bias.chao.19$ls$estimates, S.est.clump>atdn.19$Chao[4] & S.est.clump<atdn.19$Chao[5]),
##      boxplot(S, axes=FALSE, ylim = c(5e3,1.6e4), border="red"))
## with(subset(bias.chao.19$ls$estimates, S.est.rnd>atdn.19$Chao[4] & S.est.rnd<atdn.19$Chao[5]),
##      boxplot(S, axes=FALSE, ylim = c(5e3,1.6e4), border="blue"))
## par(mar = c(3, 2, 2, 4))
## plot(S ~ S.est.rnd, data = bias.chao.19$ls$estimates, ylim = c(5e3,1.6e4),
##      xlim = atdn.19$Chao[4:5]*c(.95,1.05), col="blue", ylab = "", xlab = "", type="n", main = "")
## rect(atdn.19$Chao[4],0,atdn.19$Chao[5],max(bias.chao.19$ls$estimates$S),col="lightgrey", border=NA)
## box()
## axis(1)
## points(S ~ S.est.rnd, data = bias.chao.19$ls$estimates, col="blue")
## points(S ~ S.est.clump, data = bias.chao.19$ls$estimates, col="red")
## abline(0,1)
## segments(x0 = c(atdn.19$Chao[4:5]),
##        y0 = c(atdn.19$Chao[4:5]),
##        x1 = c(0,0),
##        y1 = c(atdn.19$Chao[4:5]),
##        lty=2)
## par(las = 1)
## mtext("CHAO", side = 4, line = 0.1, cex = 1)
## par(las = 0)
## ## X axis
## par(mar=c(2,0,1,0))
## plot(0,0, xlim=c(0,1), ylim =c(0,1), axes=FALSE, xlab="", ylab="", type="n")
## mtext(text = "Estimated number of species from simulated samples", side=1, cex = 1.5)
## dev.off()
################################################################################


################################################################################
## Dotplots of estimated Species Richness  and ABC estimates##
################################################################################

pdf("figs_and_tables/bias_corrected_estimates_dotchart.pdf")
par(cex.main = 1.5,
    cex.lab = 1.25,
    font.lab = 2,
    cex.axis = 1.25)
## Bias-corrected estimates (assuming clumped sampling) ##
tmp1 <- S.estimates.bc[S.estimates.bc$sampling=="clump"&S.estimates.bc$type!="LSE TNB",]
tmp1$type[tmp1$type=="LSE LS"] <- "LSE"
tmp1$type <- factor(tmp1$type, levels=c("TNB", "CHAO","LS","ABC","LSE"))
tmp1 <- tmp1[order(tmp1$type,tmp1$dataset),]
with(tmp1,
     dotchart(mean, labels = dataset, groups = type, color=2:4, pch=1,
              pt.cex = 1.5,
              xlab = "Estimated species richness",
              xlim=range(c(IC.low, IC.up),na.rm=TRUE))
     )
Ys <- 21:23 - rep(seq(0,20, by=5), each=3)
## ## Adds uncorrected values
## points(S.estimates$mean, Ys[-(10:12)], col=rep(2:4,3), pch="|")
## Adds error bars
cores <- rep(2:4,5)
for(i in 1:nrow(tmp1))
    segments(x0=tmp1$IC.low[i], x1=tmp1$IC.up[i], y0=Ys[i], y1=Ys[i],
             col=cores[i], lwd =1.5)
dev.off()


## Original (biased) estimates ##
pdf("figs_and_tables/uncorrected_estimates_dotchart.pdf")
par(cex.main = 1.5,
    cex.lab = 1.25,
    font.lab = 2,
    cex.axis = 1.25)
with(S.estimates,
     dotchart(mean, labels = dataset, groups = type, color=2:4, pch=1,
              pt.cex = 1.5,
              xlim=range(c(S.estimates$IC.low, S.estimates$IC.up),na.rm=TRUE),
              xlab = "Estimated species richness")
     )
Ys <- 16:18 - rep(seq(0,15, by=5), each=3)
cores <- rep(2:4,4)
for(i in 1:nrow(S.estimates))
    segments(x0=S.estimates$IC.low[i], x1=S.estimates$IC.up[i], y0=Ys[i], y1=Ys[i],
             col=cores[i], lwd =1.5)
dev.off()

## Bias corrected with both types of sampling (random and clumped) ##
## Auxiliary function
f11 <- function(x, type, ...){
    tmp3 <- x[x$type==type,]
    with(tmp3,
         dotchart(mean, labels = sampling, groups = factor(dataset),
                  color=2:4,
                  pch=1,
                  pt.cex = 1.5,
                  xlim=range(c(IC.low, IC.up),na.rm=TRUE),
                  ...
                  )
     )
    ## Errors bars
    Ys <- 11:13 - rep(seq(0,10, by=5), each=3)
    cores <- rep(2:4,4)
    for(i in 1:nrow(tmp3))
        segments(x0=tmp3$IC.low[i], x1=tmp3$IC.up[i], y0=Ys[i], y1=Ys[i],
                 col=cores[i], lwd =1.5)
    }
#####################################################################################
tmp2 <- S.estimates.all[S.estimates.all$type!="LSE TNB"&S.estimates.all$type!="ABC",]
tmp2$type[tmp2$type=="LSE LS"] <- "LSE"
tmp2$sampling[is.na(tmp2$sampling)] <- "Not corrected"
tmp2 <- tmp2[order(tmp2$dataset,tmp2$sampling),]

pdf("figs_and_tables/bias_corrected_estimates_rnd_clump_dotchart.pdf", width = 9, height = 9)
par(cex.main = 1.5,
    cex.lab = 1.25,
    font.lab = 2,
    cex.axis = 1.5,
    mfrow = c(2,2))
f11(tmp2, "LS", xlab = "Estimated species richness", main= "Logseries")
f11(tmp2, "TNB", xlab = "Estimated species richness", main= "Truncated Negative Binomial")
f11(tmp2, "LSE", xlab = "Estimated species richness", main= "RAD linear extension")
f11(tmp2, "CHAO", xlab = "Estimated species richness", main= "CHAO")
dev.off()


################################################################################
## Plots of Population rad expected by ABC ##
################################################################################

## 2013
S1 <- S.estimates.all[S.estimates.all$type=="ABC"&S.estimates.all$dataset=="2013"&S.estimates.all$sampling=="clump",4:6]
abc13.sim.pop <- with(atdn.13,
                 lapply(S1, sim.abc, sad = "ls",
                        N = Tot.t,  tot.area = Tot.A, n.plots = N.plots,
                        lm.sd.fit = lm.sd, lmk.fit = lm.k, nb.fit =y.nb2,
                        summary = FALSE, upper=1e16) )
## 2013 with updated taxonomy
S1 <- S.estimates.all[S.estimates.all$type=="ABC"&S.estimates.all$dataset=="2013 updated"&S.estimates.all$sampling=="clump",4:6]
abc13t.sim.pop <- with(atdn.13.tax,
                      lapply(S1,
                             sim.abc, sad = "ls",
                             N = Tot.t,  tot.area = Tot.A, n.plots = N.plots,
                             lm.sd.fit = lm.sd, lmk.fit = lm.k, nb.fit =y.nb2,
                             summary = FALSE, upper=1e16) )
## 2019
S1 <- S.estimates.all[S.estimates.all$type=="ABC"&S.estimates.all$dataset=="2019"&S.estimates.all$sampling=="clump",4:6]
abc19.sim.pop <- with(atdn.19,
                      lapply(S1,
                             sim.abc, sad = "ls",
                             N = Tot.t,  tot.area = Tot.A, n.plots = N.plots,
                             lm.sd.fit = lm.sd, lmk.fit = lm.k, nb.fit =y.nb2,
                             summary = FALSE, upper=1e16) )


## The plots for  ##
pdf("figs_and_tables/population_RADS_and_ABC_predictions%d.pdf", onefile=TRUE)
par(cex.main = 1.5, lwd=2,
    mar = c(5, 4, 4, 2) + 0.1,
    mgp = c(3.5, 1, 0), 
    cex.lab = 1.25, font.lab = 2, cex.axis = 1.2, bty = "l", las = 1,
    oma=c(3,3,0,0))
## 2013
plot(rad(atdn.13$data$population), xlim=c(1,6e3), col="grey", ylab="")
par(las=0)
mtext("Species abundances", 2, cex =1.25, font = 2, line = 4.5)
par(las=1)
lines(rad(abc13.sim.pop$mean$clump.samp[,2]))
lines(rad(abc13.sim.pop$IC.low$clump.samp[,2]), lty=2)
lines(rad(abc13.sim.pop$IC.up$clump.samp[,2]), lty=2)
## 2013 updated
plot(rad(atdn.13.tax$data$population), xlim=c(1,6e3), col="grey", ylab = "")
par(las=0)
mtext("Species abundances", 2, cex =1.25, font = 2, line = 4.5)
par(las=1)
lines(rad(abc13t.sim.pop$mean$clump.samp[,2]))
lines(rad(abc13t.sim.pop$IC.low$clump.samp[,2]), lty=2)
lines(rad(abc13t.sim.pop$IC.up$clump.samp[,2]), lty=2)
## 2019
plot(rad(atdn.19$data$population), xlim=c(1,6e3), col="grey", ylab = "")
par(las=0)
mtext("Species abundances", 2, cex =1.25, font = 2, line = 4.5)
par(las=1)
lines(rad(abc19.sim.pop$mean$clump.samp[,2]))
lines(rad(abc19.sim.pop$IC.low$clump.samp[,2]), lty=2)
lines(rad(abc19.sim.pop$IC.up$clump.samp[,2]), lty=2)
dev.off()

################################################################################
## Plots of aggregation coefficients x mean population abundances
################################################################################
##Utility function
f12 <- function(obj, ...){
    with(obj,{
        cf1 <- coef(lm.k)
        plot(k ~ dens.ha, data= data, log = "xy", ...)
        curve(exp(cf1[1])*x^cf1[2], add=TRUE, col = 4)
    }
    )
}
###################################
yl1 <- range(atdn.13$data$k,atdn.13.tax$data$k,atdn.19$data$k[atdn.19$data$k<1])
xl1 <- range(atdn.13$data$dens.ha,atdn.13.tax$data$dens.ha,atdn.19$data$dens.ha)
pdf("figs_and_tables/kXdens.pdf", width = 12, height = 4.25)
par(mar = c(5, 5, 4, 2) + 0.1,
    mgp = c(3.5, 1, 0),
    oma=c(3,3,0,0),
    las = 0,
    bty = "l", 
    cex.main = 1.5,  
    cex.lab = 1.4, font.lab = 2, cex.axis = 1.25,
    lwd = 2,
    mfrow=c(1,3))
f12(atdn.13,
    ylab = "",
    xlab="",
    xlim = xl1, ylim = yl1,
    main = "2013", col="grey", cex=0.25
    )
mtext("Dispersion coefficient", cex=1.5, font =2, side = 2, line=5)
f12(atdn.13.tax, ylab = "",
    xlab="",
    xlim = xl1, ylim = yl1,
    main = "2013 updated", col="grey", cex=0.25)
mtext(expression(paste("Density (",ha^{-1},")")), cex=1.5, font = 2, side = 1, line= 6)
f12(atdn.19,
    ylab = "",
    xlab="",
    xlim = xl1, ylim = yl1,
    main = "2019", col="grey", cex=0.25)
dev.off()


################################################################################
## Plots of sd deviations  x mean of population abundances
################################################################################
f13 <- function(obj, ...){
    with(obj,{
        cf1 <- coef(lm.sd)
        plot(pop.sd ~ population, data= data, log = "xy", ...)
        curve(exp(cf1[1])*x^cf1[2], add=TRUE, col = 4)
    }
    )
}

yl1 <- range(atdn.13$data$pop.sd,atdn.13.tax$data$pop.sd,atdn.19$data$pop.sd[atdn.19$data$pop.sd<1])
xl1 <- range(atdn.13$data$population,atdn.13.tax$data$population,atdn.19$data$population)
pdf("figs_and_tables/sdXpopsize.pdf", width = 12, height = 4.25)
par(mar = c(5, 5, 4, 2) + 0.1,
    mgp = c(3.5, 1, 0),
    oma=c(3,3,0,0),
    las = 0,
    bty = "l", 
    cex.main = 1.5,  
    cex.lab = 1.4, font.lab = 2, cex.axis = 1.25,
    lwd = 2,
    mfrow=c(1,3))
f13(atdn.13,
    ylab = "",
    xlab="",
    xlim = xl1, ylim = yl1,
    main = "2013", col="grey", cex=0.25
    )
mtext("Standard deviation of estimated pop size", cex=1.5, font =2, side = 2, line=5)
f13(atdn.13.tax, ylab = "",
    xlab="",
    xlim = xl1, ylim = yl1,
    main = "2013 updated", col="grey", cex=0.25)
mtext("Estimated population size", cex=1.5, font = 2, side = 1, line= 6)
f13(atdn.19,
    ylab = "",
    xlab="",
    xlim = xl1, ylim = yl1,
    main = "2019", col="grey", cex=0.25)
dev.off()


################################################################################
## Plots of model selection biases
################################################################################

pdf("figs_and_tables/bias_model_selection.pdf", width = 9, height = 5)
par(mar = c(5, 5, 4, 2) + 0.1,
    mgp = c(3.5, 1, 0),
    oma=c(3,3,0,0),
    las = 1,
    bty = "l", 
    cex.main = 1.15,  
    cex.lab = 1, font.lab = 2, cex.axis = 1,
    lwd = 2,
    mfrow=c(1,2))
plot(ls.rnd.w ~ S, data = bias.msel.19$ls, ylim = c(0,1),
     xlab = "Species richness in the community",
     ylab = "Evidence weight of logseries fit",
     main = "Logseries community")
points(ls.clump.w ~ S, data = bias.msel.19$ls, col = 2)
plot(ls.rnd.w ~ S, data = bias.msel.19$tnb, ylim = c(0,1),
     xlab = "Species richness in the community",
     ylab = "Evidence weight of logseries fit",
     main = "Negative binomial community")
points(ls.clump.w ~ S, data = bias.msel.19$tnb, col = 2)
legend("topright", c("Random", "Clumped"), col=c(1,2), cex = 1, pch =1, bty="n")
dev.off()


################################################################################
## Population RAD and regional LS
################################################################################
## Species richness estimates and CI's to use:
## from ABC
## S1 <- unlist(S.estimates.all[S.estimates.all$type=="ABC"&S.estimates.all$dataset=="2019"&S.estimates.all$sampling=="clump",4:6])
## Mean of estimates weighted by uncertainty, excluding TNB. ICs as the range of IC
S1  <- 
    read.csv("figs_and_tables/estimates_S_table.csv") %>%
    filter(bias.corrected==TRUE&type!="LSE TNB"&sampling=="clump"&type!="TNB"&dataset=="2019") %>%
    mutate(se = (IC.up-IC.low)/4)%>%
    summarise(w.mean = weighted.mean(mean, w=1/se), min=min(IC.low), max=max(IC.up)) %>%    
    unlist()

## Simulated regional rads
ls.m <- rad.ls(S = S1[1], N = atdn.19$Tot.t)
ls.low <- rad.ls(S = S1[2], N = atdn.19$Tot.t)
ls.up <- rad.ls(S = S1[3], N = atdn.19$Tot.t)
## Simulated population sizes
abc19.sim.pop <- with(atdn.19,
                            lapply(S1, sim.abc, N = Tot.t, sad = "ls", tot.area = Tot.A, n.plots = N.plots,
                                   lm.sd.fit = lm.sd, lmk.fit = lm.k, nb.fit = y.nb2, summary = FALSE))
## Hyperdominants
hd <- round(atdn.19$HD$summary[c(2,4:5)])
## The figure
pdf("figs_and_tables/pop_rad_with_predicted_rad.pdf")
par(bty = "l", lwd = 2)
plot(rad(atdn.19$data$population), xlim = c(1, S1[3]),
     ##ylim=c(min(c(ls.m$y, ls.up$y, ls.low$y)), max(atdn.19$data$population)),
     ylim = c(1, max(atdn.19$data$population)),
     axes=FALSE , #yaxs="i",
     col=rep(c("red","darkgrey"), c(hd[1], length(atdn.19$data$population)-hd[1])),
     cex.lab = 1.25, font.lab = 2, cex.axis = 1.2)
axis(1, at=c(1, seq(2500,17500, by=2500)))
axis(2)
lines(rad(ceiling(ls.m$y)), col="darkblue")
lines(rad(ls.low$y), lty=2, col="darkblue")
lines(rad(ls.up$y), lty=2, col="darkblue")
text(1.7e4, 1e9, paste("Nr of hyperdominant species = ",hd[1]," (",hd[2]," - ",hd[3],")", sep=""), font=2, pos=2)
text(1.7e4, 0.25e9, paste("Estimated tree species = ",round(S1[1])," (",round(S1[2])," - ",round(S1[3]),")", sep=""), font=2, pos=2)
par(fig=c(0.09,0.49,0.05,0.50),
    mgp = c(2.5,0.5,0),
    ##mar = c(5,4,4,2),
    new = T, cex.axis = 0.8, cex.lab=1,
    ##yaxp=c(1,3,3),
    bty="o")
plot(rad(atdn.19$data$population), xlim=c(1,6e3),
     col="grey", ylab = "", xlab="", cex=0.75, axes = F)
axis(2, at=c(1e6, 1e9), labels=c("","1e+09"))
axis(1, at=c(1,2500,5000), labels=c("","","5000"))
lines(rad(abc19.sim.pop[[1]]$clump.samp[,2]))
lines(rad(abc19.sim.pop[[2]]$clump.samp[,2]), lty=2)
lines(rad(abc19.sim.pop[[3]]$clump.samp[,2]), lty=2)
par(fig=c(0,1,0,1))
dev.off()

################################################################################
## Simulating number of species recorded with sample increasing
################################################################################
## assembles the data for clumped sampling in a dataframe
sim.rich <- with(S.proj.19,
                 data.frame(
                     n.plots = n.plots, 
                     mean = matrix(unlist(S.est.mean), byrow=TRUE, ncol=2)[,2],
                     low = matrix(unlist(S.est.low), byrow=TRUE, ncol=2)[,2],
                     upp= matrix(unlist(S.est.upp), byrow=TRUE, ncol=2)[,2])
                 )
## Estimated species richness from Fisher's logseries
a1 <- fishers.alpha(N=atdn.19$Tot.t, S = S1[1])
## Density per hectare
d1 <- with(atdn.19, Tot.t/Tot.A)
a1*log(1 + (max(S.proj.19$n.plots)*d1)/a1)
## The plot
pdf("figs_and_tables/richness_extrapolations.pdf")
par(mar = c(5, 5, 4, 2) + 0.1,
    mgp = c(3.5, 1, 0),
    oma=c(3,3,0,0),
    las = 0,
    bty = "l", 
    cex.main = 1.5,  
    cex.lab = 1.4, font.lab = 2, cex.axis = 1.25,
    lwd = 2)
plot(mean ~ n.plots, data = sim.rich,
     ylim = range(sim.rich[, -1]), type = "l",
     xlab = "Sample size (1-ha plots)", ylab = "Expected Nr. of species in the sample")
lines(low ~ n.plots, data=sim.rich, lty = 2)
lines(upp ~ n.plots, data=sim.rich, lty = 2)
points(c(atdn.13$N.plots, atdn.13.tax$N.plots, atdn.19$N.plots),
       c(atdn.13$Sobs, atdn.13.tax$Sobs, atdn.19$Sobs), cex=1.5, pch=19, col="blue")
text(c(atdn.13$N.plots, atdn.13.tax$N.plots, atdn.19$N.plots),
     c(atdn.13$Sobs, atdn.13.tax$Sobs, atdn.19$Sobs), c("2013", "2013 updated", "2019"),
     pos = 4, col="blue")
dev.off()       

################################################################################
## Fits of sads to abundances in the sample - Diagnostic plots
################################################################################
## Auxiliary functions
f19 <- function(obj, legend=FALSE, ...){
    obs.oc <- octav(obj$data$N.ind)
    ls.op <- octavpred(obj$y.ls)
    tnb.op <- octavpred(obj$y.nb2)
    pln.op <- octavpred(obj$pln)
    plot(obs.oc, ylim = range(c(ls.op$Freq, tnb.op$Freq, pln.op$Freq, obs.oc$Freq)), ...)
    lines(ls.op, col= 2)
    lines(tnb.op, col=3)
    lines(pln.op, col =4)
    if(legend)
        legend("topright", c("LS","TNB", "PLN"), pch=1, lty=1, col=c(2:4), bty="n")
}

f20 <- function(obj, legend=FALSE, ...){
    obs.r <- rad(obj$data$N.ind)
    ls.rp <- radpred(obj$y.ls)
    tnb.rp <- radpred(obj$y.nb2)
    pln.rp <- radpred(obj$pln)
    plot(obs.r, ylim = range(c(ls.rp$abund, tnb.rp$abund, pln.rp$abund, obs.r$abund), na.rm=TRUE), ...)
    lines(ls.rp, col= 2)
    lines(tnb.rp, col=3)
    lines(pln.rp, col =4)
    if(legend)
        legend("topright", c("LS","TNB", "PLN"), lty=1, col=c(2:4), bty="n")
}

## Octav and rad plots
##png("figs_and_tables/sads_fit_to_samples.png", width=9, height=7, res=150, units="in")
pdf("figs_and_tables/sads_fit_to_samples.pdf", width=12, height=7)
par(##mar = c(5, 5, 4, 2) + 0.1,
    ##mgp = c(3.25, 1, 0),
    ##oma=c(3,3,0,0),
    ##las = 1,
    bty = "l", 
    cex.main = 1.15,  
    cex.lab = 1, font.lab = 2, cex.axis = 1,
    lwd = 2,
    mfrow=c(2,3))
f19(atdn.13, main = "2013", xlab="", cex.lab=1.5)
f19(atdn.13.tax, main = "2013 updated", ylab="", cex.lab=1.5)
f19(atdn.19, legend=TRUE, main = "2019", xlab="", ylab = "")
f20(atdn.13, log="xy",  xlab="", cex.lab=1.5, col="grey")
f20(atdn.13.tax, log="xy", ylab="", cex.lab=1.5, col="grey")
f20(atdn.19, log="xy", xlab="", ylab = "", col="grey")
dev.off()

## qq plots
pdf("figs_and_tables/samples_qqplots.pdf", width=9, height=9)
par(mar = c(6, 6, 4, 2) + 0.1,
    ##mgp = c(3.25, 1, 0),
    ##oma=c(3,3,0,0),
    ##las = 1,
    bty = "l", 
    cex.main = 1.5,  
    cex.lab = 1, font.lab = 2, cex.axis = 1,
    lwd = 2,
    mfrow=c(3,3))
qqsad(atdn.13$y.ls, main="2013", xlab="", ylab="", col="grey")
qqsad(atdn.13.tax$y.ls, main="2013 updated", xlab="", ylab="", col="grey")
qqsad(atdn.19$y.ls, main="2019", xlab="", ylab="", col="grey")
mtext("LS", at=8e3)
qqsad(atdn.13$y.nb2, main="", xlab="", ylab="", cex.lab=1.7, col="grey")
mtext("Sample quantiles", cex=2, side= 2, line = 3)
qqsad(atdn.13.tax$y.nb2, main="", xlab="", ylab="", col="grey")
qqsad(atdn.19$y.nb2, main="", xlab="", ylab="", col="grey")
mtext("TNB", at=7e3)
qqsad(atdn.13$pln, main="", xlab="", ylab="", col="grey")
qqsad(atdn.13.tax$pln, main="", xlab="", cex.lab=1.7, ylab="", col="grey")
mtext("Theoretical quantiles", cex=2, side =1, line = 3.5)
qqsad(atdn.19$pln, main="", xlab="", ylab="", col="grey")
mtext("PLN", at=1.8e5)
dev.off()

## ## Sample RAD and predcited by LS, TNB, PLN with R inset of the same plot in log-log scale
## pdf("figs_and_tables/samp_sad_with_inset_log_log.pdf")
## par(mar = c(5, 5, 4, 2) + 0.1,
##     mgp = c(3.5, 1, 0),
##     oma=c(3,3,0,0),
##     las = 1,
##     bty = "l", 
##     cex.main = 1.15,  
##     cex.lab = 1, font.lab = 2, cex.axis = 1,
##     lwd = 2,
##     mfrow=c(1,1))
## plot(obs.r, col="grey", ylim=c(1, max(pln.r$abund)))
## lines(ls.r)
## lines(tnb.r, col="red")
## lines(pln.r, col="green")
## par(fig=c(0.39,0.99,0.39,0.99),
##     mgp = c(2.5,0.5,0),
##     ##mar = c(5,4,4,2),
##     new = T, cex.axis = 0.8, cex.lab=1,
##     ##yaxp=c(1,3,3),
##     bty="o")
## plot(obs.r, col="grey", log = "xy", xlab="", ylab="",
##      ylim=c(1, max(pln.r$abund)), axes=FALSE)
## lines(ls.r)
## lines(tnb.r, col="red")
## lines(pln.r, col="green")
## par(fig=c(0,1,0,1))
## dev.off()
