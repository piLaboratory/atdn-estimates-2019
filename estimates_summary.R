## ----setup, echo=FALSE, include=FALSE------------------------------------
library(VGAM)
library(untb)
library(fitdistrplus)
library(sads)
library(knitr)
library(SPECIES)
library(xtable)
source("functions.R")
opts_chunk$set(fig.align = 'center', fig.show = 'hold',
               fig.height = 6.5, warning = FALSE, message = FALSE,
               error = FALSE, echo=TRUE)
options(formatR.arrow = TRUE, width = 90, cache=TRUE)

## ----Data prep-----------------------------------------------------------
dados <- read.csv2("data.csv", as.is=TRUE)
y <- dados$N.ind
Sobs <- length(y)
## Total number of trees (average density x area)
Tot.t <- 567*5.5e8
## Proportion of total trees in the sample
p1 <- sum(dados$N.ind)/Tot.t
## Total number of plots
N.plots <- 1945
## Total area hectares
Tot.A <- 5.79e8
## Sampled area ha
Samp.A <- 2.048e3
## Amazon RAD sent by Hans
load("steege_files/AmazonRAD.RData")
## Truncated Negative binomial 
## (already generated, too slow to rerun)
nb.pred.full <- read.csv("NB_RAD.csv")$x

## ----fit pln, cache=TRUE-------------------------------------------------
pln <- fitpoilog(y)
par(mfrow=c(2,2))
plot(pln)
par(mfrow=c(1,1))

## ----PLN species richness estimate 2-------------------------------------
(pln.d0 <- dpoilog(0, mu = pln.cf[1], sig=pln.cf[2]))

## ----fit ls--------------------------------------------------------------
y.ls <- fitls(y)
par(mfrow=c(2,2))
plot(y.ls)
par(mfrow=c(1,1))

## ----estimated S ls------------------------------------------------------
alpha <- coef(y.ls)[[2]]
(S.ls <- alpha*log(1 + Tot.t/alpha))

## ----ls and S est for ls-------------------------------------------------
(ls.ci <- confint(y.ls))
## Estimated species richness for lower bound of alpha's IC
ls.ci[1]*log(1 + Tot.t/ls.ci[1])
ls.ci[2]*log(1 + Tot.t/ls.ci[2])

## ----fit NB--------------------------------------------------------------
## With VGAM
y.nb <- vglm(y ~ 1, posnegbinomial)
## With sads
y.nb2 <- fitnbinom(y, start.value=c(size=0.3, mu=mean(y)))
## Comparing: 
exp(coef(y.nb))
coef(y.nb2)

## ----nb plots------------------------------------------------------------
par(mfrow=c(2,2))
plot(y.nb2)
par(mfrow=c(1,1))

## ----Tovo S estimate-----------------------------------------------------
cf.nb <- coef(y.nb2)
csi.p <- unname(cf.nb[2]/(sum(cf.nb)))
csi <- csi.p/(p1+(1-p1)*csi.p)
## Estimated number of species 
S.nb <- Sobs*(1-(1-csi)^cf.nb[1]) / (1-(1-csi.p)^cf.nb[1])
(S.nb <- unname(S.nb))

## ----NB S est CI---------------------------------------------------------
(tovo.S <- tovo(fit = y.nb2, p = p1, CI=TRUE))

## ----model selection-----------------------------------------------------
AICtab(pln, y.nb2, y.ls, base=TRUE)

## ----Comparing LS and NB, echo=FALSE-------------------------------------
par(mfrow=c(2,2))
ls.rad <- radpred(y.ls)
nb.rad <- radpred(y.nb2)
plot(ls.rad$abund, nb.rad$abund, log="xy",
     xlab="Logseries theor. quantiles",
     ylab="TNB theor. quantiles", 
     col="grey", cex=0.25)
abline(0,1, col="blue")
plot(rad(dados$N.ind), col="grey", cex=0.5)
lines(ls.rad)
lines(nb.rad, col="red")
legend("topright", c("LS", "NB"), col=c("blue","red"), 
       lty=1, bty="n")
ls.oc <- octavpred(y.ls)
nb.oc <- octavpred(y.nb2)
plot(octav(dados$N.ind), ylim=range(c(ls.oc$Freq, nb.oc$Freq)))
lines(ls.oc)
lines(nb.oc, col="red")
legend("topright", c("LS", "NB"), col=c("blue","red"), 
       lty=1, bty="n")
par(mfrow=c(1,1))

## ----Linear extrapolation from regional RAD------------------------------
## Regional RAD
pop.rad <- rad(dados$population)
## Linear regression through central 50% quantiles of the RAD
p.lm <- lm(log(abund)~rank, data=data.frame(pop.rad), 
           subset=rank>max(rank)*.25&rank<max(rank)*.75)
## Regression coefficients
cf.p.lm <- unname(coef(p.lm))
## Logseries projection (upper bound)
(S.reg1 <- abs(cf.p.lm[1]/cf.p.lm[2]))

## ----regional rad lognormal----------------------------------------------
d <- log(max(dados$population))-cf.p.lm[1]
(S.reg2 <- abs((cf.p.lm[1]-d)/cf.p.lm[2]))

## ----nbinom rad, echo=FALSE----------------------------------------------
plot(rad(Amazon.rad), col="red", lwd=2, type="n",
     ylim=c(1, max(dados$population)))
points(rad(dados$population), col="grey")
curve(exp(cf.p.lm[1]+cf.p.lm[2]*x), add=TRUE)
curve(exp(cf.p.lm[1]-d+cf.p.lm[2]*x), add=TRUE)
lines(rad(nb.pred.full), col="blue", lwd=2)
lines(rad(Amazon.rad), col="red", lwd=2)
legend("topright", c("LS", "TNB", "Linear upper/lower bounds"), 
       col=c("red", "blue", "black"), lty=1, lwd=2, bty="n")

## ----k x dens regression-------------------------------------------------
## estimating k parameter of a NB for each species 
dados$dens.ha <- dados$N.ind/Samp.A
dados$k <- est.kv(mu=dados$dens.ha, 
                  nzeroes=N.plots-dados$N.plots, 
                  Nplots=N.plots)
lm.k <-lm(log(k)~log(dens.ha), 
          data=dados, subset=k<1)
## Estimated regression standard error
lm.k.sigma <- summary(lm.k)$sigma
summary(lm.k)

## ----lok k x log density plot, echo=FALSE--------------------------------
plot(log(k)~log(dens.ha), data=dados, 
     subset=k<1, xlab="Density (ind/ha)",
     ylab="Aggregation parameter of NB",
     col="grey")
abline(lm.k, col="blue", lwd=2)

## ----sim popsizes LS and NB, cache=TRUE----------------------------------
## Simulation of population sizes from samples of 
## the regional RADs
nrep <- 100 #number of repetitions
## Samples from LS RAD
## Matrices to store samples
tmp1 <- tmp2 <- matrix(0,nrow=length(Amazon.rad), ncol=nrep)
for(j in 1:nrep){
    ## Samples aggregation parameters for each species according to 
    ## the linear model + standard errors of this model
    k1 <- exp(rnorm(length(Amazon.rad), 
                    mean=Amazon.rad.lk, sd=lm.k.sigma))
    ## Random (Poisson) sample
    y1 <- mapply(sim.occ, mu = Amazon.rad/Tot.A,
                   MoreArgs=list(N = N.plots))
    ## Clumped (negative binomial) sample
    y2 <- mapply(sim.occ, mu = Amazon.rad/Tot.A, size = k1,
                   MoreArgs=list(N = N.plots, pois.samp=FALSE))
    ## Pick the population sizes for species
    ## recorded in the sample
    tmp1[1:sum(y1),j] <- Amazon.rad[y1>0]
    tmp2[1:sum(y2),j] <- Amazon.rad[y2>0]
}
## Simulation of samples from TNB 
tmp4 <- tmp3 <- matrix(0,nrow=length(nb.pred.full), ncol=nrep)
for(j in 1:nrep){
    k2 <- exp(rnorm(length(nb.pred.full), 
                    mean=nb.pred.full.lk, sd=lm.k.sigma))
    y3 <- mapply(sim.occ, mu = nb.pred.full/Tot.A,
                   MoreArgs=list(N = N.plots))
    y4 <- mapply(sim.occ, mu = nb.pred.full/Tot.A, size = k2,
                   MoreArgs=list(N = N.plots, pois.samp=FALSE))
    tmp3[1:sum(y3),j] <- nb.pred.full[y3>0]
    tmp4[1:sum(y4),j] <- nb.pred.full[y4>0]
}
## Mean RADs for each simulation
ls.rnd <- apply(tmp1,1,mean)
ls.clump <- apply(tmp2,1,mean)
nb.rnd <- apply(tmp3,1,mean)
nb.clump <- apply(tmp4,1,mean)

## ----plot sampled RADS, echo=FALSE---------------------------------------
plot(rad(dados$population), col="grey", log="y", xlim=c(1,5500),
     ylab = "Population size")
lines(rad(ls.rnd), lwd=2)
lines(rad(ls.clump), lwd=2, col="red")
lines(rad(nb.rnd), lwd=2, col="green")
lines(rad(nb.clump), lwd=2, col="black")
legend("topright", 
       c("LS, random sample", "LS, clumped sample", "NB, random sample", "NB, clumped sample"), 
       col=c("blue","red","green", "black"), lty=1, lwd=2, bty="n")

## ----mean squared errors, echo=FALSE, results="asis"---------------------
## Mean squared deviation of log abundances
MS.tab <- xtable(
    data.frame(
        MS= c(
            mean((log(dados$population) - log(ls.rnd[1:nrow(dados)]))^2),
            mean((log(dados$population) - log(nb.rnd[1:nrow(dados)]))^2),
            mean((log(dados$population) - log(ls.clump[1:nrow(dados)]))^2)
        ),
        row.names=c("LS random", "NB random", "LS clumped")
    ), digits=4)
## Inf because of many species with zero abundance in the sample
## mean((log(dados$population) - log(nb.clump[1:nrow(dados)]))^2)
print(MS.tab)    

