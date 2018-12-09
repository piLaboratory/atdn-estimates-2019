## ----setup, echo=FALSE, include=FALSE------------------------------------
library(VGAM)
library(untb)
library(fitdistrplus)
library(sads)
library(knitr)
library(SPECIES)
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
## Amazon RAD sent by Hans
load("steege_files/AmazonRAD.RData")

## ----fit pln-------------------------------------------------------------
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

## ----nbinom rad, cache=TRUE----------------------------------------------
## ppoints for estimated species richness
S.nb.p1 <- rev(ppoints(tovo.S$S.est))
## Predicted Negative Binomial for 100 points over the rad using 
## richness predicted from negbib lower limit mu
## Estimated richness value
nb.pred <- rad.posnegbin(S = tovo.S$S.est, size = cf.nb[1], 
                         prob = 1-tovo.csi(cf.nb[1], cf.nb[2], p1))
## Lower and upper limits
nb.pred.l <- rad.posnegbin(S = tovo.S$CIs[4,1], size = tovo.S$CIs[1,1], 
                           prob = tovo.S$CIs[3,1])
nb.pred.u <- rad.posnegbin(S = tovo.S$CIs[4,2], size = tovo.S$CIs[1,2], 
                           prob = tovo.S$CIs[3,2])
## The plots
plot(abund~rank, data=data.frame(pop.rad), xlim=c(1,S.reg1*1.1), 
     ylim=c(1,max(dados$population)), log="y", col="grey")
lines(nb.pred, col="blue")
lines(nb.pred.l, col="blue", lty=2)
lines(nb.pred.u, col="blue", lty=2)
lines(Amazon.rad, col="red")
##curve(exp(cf.p.lm[1]+cf.p.lm[2]*x), add=TRUE, 
##      col="red", lty=2)
##curve(exp(cf.p.lm[1]-d + cf.p.lm[2]*x), add=TRUE, 
##      col="red", lty=2)

