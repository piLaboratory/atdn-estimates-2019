## ----setup, echo=FALSE---------------------------------------------------
library(VGAM)
library(fitdistrplus)
library(sads)
library(knitr)
library(SPECIES)
source("functions.R")
opts_chunk$set(fig.align = 'center', fig.show = 'hold',
               fig.height = 6, warning = FALSE, message = FALSE,
               error = FALSE, echo=TRUE)
options(formatR.arrow = TRUE, width = 90, cache=TRUE)    
y <- read.table("abundVector.txt")[,1]
Sobs <- length(y)

## ----fit pln, cache=TRUE-------------------------------------------------
pln <- fitpoilog(y)
par(mfrow=c(2,2))
plot(pln)
par(mfrow=c(1,1))

## ----poilog coefficients-------------------------------------------------
(pln.cf <- coef(pln))

## ----mean lognormal------------------------------------------------------
pln.cf <- unname(pln.cf)
Mu <- pln.cf[1]-log(1.9e3/5.5e8)
(mean.ln <- exp(Mu + pln.cf[2]^2/2))

## ----lognorm S estimate--------------------------------------------------
## Total number of trees (average density x area)
Tot.t <- 567*5.5e8
## Estimated number of species
(pln.S <- Tot.t/mean.ln)

## ----PLN species richness estimate 2-------------------------------------
(pln.d0 <- dpoilog(0, mu = pln.cf[1], sig=pln.cf[2]))

## ----confindence and likelihood limits of PLN, cache=TRUE----------------
pln.prf <- profile(pln)
par(mfrow=c(1,2))
plotprofmle(pln.prf)
par(mfrow=c(1,1))
(pln.ci <- confint(pln.prf))
(pln.li <- likelregions(pln.prf))

## ----PLN richness ci, echo=FALSE, cache=TRUE-----------------------------
pln.smu <- seq(pln.li$mu[[1]][1], pln.li$mu[[1]][2], length=25)
pln.ssig <- seq(pln.li$sig[[1]][1], pln.li$sig[[1]][2], length=25)
LL1 <- Vectorize( function(mu,sig)
    -sum(dtrunc("poilog", x=y, trunc=0, coef=list(mu=mu, sig=sig), log=TRUE)) )
ER1 <- Vectorize(function(mu,sig) Sobs/ (1 - dpoilog(0, mu,sig)))
pln.ll <- outer(pln.smu, pln.ssig, LL1)
pln.llr <- pln.ll - LL1(pln.cf[1],pln.cf[2])
Sest <- outer(pln.smu, pln.ssig, ER1)
contour(x=pln.smu, y=pln.ssig, z=pln.llr, xlab="mu", ylab="sigma",
        levels=3:5)
contour(x=pln.smu, y=pln.ssig, z=pln.llr, levels=2, add=TRUE, lwd=2, labcex=1.5)
contour(x=pln.smu, y=pln.ssig, z=Sest, add=TRUE, col="blue", lwd=2, labcex=1,
        levels= c(5501,5650, round(Sobs/(1-pln.d0))))

## ----<pln S ci-----------------------------------------------------------
## Lower bound of S
Sobs / (1 - dpoilog(0, pln.ci[1,2], pln.ci[2,1]))
## Upper bound
Sobs / (1 - dpoilog(0, pln.ci[1,1], pln.ci[2,2]))

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

## ----fit NB, cache=TRUE--------------------------------------------------
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
(S.est <- unname(Sobs*(1-(1-csi)^cf.nb[1]) / (1-(1-csi.p)^cf.nb[1])))

## ----NB S est CI---------------------------------------------------------
tovo(fit = y.nb2, p = p1, CI=TRUE)

## ----model selection-----------------------------------------------------
AICtab(pln, y.nb2, y.ls)

## ----beta fit, echo=FALSE, eval=FALSE------------------------------------
## dados <- read.csv2("data.csv", as.is=TRUE)
## ########################################################
## ## Fit a beta and truncated beta to observed proportions
## ########################################################
## ## Initial fit to get start values
## f.st <- fitdist(dados$N.plots/1900, "beta")
## cf.st <- coef(f.st)
## ## likelihood functions
## LL.beta <- function(s1, s2){
##     -sum ( dbeta(x = dados$N.plots/1900, shape1=s1, shape2=s2, log=TRUE) )
## }
## LL.betat.1 <- function(s1, s2){
##     -sum ( sads::dtrunc("beta", x = dados$N.plots/1900, trunc=0.1/1900, coef = list(shape1=s1, shape2=s2), log=TRUE) )
## }
## LL.betat.5 <- function(s1, s2){
##     -sum ( sads::dtrunc("beta", x = dados$N.plots/1900, trunc=0.5/1900, coef = list(shape1=s1, shape2=s2), log=TRUE) )
## }
## LL.betat.9 <- function(s1, s2){
##     -sum ( sads::dtrunc("beta", x = dados$N.plots/1900, trunc=0.9/1900, coef = list(shape1=s1, shape2=s2), log=TRUE) )
## }
## ## Fits: ran first with SANN and the with default method
## f.beta <- mle2(LL.beta, start=list(s1=cf.st[1], s2=cf.st[2]))
## f.betat.1 <- mle2(LL.betat.1, start=list(s1=cf.st[1], s2=cf.st[2]),
##                   method="SANN")
## f.betat.1 <- mle2(LL.betat.1, start=as.list(coef(f.betat.1)))
## f.betat.5 <- mle2(LL.betat.5, start=list(s1=cf.st[1], s2=cf.st[2]),
##                   method="SANN")
## f.betat.5 <- mle2(LL.betat.5, start=as.list(coef(f.betat.5)))
## f.betat.9 <- mle2(LL.betat.9, start=list(s1=cf.st[1], s2=cf.st[2]),
##                   method="SANN")
## f.betat.9 <- mle2(LL.betat.9, start=as.list(coef(f.betat.9))) ## did not converge
## ## Coefficients
## (cf.beta <- coef(f.beta))
## (cf.betat.1 <- coef(f.betat.1))
## (cf.betat.5 <- coef(f.betat.5))
## (cf.betat.9 <- coef(f.betat.9))

## ----beta qqplots, echo=FALSE--------------------------------------------
## Checking fits with QQ-plots
par(mfrow=c(2,2))
QQ.plot(dados$N.plots/1900, distr = "beta",
        coef = list(shape1 = cf.beta[1], shape2 = cf.beta[2]),
        main="No Truncation")
QQ.plot(dados$N.plots/1900, distr = "beta", trunc=0.1/1900,
        coef = list(shape1 = cf.betat.1[1], shape2 = cf.betat.1[2]),
        main="Truncation at 0.1/1900")
QQ.plot(dados$N.plots/1900, distr = "beta", trunc=0.5/1900,
        coef = list(shape1 = cf.betat.5[1], shape2 = cf.betat.5[2]),
        main="Truncation at 0.5/1900")
QQ.plot(dados$N.plots/1900, distr = "beta", trunc=0.9/1900,
        coef = list(shape1 = cf.betat.9[1], shape2 = cf.betat.9[2]),
        main="Truncation at 0.9/1900")
par(mfrow=c(1,1))

## ----beta model selection------------------------------------------------
AICtab(f.beta, f.betat.1, f.betat.5, f.betat.9)

## ----Estimated number of species, echo=FALSE-----------------------------
## Estimated number of species using the beta binom truncated
(bb.S <- bb.Sest(shape1=cf.beta[1], shape2=cf.beta[2], Sobs=S.obs, size=1900))
(bbt.1.S <- bb.Sest(shape1=cf.betat.1[1], shape2=cf.betat.1[2], Sobs=S.obs, size=1900))
(bbt.5.S <- bb.Sest(shape1=cf.betat.5[1], shape2=cf.betat.5[2], Sobs=S.obs, size=1900))
(bbt.9.S <- bb.Sest(shape1=cf.betat.9[1], shape2=cf.betat.9[2], Sobs=S.obs, size=1900)) ## unrealistic high

## ----table of estimated richness, echo=FALSE-----------------------------
kable(
    data.frame(
        trunc=c("None","0.1/1900","0.5/1900","0.9/900"),
        S = c(bb.S, bbt.1.S, bbt.5.S, bbt.9.S)),
    digits=0, col.names=c("Truncation","Estimated richness")
    )

## ----beta likelihood surface, echo=FALSE, cache=TRUE---------------------
## Profile
betat.5.prf <- profile(f.betat.5)
betat.5.ci <- confint(betat.5.prf)
# plotprofmle(betat.5.prf)
## Likelihood surface
## Shows that likelihood interval of species richness is ca 9880 - 10280
x <- seq(betat.5.ci[1,1]*.95, betat.5.ci[1,2]*1.05, length=30)
y <- seq(betat.5.ci[2,1]*.95, betat.5.ci[2,2]*1.05, length=30)
Vbb.Sest <- Vectorize(function(s1, s2,...) bb.Sest(shape1=s1, shape2=s2, ...), c("s1", "s2"))
VLL.betat.5 <- Vectorize(LL.betat.5)
z1 <- outer(x, y, VLL.betat.5)
z1 <- z1 - min(z1)
z2 <- outer(x, y, Vbb.Sest, Sobs=S.obs, size=1900)
contour(x,y,z1, levels=2, col="red", xlab="shape1", ylab="shape2")
contour(x, y, z2, add=TRUE, col="blue", levels=c(12550,20800) )
contour(x, y, z2, add=TRUE, col="blue", levels=round(bbt.5.S))

## ----fit truncated beta-binomial, cache=TRUE, echo=FALSE, eval=FALSE-----
## ##start values converted for rho and prob,
## ## both bounded between zero and one
## ## See help(VGAM::betabinom)
## mu.st <- cf.st[1]/sum(cf.st)
## rho.st <- 1/(1 + sum(cf.st))
## ### With logit of parameters
## ## Likelihood function
## unlogit <- function(x) exp(x)/(1+exp(x))
## logit <- function(x) log(x/(1-x))
## LL2.bbt0 <- function(lmu, lrho){
##     -sum(dbetabinom.t0(x = dados$N.plots, size = 1900, mu = unlogit(lmu), rho = unlogit(lrho), log=TRUE))
## }
## ## fit with mle2
## beta.bin.ML2 <- bbmle::mle2(LL2.bbt0, start=list(lmu=logit(mu.st), lrho=logit(rho.st)), method="SANN")
## beta.bin.ML2 <- bbmle::mle2(LL2.bbt0, start=as.list(coef(beta.bin.ML2)))
## bb.mle2 <- unlogit(coef(beta.bin.ML2))
## beta.bin.prf2 <- profile(beta.bin.ML2)
## bb.ci2 <- confint(beta.bin.prf2)
## ## Total number of species
## bb2.S <- bb.Sest(mu = bb.mle2[1], rho = bb.mle2[2], size = 1900, Sobs=S.obs)

## ----beta-binomial likelihood surface, echo=FALSE------------------------
x <- seq(unlogit(bb.ci2[1,1]), unlogit(bb.ci2[1,2]), length=30)
y <- seq(unlogit(bb.ci2[2,1]), unlogit(bb.ci2[2,2]), length=30)
Vbb.Sest <- Vectorize(bb.Sest, c("mu", "rho"))
VLL.bbt0 <- Vectorize(LL.bbt0)
z1 <- outer(x, y, VLL.bbt0)
z1 <- z1 - min(z1)
z2 <- outer(x, y, Vbb.Sest, Sobs=S.obs, size=1900)
contour(x,y,z1, levels=2, col="red", xlab="mu", ylab="rho")
contour(x, y, z2, add=TRUE, col="blue", levels=c(9878,10280)) 
contour(x, y, z2, add=TRUE, col="blue", levels=round(bb2.S))

## ----SPECIES package, cache=TRUE, eval=TRUE------------------------------
Ab <- as.data.frame(table(dados$N.ind))
names(Ab) <- c("j", "n_j")
Ab$j <- as.integer(Ab$j)
##jackknife method
(Sj <- jackknife(Ab,k=5))
##ACE coverage method
(SChao92 <- ChaoLee1992(Ab,t=10, method="all"))
##Chao1984 lower bound estimator
(SChao84 <- chao1984(Ab))
##Chao and Bunge coverage-duplication method
(SChaoB <- ChaoBunge(Ab,t=10))
##Penalized NPMLE method
(SNPMLE <- pnpmle(Ab,t=15))
##Unconditonal NPMLE method
(SUNPMLE <-unpmle(Ab,t=10))
##Poisson-compound Gamma method
(SPG <- pcg(Ab,t=20))

