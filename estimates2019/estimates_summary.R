## ----setup, echo=FALSE, include=FALSE------------------------------------
library(VGAM)
library(untb)
#library(fitdistrplus)
library(sads)
library(knitr)
#library(SPECIES)
library(xtable)
library(abc)
source("functions.R")
opts_chunk$set(fig.align = 'center', fig.show = 'hold',
               fig.height = 6.5, warning = FALSE, message = FALSE,
               error = FALSE, echo=TRUE)
options(formatR.arrow = TRUE, width = 90, cache=TRUE, scipen = 1, digits = 2)

## ----Data prep-----------------------------------------------------------
atdn.2019 <- read.csv2("Populations_2019_V1.csv", as.is=TRUE)
N.ind <- atdn.2019$N.ind
Sobs <- length(N.ind)
## Total number of trees (average density x area)
## Tot.t <- 567*5.5e8
Tot.t <- 3.06e11
## Proportion of total trees in the sample
p1 <- sum(atdn.2019$N.ind)/Tot.t
## Total number of plots
N.plots <- 1946
## Total area hectares
#Tot.A <- 5.79e8
Tot.A <- 6.29e8
## Sampled area ha
Samp.A <- 2.038e3

## ----fit pln, cache=TRUE-------------------------------------------------
pln <- fitpoilog(N.ind)
par(mfrow=c(2,2))
plot(pln)
par(mfrow=c(1,1))

## ----PLN species richness estimate 2-------------------------------------
pln.cf <- coef(pln)
(pln.d0 <- dpoilog(0, mu = pln.cf[1], sig=pln.cf[2]))

## ----fit ls--------------------------------------------------------------
y.ls <- fitls(N.ind)
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
y.nb <- vglm(N.ind ~ 1, posnegbinomial)
## With sads
y.nb2 <- fitnbinom(N.ind, 
                   start.value=c(size=0.03, mu=mean(N.ind)))
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
par(mfrow=c(1,2))
ls.rad <- radpred(y.ls)
nb.rad <- radpred(y.nb2)
## plot(ls.rad$abund, nb.rad$abund, log="xy",
##      xlab="Logseries theor. quantiles",
##      ylab="TNB theor. quantiles", 
##      col="grey", cex=0.25)
## abline(0,1, col="blue")
plot(rad(atdn.2019$N.ind), col="grey", cex=0.5)
lines(ls.rad)
lines(nb.rad, col="red")
legend("topright", c("LS", "NB"), col=c("blue","red"), 
       lty=1, bty="n")
ls.oc <- octavpred(y.ls)
nb.oc <- octavpred(y.nb2)
plot(octav(atdn.2019$N.ind), ylim=range(c(ls.oc$Freq, nb.oc$Freq)))
lines(ls.oc)
lines(nb.oc, col="red")
legend("topright", c("LS", "NB"), col=c("blue","red"), 
       lty=1, bty="n")
par(mfrow=c(1,1))

## ----Linear extrapolation from regional RAD------------------------------
S.ulrich <- ulrich(atdn.2019$population)
(S.r.ls <- S.ulrich$S[1])

## ----Amazon alpha--------------------------------------------------------
(alpha.r <- fishers.alpha(N = Tot.t, S = S.r.ls))

## ----amazon LS rad-------------------------------------------------------
reg.ls.rad <- ceiling(
    rad.ls(S = S.r.ls, N = Tot.t, alpha = alpha.r)$y
)

## ----regional rad lognormal----------------------------------------------
(S.ulrich$S[2])

## ----TNB regionl RAD-----------------------------------------------------
reg.nb.rad <- rad.posnegbin(S = S.nb, size = cf.nb[1], 
                            prob = 1-csi)$y

## ----nbinom rad, echo=FALSE----------------------------------------------
cf.u <- S.ulrich$coefs
plot(rad(reg.ls.rad), col="red", lwd=2, type="n",
     ylim=c(1, max(atdn.2019$population)))
points(rad(atdn.2019$population), col="grey")
curve(exp(cf.u[1]+cf.u[2]*x), add=TRUE)
curve(exp(cf.u[1]-cf.u[3]+cf.u[2]*x), add=TRUE)
lines(rad(reg.nb.rad), col="blue", lwd=2)
lines(rad(reg.ls.rad), col="red", lwd=2)
legend("topright", c("LS", "TNB", "Linear upper/lower bounds"), 
       col=c("red", "blue", "black"), lty=1, lwd=2, bty="n")

## ----k x dens regression, cache=TRUE-------------------------------------
## estimating k parameter of a NB for each species 
atdn.2019$dens.ha <- atdn.2019$N.ind/Samp.A
atdn.2019$k <- est.kv(mu=atdn.2019$dens.ha, 
                  nzeroes=N.plots-atdn.2019$N.plots, 
                  Nplots=N.plots)
lm.k <-lm(log(k)~log(dens.ha), 
          data=atdn.2019, subset=k<1)
## Estimated regression standard error
lm.k.sigma <- summary(lm.k)$sigma
## Model summary
summary(lm.k)

## ----lok k x log density plot, echo=FALSE--------------------------------
plot(log(k)~log(dens.ha), data=atdn.2019, 
     subset=k<1, xlab="Density (ind/ha)",
     ylab="Aggregation parameter of NB",
     col="grey")
abline(lm.k, col="blue", lwd=2)

## ----loads ABC simulation data-------------------------------------------
load("abc_simulations/abc2019.RData")

## ----abc model selection-------------------------------------------------
## Target: observed number of species, Simpson's 1/D, 
## lmean, sdmean             
target <- c(Sobs, 
            D(atdn.2019$population), 
            mean(log(atdn.2019$population)), 
            sd(log(atdn.2019$population)))
## Model selection
model.sel <- postpr(target = target,
                    index=abc2019$labels,
                    sumstat = abc2019$sims,
                    tol=0.025, method="rejection",
                    corr=TRUE)
msel.s <- summary(model.sel)

## ----posterior species richness------------------------------------------
## Posterior distribution of Species richness from the selected model
S.post1 <- abc(target = target, param=data.frame(S=abc2019$params[abc2019$labels=="LSclump"]),
              sumstat = abc2019$sims[abc2019$labels=="LSclump",],
              tol=0.025, method="rejection")
S.post1.s <- summary(S.post1)
hist(S.post1)

## ----clumped LS RAD and CI, cache=TRUE-----------------------------------
## Predicted log(k) values for LS rad
reg.ls.rad.lk <- predict(lm.k, 
                         newdata=data.frame(dens.ha=reg.ls.rad/Tot.A))

## LS RAD with the lower CI
abc.ls.c.rad.l <- rad.ls(S = S.post1.s[2,],
                           N = Tot.t, 
                           alpha = fishers.alpha(Tot.t, S.post1.s[2,]))$y
## LS RAD with the upper CI
abc.ls.c.rad.u <- rad.ls(S = S.post1.s[6,],
                           N = Tot.t, 
                           alpha = fishers.alpha(Tot.t, S.post1.s[6,]))$y
## Simulated clumped samples
## from LS with species richness estimated from linear extrapolation
ls.clump <- NB.samp(rad = reg.ls.rad, tot.area = Tot.A, 
                    n.plots = N.plots, 
                    lmean.k = reg.ls.rad.lk, 
                    lsd.k = lm.k.sigma, 
                    nrep=100)
## From LS with richeness from posterior cerdible intervals
ls.s2 <- NB.samp(rad = abc.ls.c.rad.l, 
                 tot.area = Tot.A, 
                 n.plots = N.plots, 
                 lmean.k =  
                     predict(lm.k, 
                             newdata=data.frame(dens.ha=abc.ls.c.rad.l/Tot.A)),
                 lsd.k = lm.k.sigma, 
                 nrep=100)
ls.s3 <- NB.samp(rad = abc.ls.c.rad.u, 
                 tot.area = Tot.A, 
                 n.plots = N.plots, 
                 lmean.k =  
                     predict(lm.k, 
                             newdata=data.frame(dens.ha=abc.ls.c.rad.u/Tot.A)),
                 lsd.k = lm.k.sigma, 
                 nrep=100)
## Plot
plot(rad(atdn.2019$population), col="grey",
     ylab = "Population size", xlim=c(1,sum(ls.s3>=1)))

lines(rad(ls.clump), lwd=2, col="blue")
lines(rad(ls.s2), lwd=2, col="blue", lty=2)
lines(rad(ls.s3), lwd=2, col="blue", lty=2)

## ----sim popsizes LS and NB, cache=TRUE----------------------------------
## Predicted log(k) values for LS rad
reg.ls.rad.lk <- predict(lm.k, 
                         newdata=data.frame(dens.ha=reg.ls.rad/Tot.A))
## Predicted log(k) values for TNB rad
reg.nb.rad.lk <- predict(lm.k, 
                           newdata=data.frame(dens.ha=reg.nb.rad/Tot.A))
## Simulation of population sizes from samples of each RAD
ls.rnd <- Pois.samp(rad = reg.ls.rad, tot.area = Tot.A, 
                    n.plots = N.plots, nrep=100)
nb.rnd <- Pois.samp(rad = reg.nb.rad, tot.area = Tot.A, 
                    n.plots = N.plots, nrep = 100)
nb.clump <- NB.samp(rad = reg.nb.rad, tot.area = Tot.A, 
                    n.plots = N.plots, 
                    lmean.k = reg.nb.rad.lk, 
                    lsd.k = lm.k.sigma, 
                    nrep = 100)

## ----plot sampled RADS, echo=FALSE---------------------------------------
plot(rad(atdn.2019$population), col="grey", log="y", xlim=c(1,5500),
     ylab = "Population size")
lines(rad(ls.rnd), lwd=2)
lines(rad(ls.clump), lwd=2, col="red")
lines(rad(nb.rnd), lwd=2, col="green")
lines(rad(nb.clump), lwd=2, col="black")
legend("topright", 
       c("LS, random sample", "LS, clumped sample", "NB, random sample", "NB, clumped sample"), 
       col=c("blue","red","green", "black"), lty=1, lwd=2, bty="n")

## ----shen estimate, eval=FALSE-------------------------------------------
## ## table of frequencies of occurrences
## Y <- data.frame(table(atdn.2019$N.plots))
## Y[,1] <- as.integer(as.character(Y[,1]))
## ## Estimate of alfa and beta (Eq.6)
## ## To be use as starting values for the unconditional estimation below
## ab.est <- shen.ab(Y = Y, t = N.plots, T = Tot.A,
##                 start=list(lalpha=-5, lbeta=3), method="SANN")
## ## estimated coeficients
## cf.st1 <- coef(ab.est)
## ## Estimate with uncoditional likelihood (Eq.3)
## ## restricted to species richness between 1e4 and 2e4
## ShenHe <- shen.S( Y = Y, t = N.plots, T = Tot.A,
##                 start=c(list(lS = log(1.5e4)), as.list(cf.st1)),
##                 method="L-BFGS-B",
##                 upper=c(lS=log(2e4), lalpha=Inf, lbeta=Inf),
##                 lower=c(lS=log(1e3), lalpha=-Inf, lbeta=-Inf))
## cf.st2 <- coef(ShenHe)
## (S.sh <- unname(exp(cf.st2[1])))

## ----Shen He confint-----------------------------------------------------
ShenHe.prf <- profile(ShenHe, which=1)
exp(confint(ShenHe.prf))

## ----Hui ORC estimate----------------------------------------------------
S.orc <- hui.orc(atdn.2019$N.plots, effort=Samp.A/Tot.A)
orc.cf <- coef(S.orc$model)

## ----ORC plots with Hui function-----------------------------------------
x <- 1:nrow(atdn.2019)
y <- exp(orc.cf[1])*exp(orc.cf[2]*x)*(x^orc.cf[3])
#plot(rad(atdn.2019$N.plots), ylim=range(c(atdn.2019$N.plots, y)))
#lines(rad(y))
x2 <- 1:S.orc$S.est
y2 <- exp(orc.cf[1]-log(eff18))*exp(orc.cf[2]*x2)*(x2^orc.cf[3])
plot(rad(y2), type="n")
points(rad(atdn.2019$N.plots/eff18), col="grey")
lines(rad(y2), type="l")

