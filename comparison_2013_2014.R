## ----setup, echo=FALSE, include=FALSE------------------------------------
library(VGAM)
library(untb)
library(sads)
library(abc)
source("../functions.R")

## ----Data prep-----------------------------------------------------------
atdn.2013 <- read.csv2("estimates2013/science_appendix2013.csv", as.is=TRUE)
atdn.2019 <- read.csv2("estimates2019/Populations_2019_V1.csv", as.is=TRUE)
N.ind.13 <- atdn.2013$N.ind
Sobs.13 <- length(N.ind.13)
## Total number of trees (average density x area)
Tot.t.13 <- 3.9e11
## Proportion of total trees in the sample
p1.13 <- sum(atdn.2013$N.ind)/Tot.t.13
## Total number of plots
N.plots.13 <- 1170
## Total area hectares
#Tot.A <- 5.79e8
Tot.A <- 6.29e8
## Sampled area ha
Samp.A.13 <- 1.10931e3

## ----fit pln, cache=TRUE-------------------------------------------------
pln.13 <- fitpoilog(N.ind.13)

## ----PLN species richness estimate 2-------------------------------------
pln.13.cf <- coef(pln.13)
(pln.13.d0 <- dpoilog(0, mu = pln.13.cf[1], sig=pln.13.cf[2]))

## ----fit ls--------------------------------------------------------------
ls.13 <- fitls(N.ind.13)

## ----estimated S ls------------------------------------------------------
ls.alpha.13 <- coef(ls.13)[[2]]
(S.ls.13 <- ls.alpha.13*log(1 + Tot.t.13/ls.alpha.13))

## ----ls and S est for ls-------------------------------------------------
(ls.ci.13 <- confint(ls.13))
## Estimated species richness for lower bound of alpha's IC
ls.ci.13[1]*log(1 + Tot.t.13/ls.ci.13[1])
ls.ci.13[2]*log(1 + Tot.t.13/ls.ci.13[2])

## ----fit NB, cache=TRUE--------------------------------------------------
## With sads
tnb.13 <- fitnbinom(N.ind.13, 
                   start.value=c(size=0.3, mu=mean(N.ind.13)))

## ----Tovo S estimate-----------------------------------------------------
tnb.cf.13 <- coef(tnb.13)
csi.p.13 <- unname(tnb.cf.13[2]/(sum(tnb.cf.13)))
csi.13 <- csi.p.13/(p1.13+(1-p1.13)*csi.p.13)
## Estimated number of species 
S.nb.13 <- Sobs.13*(1-(1-csi.13)^tnb.cf.13[1]) / (1-(1-csi.p.13)^tnb.cf.13[1])
(S.nb.13 <- unname(S.nb.13))

## ----NB S est CI---------------------------------------------------------
(S.tovo.13 <- tovo(fit = tnb.13, p = p1.13, CI=TRUE))

## ----model selection-----------------------------------------------------
AICtab(pln.13, tnb.13, ls.13, base=TRUE)


## ----Linear extrapolation from regional RAD------------------------------
S.ulrich.13 <- ulrich(atdn.2013$population)
(S.r.ls.13 <- S.ulrich.13$S[1])

## ----Amazon alpha--------------------------------------------------------
(ls.alpha.reg.13 <- fishers.alpha(N = Tot.t.13, S = S.r.ls.13))

## ----amazon LS rad-------------------------------------------------------
reg.rad.ls.13 <- ceiling(
    rad.ls(S = S.r.ls.13, N = Tot.t.13, alpha = ls.alpha.reg.13)$y
)

## ----regional rad lognormal----------------------------------------------
(S.ulrich.13$S[2])

## ----TNB regionl RAD-----------------------------------------------------
reg.rad.tnb.13 <- rad.posnegbin(S = S.nb.13, size = tnb.cf.13[1], 
                            prob = 1-csi.13)$y

## ----k x dens regression, cache=TRUE-------------------------------------
## estimating k parameter of a NB for each species 
atdn.2013$dens.ha <- atdn.2013$N.ind/Samp.A.13
atdn.2013$k <- est.kv(mu=atdn.2013$dens.ha, 
                  nzeroes=N.plots.13-atdn.2013$N.plots, 
                  Nplots=N.plots.13)
lm.k.13 <-lm(log(k)~log(dens.ha), 
          data=atdn.2013)
## Estimated regression standard error
lm.k.sigma.13 <- summary(lm.k)$sigma
## Model summary
summary(lm.k.13)

## ----lok k x log density plot, echo=FALSE--------------------------------
plot(log(k)~log(dens.ha), data=atdn.2013, 
     #subset=k<1,
     xlab="Density (ind/ha)",
     ylab="Aggregation parameter of NB",
     col="grey")
abline(lm.k.13, col="blue", lwd=2)

## ----loads ABC simulation data-------------------------------------------
load("abc_simulations/abc2013.RData")

## ----abc model selection-------------------------------------------------
## Target: observed number of species, Simpson's 1/D, 
## lmean, sdmean             
target <- c(Sobs.13, 
            D(atdn.2013$population), 
            mean(log(atdn.2013$population)), 
            sd(log(atdn.2013$population)))
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
reg.rad.ls.13.lk <- predict(lm.k.13, 
                         newdata=data.frame(dens.ha=reg.rad.ls.13/Tot.A))

## LS RAD with the lower CI
abc.ls.c.rad.l <- rad.ls(S = S.post1.s[2,],
                           N = Tot.t.13, 
                           alpha = fishers.alpha(Tot.t.13, S.post1.s[2,]))$y
## LS RAD with the upper CI
abc.ls.c.rad.u <- rad.ls(S = S.post1.s[6,],
                           N = Tot.t.13, 
                           alpha = fishers.alpha(Tot.t.13, S.post1.s[6,]))$y
## Simulated clumped samples
## from LS with species richness estimated from linear extrapolation
ls.clump <- NB.samp(rad = reg.rad.ls.13, tot.area = Tot.A, 
                    n.plots = N.plots.13, 
                    lmean.k = reg.rad.ls.13.lk, 
                    lsd.k = lm.k.sigma.13, 
                    nrep=100)
## From LS with richeness from posterior cerdible intervals
ls.s2 <- NB.samp(rad = abc.ls.c.rad.l, 
                 tot.area = Tot.A, 
                 n.plots = N.plots.13, 
                 lmean.k =  
                     predict(lm.k.13, 
                             newdata=data.frame(dens.ha=abc.ls.c.rad.l/Tot.A)),
                 lsd.k = lm.k.sigma.13, 
                 nrep=100)
ls.s3 <- NB.samp(rad = abc.ls.c.rad.u, 
                 tot.area = Tot.A, 
                 n.plots = N.plots.13, 
                 lmean.k =  
                     predict(lm.k.13, 
                             newdata=data.frame(dens.ha=abc.ls.c.rad.u/Tot.A)),
                 lsd.k = lm.k.sigma.13, 
                 nrep=100)
## Plot
plot(rad(atdn.2013$population), col="grey",
     ylab = "Population size", xlim=c(1,sum(ls.s3>=1)))

lines(rad(ls.clump), lwd=2, col="blue")
lines(rad(ls.s2), lwd=2, col="blue", lty=2)
lines(rad(ls.s3), lwd=2, col="blue", lty=2)

## ----sim popsizes LS and NB, cache=TRUE----------------------------------
## Predicted log(k) values for LS rad
reg.rad.ls.13.lk <- predict(lm.k.13, 
                         newdata=data.frame(dens.ha=reg.rad.ls.13/Tot.A))
## Predicted log(k) values for TNB rad
reg.rad.tnb.13.lk <- predict(lm.k.13, 
                             newdata=data.frame(dens.ha=reg.rad.tnb.13/Tot.A))
### PAREI aqui de renomear objetos ##

## Simulation of population sizes from samples of each RAD
ls.rnd <- Pois.samp(rad = reg.rad.ls.13, tot.area = Tot.A, 
                    n.plots = N.plots.13, nrep=100)
nb.rnd <- Pois.samp(rad = reg.rad.tnb.13, tot.area = Tot.A, 
                    n.plots = N.plots.13, nrep = 100)
nb.clump <- NB.samp(rad = reg.rad.tnb.13, tot.area = Tot.A, 
                    n.plots = N.plots.13, 
                    lmean.k = reg.rad.tnb.13.lk, 
                    lsd.k = lm.k.sigma.13, 
                    nrep = 100)

## ----plot sampled RADS, echo=FALSE---------------------------------------
plot(rad(atdn.2013$population), col="grey", log="y", xlim=c(1,5500),
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
Y <- data.frame(table(atdn.2013$N.plots))
Y[,1] <- as.integer(as.character(Y[,1]))
## Estimate of alfa and beta (Eq.6)
## To be use as starting values for the unconditional estimation below
ab.est <- shen.ab(Y = Y, t = N.plots.13, T = Tot.A,
                start=list(lalpha=-5, lbeta=3), method="SANN")
## estimated coeficients
cf.st1 <- coef(ab.est)
## Estimate with uncoditional likelihood (Eq.3)
## restricted to species richness between 1e4 and 2e4
ShenHe <- shen.S( Y = Y, t = N.plots.13, T = Tot.A,
                start=c(list(lS = log(1.5e4)), as.list(cf.st1)),
                method="L-BFGS-B",
                upper=c(lS=log(2e4), lalpha=Inf, lbeta=Inf),
                lower=c(lS=log(1e3), lalpha=-Inf, lbeta=-Inf))
cf.st2 <- coef(ShenHe)
(S.sh <- unname(exp(cf.st2[1])))

## ----Shen He confint-----------------------------------------------------
ShenHe.prf <- profile(ShenHe, which=1) # not working
exp(confint(ShenHe.prf))

## ----Hui ORC estimate----------------------------------------------------
S.orc <- hui.orc(atdn.2013$N.plots, effort=Samp.A.13/Tot.A)
orc.cf <- coef(S.orc$model)
S.orc$S.est
## ----ORC plots with Hui function-----------------------------------------
x <- 1:nrow(atdn.2013)
y <- exp(orc.cf[1])*exp(orc.cf[2]*x)*(x^orc.cf[3])
##plot(rad(atdn.2013$N.plots), ylim=range(c(atdn.2013$N.plots, y)))
##lines(rad(y))
eff13 <- Samp.A.13/Tot.A
x2 <- 1:S.orc$S.est
y2 <- exp(orc.cf[1]-log(eff13))*exp(orc.cf[2]*x2)*(x2^orc.cf[3])
plot(rad(y2), type="n")
points(rad(atdn.2013$N.plots/eff13), col="grey")
lines(rad(y2), type="l")

