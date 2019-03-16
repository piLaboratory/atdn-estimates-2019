################################################################################
## Following instructions at https://github.com/stan-dev/rstan/wiki/Installing-RStan-on-Mac-or-Linux
## To set compiler to install rstan in linux, check  if needed in oter OS (runned once)
dotR <- file.path(Sys.getenv("HOME"), ".R")
if (!file.exists(dotR)) dir.create(dotR)
M <- file.path(dotR, "Makevars")
if (!file.exists(M)) file.create(M)
cat("\nCXXFLAGS=-O3 -mtune=native -march=native -Wno-unused-variable -Wno-unused-function", 
    file = M, sep = "\n", append = TRUE)
cat("\nCXXFLAGS+=-flto -ffat-lto-objects  -Wno-unused-local-typedefs", 
    file = M, sep = "\n", append = TRUE)
################################################################################

library(rstan)
library(parallel)
library(fitdistrplus)
library(truncdist)
library(sads)
library(VGAM)
source("functions.R")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## Total number of species estimated from a Bayesian beta-binomial model
dados <- read.csv2("data.csv", as.is=TRUE)

################################################################################
## Fit a beta and truncated beta to observed proportions
################################################################################
## Fit by mle
f.st <- fitdist(dados$N.plots/1900, "beta")
cf.st <- coef(f.st)
## beta truncated at x<0.5/1900
## likelihood functions
LL.beta <- function(s1, s2){
    -sum ( dbeta(x = dados$N.plots/1900, shape1=s1, shape2=s2, log=TRUE) )
}
LL.betat.1 <- function(s1, s2){
    -sum ( dtrunc("beta", x = dados$N.plots/1900, trunc=0.1/1900, coef = list(shape1=s1, shape2=s2), log=TRUE) )
}
LL.betat.5 <- function(s1, s2){
    -sum ( dtrunc("beta", x = dados$N.plots/1900, trunc=0.5/1900, coef = list(shape1=s1, shape2=s2), log=TRUE) )
}
LL.betat.9 <- function(s1, s2){
    -sum ( dtrunc("beta", x = dados$N.plots/1900, trunc=0.9/1900, coef = list(shape1=s1, shape2=s2), log=TRUE) )
}
f.beta <- mle2(LL.beta, start=list(s1=cf.st[1], s2=cf.st[2]))
f.betat.1 <- mle2(LL.betat.1, start=list(s1=cf.st[1], s2=cf.st[2]), method="SANN")
f.betat.1 <- mle2(LL.betat.1, start=as.list(coef(f.betat.1)))
f.betat.5 <- mle2(LL.betat.5, start=list(s1=cf.st[1], s2=cf.st[2]), method="SANN")
f.betat.5 <- mle2(LL.betat.5, start=as.list(coef(f.betat.5)))
f.betat.9 <- mle2(LL.betat.9, start=list(s1=cf.st[1], s2=cf.st[2]), method="SANN")
f.betat.9 <- mle2(LL.betat.9, start=as.list(coef(f.betat.9))) ## only with SANN
## Model selection
AICtab(f.beta, f.betat.1, f.betat.5, f.betat.9)
## Coefficients
(cf.beta <- coef(f.beta))
(cf.betat.1 <- coef(f.betat.1))
(cf.betat.5 <- coef(f.betat.5))
(cf.betat.9 <- coef(f.betat.9))
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

## Estimated number of species using the beta binom truncated
(bb.S <- bb.Sest(shape1=cf.beta[1], shape2=cf.beta[2], Sobs=S.obs, size=1900))
(bbt.1.S <- bb.Sest(shape1=cf.betat.1[1], shape2=cf.betat.1[2], Sobs=S.obs, size=1900))
(bbt.5.S <- bb.Sest(shape1=cf.betat.5[1], shape2=cf.betat.5[2], Sobs=S.obs, size=1900))
(bbt.9.S <- bb.Sest(shape1=cf.betat.9[1], shape2=cf.betat.9[2], Sobs=S.obs, size=1900)) ## unrealistic high

## Profile
betat.5.prf <- profile(f.betat.5)
betat.5.ci <- confint(betat.5.prf)
plotprofmle(betat.5.prf)
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


################################################################################
## Fit truncated beta-binomial with maximum likelihood
################################################################################
##start values converted for rho and prob, both bounded between zero and one
## See help(VGAM::betabinom)
mu.st <- cf.st[1]/sum(cf.st)
rho.st <- 1/(1 + sum(cf.st))


## Likelihood function
LL.bbt0 <- function(mu, rho){
    -sum(dbetabinom.t0(x = dados$N.plots, size = 1900, mu = mu, rho = rho, log=TRUE))
}
## fit with mle2
beta.bin.ML <- mle2(LL.bbt0, start=list(mu=mu.st, rho=rho.st),
                    method="L-BFGS-B", lower= c(mu=1e-9, rho=1e-9),
                    upper = c(mu=1-1e-9, rho=1-1e-9))

beta.bin.prf <- profile(beta.bin.ML) ## returned a better fit and then does not the calculation
beta.bin.ML <- beta.bin.prf
summary(beta.bin.ML)
bb.mle <- coef(beta.bin.ML)
## Estimated total species
(bb.p0 <- bb.Sest(mu = bb.mle[1], rho = bb.mle[2], size = 1900, Sobs=S.obs))

## checking
mu <- bb.mle[1]
rho <- bb.mle[2]
## Backtransform in shape 1 and shape 2
bb.mle.s1 <- -(mu*rho-mu)/rho
bb.mle.s2 <- ((mu-1)*rho-mu+1)/rho
sum( 1- (1-rbeta(bb.p0, bb.mle.s1, bb.mle.s2))^1900) ## close, as expected

### With logit of parameters
## Likelihood function
unlogit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))
LL2.bbt0 <- function(lmu, lrho){
    -sum(dbetabinom.t0(x = dados$N.plots, size = 1900, mu = unlogit(lmu), rho = unlogit(lrho), log=TRUE))
}
## fit with mle2
beta.bin.ML2 <- bbmle::mle2(LL2.bbt0, start=list(lmu=logit(mu.st), lrho=logit(rho.st)), method="SANN")
beta.bin.ML2 <- bbmle::mle2(LL2.bbt0, start=as.list(coef(beta.bin.ML2))) ## ok
bb.mle2 <- unlogit(coef(beta.bin.ML2))
beta.bin.prf2 <- profile(beta.bin.ML2)
par(mfrow=c(1,2))
plotprofmle(beta.bin.prf2)
par(mfrow=c(1,1))
(bb.ci2 <- confint(beta.bin.prf2))
## Total number of species
(bb2.S <- bb.Sest(mu = bb.mle2[1], rho = bb.mle2[2], size = 1900, Sobs=S.obs))
## likelihood surface
## Shows that likelihood interval of species richness is ca 9880 - 10280
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

############################################################################################
## Fits zero-truncated beta binomial distribution in a explicit hierarchical Bayesian model 
############################################################################################
## Explicit hierarchical: compound distribution of a truncated binomial with
## values of parameter p from a truncated beta distribution

## Bayesian fit with Stan: estimated richness considerably low than the ML fit
## List with all data needed for the model
m1.data <- list(Nsites=1900, Sobs = nrow(dados), y = dados$N.plots)
fit.betabinom.1 <- stan(file="betabinom.stan", data=m1.data, iter=200, chains=3, cores=3)
save.image()

## Diagnostics
traceplot(fit.betabinom.1, c("alpha", "beta", "phi", "lambda", "Sest")) 
pairs(fit.betabinom, pars=c("alpha","beta"))
summary(fit.betabinom, c("alpha","beta"))$summary 
print(fit.betabinom.1, c("alpha","beta", "Sest"))
#sampler_params <- get_sampler_params(fit.betabinom, inc_warmup = TRUE)
#lapply(sampler_params, summary, digits = 2)
##  A dataframe with a posteriori distributions
fit.betabinom.df <- as.data.frame(fit.betabinomb, c("alpha","beta"))

## Derived quantities: Expected species richness

################################################################################
## Fit truncated beta binomial distribution in Bayesian model ##
## Beta-binomial distribution truncated at zero
## painfully slow!
################################################################################

## List with all data needed for the model
m1.data <- list(Nsites=1900, Sobs = nrow(dados), y = dados$N.plots)
## Run stan (much faster than the hierachical version above)
fit.betabinom.2 <- stan(file="betabinom2.stan", data=m1.data, iter=200, chains=3, cores=3)
fit.betabinom.2b <- stan(file="betabinom2.stan", data=m1.data, iter=200, chains=3, cores=3)
## Diagnostics
traceplot(fit.betabinom.2b, c("alpha", "beta", "phi", "lambda", "Sest")) 
print(fit.betabinom.2, c("alpha","beta")) ## 

