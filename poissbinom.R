## Some tries with Poisson-binomial
library(fitdistrplus)
library(sads)
library(poisbinom)
library(lme4)
##
dados <- read.csv2("../data.csv", as.is=TRUE)
S.obs <- nrow(dados)
## Probability of being undetected (log scale)
lp0 <- 1900*log(1-dados$N.plots/1900)
descdist(lp0, boot=1000)
plotdist(lp0)
p0.beta <- fitdist(lp0, "beta")
## Distribution of probabilities of being detected in the 1900 samples
pd <- 1-(1-dados$N.plots/1900)^1900
hist(pd)
## Maybe a truncated beta
plotdist(pd)
descdist(pd, boot=1000)
pd.beta <- fitdist(pd, "beta")
par(mfrow=c(2,2))
plot(pd.beta)
cf1 <- coef(pd.beta)
## A truncated negative binomial to the number of plots
hist(dados$N.plots)
fit2 <- fitnbinom(dados$N.plots, trunc=0, start.value=list(mu=mean(dados$N.plots),size=100))
par(mfrow=c(2,2))
plot(fit2)
par(mfrow=c(1,1))
## Bom fit!
fit2.cf <- coef(fit2)
## Ajuste a uma gama truncada
## Selecionando
fit3.0 <- fitgamma(dados$N.plots)
fit3.1 <- fitgamma(dados$N.plots, trunc=0.1, start.value=c(shape=0.5,scale=100))
fit3.3 <- fitgamma(dados$N.plots, trunc=0.3, start.value=c(shape=0.5,scale=100))
fit3.5 <- fitgamma(dados$N.plots, trunc=0.5, start.value=c(shape=0.5,scale=100))
AICtab(fit3.0, fit3.1, fit3.3, fit3.5)
par(mfrow=c(2,2))
plot(fit3.5)
par(mfrow=c(1,1))
cf3 <- coef(fit3.5)

## Uma simulacao tomando sorteios desta gamma para usar como vetores de probabilidades em uma Poisson-binomial
## Prob de nao ser detectada por parcela
p0 <- 1-rgamma(1.3e4, shape=cf3[1], rate=cf3[2])/1900
## Prob de ser detectada
p3 <- 1 - p0^1900
## Predicted richness
sum(p3)
## Distribuicao dos valores de S
plot(0:length(p3), dpoisbinom(0:length(p3), p3), type="l")
abline(v=S.obs)

### Ajuste a beta da probabilidade de deteccao
f.st <- fitdist(dados$N.plots/1900, "beta")
cf.st <- coef(f.st)
## Fit by mle
## Likelihood
LL.beta <- function(s1, s2) -sum(dbeta(dados$N.plots/1900, shape1=s1, shape2=s2, log=TRUE))
## Fit
f.beta <- mle2(LL.beta, start=list(s1=cf.st[1], s2=cf.st[2]))
cf.beta <- coef(f.beta)
##Sample of 1e4 values
beta.p <- 1-(1-rbeta(1e4, cf.beta[1], cf.beta[2]))^1900
## Expected richness
sum(beta.p)
## Truncated beta
## Verossimilhanca
LL.betat <- function(s1, s2){
    -sum ( dtrunc("beta", x = dados$N.plots/1900, trunc=0.5/1900, coef = list(shape1=s1, shape2=s2), log=TRUE) )
    }
f.betat <- mle2(LL.betat, start=list(s1=cf.st[1], s2=cf.st[2]))
cf.betat <- coef(f.betat)
## Trying with a tight truncation point: convergence problems
LL.betat2 <- function(s1, s2){
    -sum ( dtrunc("beta", x = dados$N.plots/1900, trunc=0.7/1900, coef = list(shape1=s1, shape2=s2), log=TRUE) )
    }
f.betat2 <- mle2(LL.betat2, start=list(s1=cf.betat[1], s2=cf.betat[2]))
## Checking all fits
AICtab(f.beta, f.betat)
##Sample of 
betat.p <- 1-(1-rbeta(1.55e4, cf.betat[1], cf.betat[2]))^1900
## Expected richness
sum(betat.p)

## A function to estimate S
estS.pb <- function(pp, S.obs, nsites, upper=10){
    if(missing(S.obs))
        S.obs <- length(pp)
    cf.st <- coef(fitdist(pp, "beta"))
    LL.betat <- function(s1, s2){
        -sum ( dtrunc("beta", x = pp, trunc=0.5/nsites, coef = list(shape1=s1, shape2=s2), log=TRUE) )
    }
    fit <- mle2(LL.betat, start = list(s1 = cf.st[1], s2 = cf.st[2]))
    cf <- coef(fit)
    ## LL.pb <- function(S){
    ##     ps <- 1- (1-rbeta(S, cf[1], cf[2]))^nsites
    ##     -sum(dpoisbinom(S.obs, ps, log=TRUE))
    ## }
    ## mle2(LL.pb, start=list(S=S.obs), method="Brent", lower=S.obs, upper=S.obs*upper)
    ## Include here sampling of beta parameters from a bivariate gaussian
    f1 <- function(S) abs(sum(1- (1-rbeta(S, cf[1], cf[2]))^nsites) - S.obs)
    optimise(f1, lower = S.obs, upper = S.obs*upper)
}

## Test
S.est <- estS.pb(pp = dados$N.plots/1900, nsites=1900)
summary(S.est)
## Checking loglikelihhod profile
## Likelihood function
nsites <- 1900
LL.pb <- function(S){
         ps <- 1- (1-rbeta(S, cf.betat[1], cf.betat[2]))^nsites
         -sum(dpoisbinom(S.obs, ps, log=TRUE))
}
S.seq <- seq(15200, 15800, by=100)
S.seq <- sort(c(S.seq,S.est$minimum))
nrep <- 10
S.prf <- matrix( ncol=length(S.seq), nrow=nrep)
for(i in 1:nrep){
    S.prf[i,] <- sapply(S.seq, LL.pb)
}
S.prfm <- apply(S.prf,2,mean)
plot(S.seq, S.prfm-min(S.prfm))
abline(h=2, lty=2)

## A parametric bootstrap: upper bounded by the estimated richness above
## Matrix with simulated proportion of sites occupied by each species,
## sampled from the observed occupancies
##pps <- dados$N.plots/1900
results <- c()
results[1] <- S.est$minimum[[1]]
nrep <- 100
nsites <- 1900
for(i in 2:nrep){
    ## Draw occupancy probabilities  from the fitted beta distribution for each species of the metacommunity (estimated above)
    pps <- rbeta(round(results[i-1]), cf.betat[1], cf.betat[2])
    ## Generate new values for the number of sites occupied by each species in the metacommunity
    occup <- rbinom(length(pps), size=nsites, prob=pps)
    ## Estimate the value of S.esp from the simulated vector of observed occupancies
    results[i] <- estS.pb(occup[occup>0]/nsites, S.obs=S.obs, nsites=nsites, upper=5)$minimum
    }
## Bootstrap confidence interval
summary(results)
quantile(results, c(0.025, 0.975))
              
### Ajuste dos logitos a uma normal com um glmm: Não é uma boa
m1 <- glmer(cbind(N.plots, 1900-N.plots) ~ 1 + (1|Species), data=dados, family=binomial)
## Estimated mean and standard deviation of the gaussian of the detectabilities
m1.sd <- attr(VarCorr(m1)$Species, "stddev")[[1]]
m1.mean <- fixef(m1)[[1]]
## Sample of 1e4 values from this gaussian
m1.lp <- rnorm(1e4, m1.mean, m1.sd)
## Convert to probabilities
m1.p <- exp(m1.lp)/(1+exp(m1.lp))
## Probabilities of occurr in at least one plot
m1.pn <- 1 - (1-m1.p)^1900
## Expected species richness: too high
sum(m1.pn)
