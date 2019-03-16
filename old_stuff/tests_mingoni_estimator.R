library(MASS)
library(VGAM)
## Viability of using Mingoni 1999 (J. Appl Stat) Bayesian estimator

## 1. with an arbitrary prior for the true number of species
## Assuming that prior for occupancies is a bet fitted to the observed occupancies
## Fitting a beta distribution to frequencies
beta.f <- fitdistr(dados$N.plots/N.plots, "beta", start=list(shape1=1, shape2=1))
(b.cf <- coef(beta.f))

## gamma0 is eq 1 , which is an integral involving the beta distribution above
g0 <- function(shape1, shape2, n){
    f1 <- function(p) ((1-p)^n) * dbeta(p, shape1, shape2)
    integrate(f1, 0,1)
}

## Using estimates of beta distribution fitted to observed occupancies (see Rnw file)
(g0.1 <- g0(cf.beta[1], cf.beta[2], n = N.plots))
## Beta truncated at 0.1/Nplots
(g0.2 <- g0(cf.betat.1[1], cf.betat.1[2], n = N.plots))

## The common term of eq1 includes a binomial coeeficiente which is computationally doable only in log scale
## Checking if we can calculate it in log scale and then get it back to arithmetic scale

m.eq1.a <- function(S,s,g0, pi.prior){
    lchoose(s+S,s) + S*log(g0) + pi.prior(S+s, log=TRUE)
}

## Two priors for the true number of species
## Uniform between observed number and twice
pi.p1 <- function(x,...) dunif(x, nrow(dados), 5*nrow(dados), ...)
## Normal with mean at midpoint between observed and twice
pi.p2 <- function(x, ...) dnorm(x, mean = 3*nrow(dados), sd = nrow(dados), ...)


m.eq1.a(1, nrow(dados), g0.1$value, pi.p1)
m.eq1.a(nrow(dados)*2, nrow(dados), g0.1$value, pi.p1)
m.eq1.a(nrow(dados)*2, nrow(dados), 1, pi.p1) ## too large
exp(m.eq1.a(0:nrow(dados)*2, nrow(dados), g0.1$value, pi.p1))

## Seems to work
## Equation 1

S.pi <- function(S,s,g0, pi.prior) {
    s +
        sum ((0:S)  * exp( m.eq1.a(0:S, nrow(dados), g0, pi.prior)) ) /
        sum ( exp( m.eq1.a(0:S, nrow(dados), g0, pi.prior) ) )
    }

## Test
S.pi(S = nrow(dados)*2, s  = nrow(dados), g0 = g0.1$value, pi.prior = pi.p1) ## Low
S.pi(S = nrow(dados)*2, s  = nrow(dados), g0 = g0.1$value, pi.prior = pi.p2) ## As above, so not sensible to S prior
S.pi(S = nrow(dados)*2, s  = nrow(dados), g0 = g0.2$value, pi.prior = pi.p1) ## NaN, seems to go to Inf for g0>0.14


## 2. Checking the solution for a zero-truncated negative binomial prior for S

S.pi.5 <- function(s, g0, R, q){
    (s + R*q*g0)/(1-q*g0)
}

## Ad hoc procedure to define R and q (section 4)

S.pi.5(s = nrow(dados), g0 = g0.1$value,
       R = ifelse(floor(N.plots/nrow(dados))==0, 1, floor(N.plots/nrow(dados))),
       q = mean(dados$N.ind==1)
       )

S.pi.5(s = nrow(dados), g0 = g0.1$value,
       R = N.plots/nrow(dados),
       q = mean(dados$N.ind==1)
       )

S.pi.5(s = nrow(dados), g0 = g0.2$value,
       R = ifelse(floor(N.plots/nrow(dados))==0, 1, floor(N.plots/nrow(dados))),
       q = mean(dados$N.ind==1)
       )
## Maximum value for a given R and q, making g0 -> 1 (section 3) ##
S.pi.5(s = nrow(dados), g0 = 1,
       R = ifelse(floor(N.plots/nrow(dados))==0, 1, floor(N.plots/nrow(dados))),
       q = mean(dados$N.ind==1)
       )
## Checking with eq6: ok
(nrow(dados) + 
 ifelse(floor(N.plots/nrow(dados))==0, 1, floor(N.plots/nrow(dados)))* mean(dados$N.ind==1)) /
    (1 - mean(dados$N.ind==1))



## Checking with a larger q: some increase but not impressive
S.pi.5(s = nrow(dados), g0 = g0.1$value,
       R = ifelse(floor(N.plots/nrow(dados))==0, 1, floor(N.plots/nrow(dados))),
       q = 0.9
       )

