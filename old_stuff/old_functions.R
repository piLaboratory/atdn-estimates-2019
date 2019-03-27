
#' Find total number of species by sampling probabilities of occurence
#' from a beta distribution
#' 
#' @param shape1 parameter shape1 of a beta distribution.
#' @param shape2 parameter shape2 of a beta distribution.
#' @param S.obs observed number of species.
#' @param size number of replicates (that is, sites where the species
#'     were surveyed). 
#' @param upper upper bound of the interval of possible richness
#'     values for optimise.
#' 
find.S <- function(shape1, shape2, S.obs, size, upper){
    f1 <- function(S) abs(sum( 1- (1-rbeta(S, shape1, shape2))^size ) - S.obs)
    optimise(f1, lower = S.obs, upper = S.obs*upper)
    }


#' Truncated beta-binomial
## pdf
dbetabinom.t0 <- function(x, size, mu, rho, log=FALSE){
    y <- VGAM::dbetabinom(x, size, mu, rho, log=TRUE) - (1-VGAM::dbetabinom(0, size, mu, rho, log=TRUE))
    if(log)
        return(y)
    else
        return(exp(y))
}
## cdf
pbetabinom.t0 <- function(q, size, mu, rho, log=FALSE){
    y <- VGAM::pbetabinom(x, size, mu, rho, log=TRUE) - (1-VGAM::dbetabinom(0, size, mu, rho, log=TRUE))
    if(log)
        return(y)
    else
        return(exp(y))
}
## quantile
qbetabinom.t0 <- function(p, size, mu, rho, ...) qtrunc("betabinom", p, trunc=0, coef=list(size=size, mu=mu, rho=rho), ...)
## random deviates
rbetabinom.t0 <- function(n, size, mu, rho) rtrunc("betabinom", n, trunc=0, coef=list(size=size, mu=mu, rho=rho))
## Poilog truncated at zero
dpoilog.t0 <- function(x, mu, sig, log=FALSE) dtrunc("poilog", x, trunc=0, coef=list(mu=mu, sig=sig), log=log)
ppoilog.t0 <- function(q, mu, sig, ...) ptrunc("poilog", q, trunc=0, coef=list(mu=mu, sig=sig), ...)
qpoilog.t0 <- function(p, mu, sig, ...) qtrunc("poilog", p, trunc=0, coef=list(mu=mu, sig=sig), ...)
rpoilog.t0 <- function(n, mu, sig) rtrunc("poilog", n, trunc=0, coef=list(mu=mu, sig=sig))

#' testing: estimated species richness from zero-truncated negative binomial lognormal
nbln.Sest <- function(mu, sig, k, Sobs, a=1){
    f1 <- function(x){
        dlnorm(x, mu, sig)*dnbinom(0, size=k, mu=x*a)
    }
    integrate(f1, 0, Inf)
}

#' estimates parameter k from NB using # of observed zeroes
est.k <- function(mu, nzeroes, Nplots){
    f1 <- function(k) {
        p <- dnbinom(0, size=k, mu=mu)
        -dbinom(nzeroes, size=Nplots, prob=p, log=TRUE)
    }
    fit <- mle2(f1, method="Brent", start=list(k=1), lower=1e-9, upper=10)
    unname(coef(fit))
}
## Vectorized version
est.kv <- Vectorize(est.k, c("mu", "nzeroes"))

#'Draws N samples from Negative binomial for a pair of parameters of the NB (mu and size) and the sums up these values
rnbinom2 <- function(mu, size, N){
        y <- rnbinom(n = N, mu  = mu, size = size)
        sum(y)
        }

#'Draws N samples from a Poisson distribution and then sums up these values
rpois2 <- function(lambda, N){
        y <- rpois(n = N, lambda  = lambda)
        sum(y)
    }
