library(VGAM)

#' Tovo et al estimate of total species richness from negative binomial
#' @param fit An object of class fitsad with the fit of a truncated negative binomial to the data
#' @param cf  vector of two elements with the coefficients of fit of a Negative Binomial do the data by fitnbinom (size, and mu)
#' @param S.obs integer, observed species richness
#' @param p real positive, the proportion of the community that has been sampled
#' @param CI logic, calculate Confidence intervals based on teh Ci of the estimated parameters?
tovo <- function(fit, cf, S.obs, p, CI=FALSE){
    if(missing(S.obs))
        S.obs <- length(fit@data$x)
    if(missing(cf))
        cf <- unname(coef(fit))
    csi <- tovo.csi(cf[1], cf[2], p=p, log=FALSE)
    ## Estimated number of species 
    S.est <- S.obs*(1-(1-csi)^cf[1]) / (1-(1-csi.p)^cf[1])
    if(CI){
        ci <- confint(fit)
        p.low <- 1 - tovo.csi(ci[1,1], ci[2,1], p=p, log=FALSE)
        p.up <- 1 - tovo.csi(ci[1,2], ci[2,2], p=p, log=FALSE)
        S.low <- tovo(cf=ci[,1], S.obs=S.obs, p=p)
        S.upp <- tovo(cf=ci[,2],S.obs=S.obs, p=p)
        CIs <- rbind(ci, c(p.low,p.up), c(S.low,S.upp)  )
        rownames(CIs)[c(3,4)] <- c("prob", "S.est")
        cat("Estimated species richness:", S.est, "\n",
            "95% CI:",CIs[3,2],"-",CIs[3,1], "\n")
        return(invisible(list(S.est=S.est, CIs=CIs)))
    }
    else
        return(S.est)
}

#' Calculate re-scaled value of mu from Tovo et al re-scaled nbinom
tovo.mu <- function(size, mu, p, log=TRUE){
    csi.p <- mu/(mu+size)
    if(log){
        ##re-scaled csi (log)
        lcsi <- tovo.csi(size, mu, p, log=TRUE)
        ##re-scaled mu
        f1 <- function(m) l.csi-(log(m)-log(m+size))
    }
    else{
        csi <- tovo.csi(size, mu, p, log=FALSE)
        f1 <- function(m) csi - (m/(m+size))
    }
    uniroot(f1, c(exp(-100), exp(100)))
}

#'Tovo et al csi calulation
tovo.csi <- function(size, mu, p, log=FALSE){
    csi.p <- mu/(mu+size)
    if(log)
        y <- log(csi.p)-log( (p+(1-p)*csi.p) )
    else
        y <- csi.p/(p1+(1-p1)*csi.p)
    return(y)
    }

#' Zero-truncated negative binomial RAD
rad.posnegbin <- function(S.est, size, prob, npoints=100){
    pp <- rev(ppoints(S.est))
    x <- seq(1,S.est, length=npoints)
    y <- qposnegbin(p = pp[x], size = size, prob = prob)
    data.frame(x, y)
    }

#' Find total number of species sampling probabilities of occurence from a beta distribution
#' @param shape1 parameter shape1 of a beta distribution
#' @param shape2 parameter shape2 of a beta distribution
#' @param S.obs observed number of species
#' @param size number of replicates (that is, sites where the species were surveyed)
#' @param upper upper bound of the interval of possible richness values for optimise
find.S <- function(shape1, shape2, S.obs, size, upper){
    f1 <- function(S) abs(sum( 1- (1-rbeta(S, shape1, shape2))^size ) - S.obs)
    optimise(f1, lower = S.obs, upper = S.obs*upper)
    }

#' Species richness estimate from a zero-truncated beta-binomial distribution
#' @param mu parameter mu of the compounding beta distribution (mu = shape1 / (shape1 + shape2))
#' @param rho parameter rho of the compounding beta distribution (rho = 1 / (1 + shape1 + shape2))
#' @param size number of trials of the compounding binomial distribution
#' @param S.obs observed number of species
#' @param shape1  alternative parametrization of the beta distribution, ignored if mu is provided
#' @param shape2  alternative parametrization of the beta distribution, ignored if rho is provided
bb.Sest <- function(mu, rho, size, Sobs, shape1, shape2){
    if(missing(mu))
        mu <- shape1 / (shape1 + shape2)
    if(missing(rho))
        rho <- 1 / (1 + shape1 + shape2)
    Sobs / (1 - VGAM::dbetabinom(0, size=size, prob = mu, rho = rho))
    }

#'qq-plot to compare to a parametric distribution
QQ.plot <- function(x, distr, coef, trunc=NA, plot=TRUE, line=TRUE, ...){
    x.sorted <- sort(x)
    S <- length(x)
    p <- ppoints(S)
    if(!is.na(trunc))
        q <- do.call(sads::qtrunc, list(distr, p = p, trunc = trunc, coef=coef))
    else{
        qdistr <- get(paste("q", distr, sep=""), mode = "function")
        q <- do.call(qdistr, c(list(p = p), coef))
    }
    if(plot){
        dots <- list(...)
        if(!"main" %in% names(dots)) dots$main = "Q-Q plot"
        if(!"xlab" %in% names(dots)) dots$xlab = "Theoretical Quantile"
        if(!"ylab" %in% names(dots)) dots$ylab = "Sample Quantiles"
        do.call(graphics::plot, c(list(x=q, y=x.sorted),dots))
        if(line) abline(0, 1, col = "red", lty = 2)
    }
    return(invisible(data.frame(theoret.q=q, sample.q=x.sorted)))
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
## rando deviates
rbetabinom.t0 <- function(n, size, mu, rho) rtrunc("betabinom", n, trunc=0, coef=list(size=size, mu=mu, rho=rho))
## Poilog truncated at zero
dpoilog.t0 <- function(x, mu, sig, log=FALSE) dtrunc("poilog", x, trunc=0, coef=list(mu=mu, sig=sig), log=log)
ppoilog.t0 <- function(q, mu, sig, ...) ptrunc("poilog", q, trunc=0, coef=list(mu=mu, sig=sig), ...)
qpoilog.t0 <- function(p, mu, sig, ...) qtrunc("poilog", p, trunc=0, coef=list(mu=mu, sig=sig), ...)
rpoilog.t0 <- function(n, mu, sig) rtrunc("poilog", n, trunc=0, coef=list(mu=mu, sig=sig))

#' testing: estimated species richnedd from zero-truncate negative binomial lognormal
nbln.Sest <- function(mu, sig, k, Sobs, a=1){
    f1 <- function(x){
        dlnorm(x, mu, sig)*dnbinom(0, size=k, mu=x*a)
    }
    integrate(f1, 0, Inf)
}

#' Fisher logseries predicted species abundance
ls.pred <- function(rank, N, alpha){
    X <- N/(alpha+N)
    alpha*(X^rank)/rank
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

est.kv <- Vectorize(est.k, c("mu", "nzeroes"))

#' analytical
est.ka <- function(mu, nzeroes, Nplots){
    (mu * nzeroes)/(Nplots - nzeroes)
    }
 
#'Drawn N samples from Negative binomial for a pair of parameters of the NB (mu and size) and the sums up these values
rnbinom2 <- function(mu, size, N){
        y <- rnbinom(n = N, mu  = mu, size = size)
        sum(y)
        }

#'Drawn N samples from a Poisson distribution and then sums up these values
rpois2 <- function(lambda, N){
        y <- rpois(n = N, lambda  = lambda)
        sum(y)
    }

#' Simulates the number of non-occupied plots for a Negative Binomial sampling or a Poisson samplig
sim.occ <- function(mu, size, N, pois.samp=TRUE){
    if(pois.samp)
        lp <- -mu
    else
        lp <- size*(log(size)- log(mu+size))
    p0 <- exp(N*lp)
    sample(0:1, size=1, prob=c(p0, 1-p0))
}
