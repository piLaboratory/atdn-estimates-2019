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
    csi.p <- cf[2]/(sum(cf))
    csi <- csi.p/(p+(1-p)*csi.p)
    ## Estimated number of species 
    S.est <- S.obs*(1-(1-csi)^cf[1]) / (1-(1-csi.p)^cf[1])
    if(CI){
        ci <- confint(fit)
        low <- tovo(cf=ci[,1], S.obs=S.obs, p=p)
        upp <- tovo(cf=ci[,2],S.obs=S.obs, p=p)
        CIs <- rbind(ci,c(low,upp))
        rownames(CIs)[3] <- "S est"
        cat("Estimated species richness:", S.est, "\n",
            "95% CI:",CIs[3,2],"-",CIs[3,1], "\n")
        return(invisible(list(S.est=S.est, CIs=CIs)))
    }
    else
        return(S.est)
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
