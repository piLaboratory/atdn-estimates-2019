library(VGAM)
library(untb)
library(sads)
library(abc)
library(parallel)

#' Fit to negative binomial with log-link
fitnbinom2 <- function (x, trunc = 0, start.value, ...) {
    dots <- list(...)
    ##if ((any(x <= 0) & !is.null(trunc)) | any(!is.wholenumber(x))) 
    ##    stop("All x must be positive integers")
    if (!is.null(trunc)) {
        if (min(x) <= trunc) 
            stop("truncation point should be lower than the lowest data value")
    }
    if (missing(start.value)) {
        muhat <- length(x)/(length(x) + mean(x))
        sizehat <- muhat * mean(x)
    }
    else {
        sizehat <- start.value[[1]]
        muhat <- start.value[[2]]
    }
    if (is.null(trunc)) {
        LL <- function(lsize, lmu) -sum(dnbinom(x, size = exp(lsize), 
            mu = exp(lmu), log = TRUE))
    }
    else {
        LL <- function(lsize, lmu) -sum(dtrunc("nbinom", x = x, 
            coef = list(size = exp(lsize), mu = exp(lmu)), trunc = trunc, 
            log = TRUE))
    }
    result <- do.call("mle2", c(list(LL, start = list(lsize = log(sizehat), 
        lmu = log(muhat)), data = list(x = x)), dots))
    new("fitsad", result, sad = "nbinom", #distr = distr.depr, 
        trunc = ifelse(is.null(trunc), NaN, trunc))
}

#' Tovo et al estimate of total species richness from negative
#'     binomial
#' 
#' @param fit An object of class fitsad with the fit of a truncated
#'     negative binomial to the data. 
#' @param cf  vector of two elements with the coefficients of fit of a
#'     Negative Binomial to the data by fitnbinom (size, and mu). 
#' @param S.obs integer, observed species richness. 
#' @param p real positive, the proportion of the community that has
#'     been sampled. 
#' @param CI logic, calculate Confidence intervals based on the 95% CI of
#'     the estimated parameters? 
tovo <- function(fit, cf, S.obs, p, CI=FALSE, loglink=FALSE){
    if(missing(S.obs))
        S.obs <- length(fit@data$x)
    if(missing(cf)){
        cf <- unname(coef(fit))
        if(loglink)
            cf <- exp(cf)
        }
    csi <- tovo.csi(cf[1], cf[2], p = p, log=FALSE)
    csi.p <- cf[2] / (cf[1] + cf[2])
    ## Estimated number of species 
    S.est <- S.obs*(1-(1-csi)^cf[1]) / (1-(1-csi.p)^cf[1])
    if(CI){
        if(loglink)
            ci <- exp(confint(fit))
        else
            ci <- confint(fit) 
        p.low <- 1 - tovo.csi(ci[1,1], ci[2,1], p=p, log=FALSE)
        p.up <- 1 - tovo.csi(ci[1,2], ci[2,2], p=p, log=FALSE)
        S.low <- tovo(cf=ci[,1], S.obs=S.obs, p=p)
        S.upp <- tovo(cf=ci[,2],S.obs=S.obs, p=p)
        CIs <- rbind(ci, c(p.low,p.up), c(S.low,S.upp)  )
        rownames(CIs)[c(3,4)] <- c("prob", "S.est")
        cat("Estimated species richness:", S.est, "\n",
            "95% CI:",CIs[4,2],"-",CIs[4,1], "\n")
        return(invisible(list(S.est=S.est, CIs=CIs)))
    }
    else
        return(S.est)
}

#' Calculates re-scaled value of mu from Tovo et al re-scaled nbinom
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

#'Tovo et al csi calculation
tovo.csi <- function(size, mu, p, log=FALSE){
    csi.p <- mu/(mu + size)
    if(log)
        y <- log(csi.p)-log( (p + (1-p)*csi.p) )
    else
        y <- csi.p/(p + (1-p)*csi.p)
    return(y)
    }

#' Find the value of Csi in Tovo's TNB given observed and total
#'     species richness
#' 
#' @param S total species richness.
#' @param S.obs species richness in the sample.
#' @param k estimated size parameter of the Truncated negative
#'     binomial. 
#' @param csi.p Csi parameter (1 - prob, see Tovo et al) estimated for
#'     the sample. 
tovo.Scsi <- function(S, S.obs, k, csi.p, prob = FALSE){
    C1 <- S.obs / (1-(1-csi.p)^k)
    if(!prob)
        return( -((C1 - S)/C1)^(1/k) + 1)
    else
      ((C1 - S)/C1)^(1/k)  
    }

#' utility function: incomplete beta function
ibeta <- function(x,a,b, log=FALSE){
    y <- pbeta(x,a,b, log.p=TRUE) + lbeta(a,b)
    if(log)
        return(y)
    else
        exp(y)
}

#' CDF of logseries, using incomplete beta function (https://en.wikipedia.org/wiki/Logarithmic_distribution)
#' 
pls2 <-function(x, alpha, N){
    p <- N/(N+alpha)
    1 + ibeta(p, x+1, 1e-12)/log(1-p) 
    }

#' Continuous approximation for quantile function for Log-series distribution
qls2 <- function(p, N, alpha, lower=3e-9, upper=3e9){
    f2 <- function(target){
        f1 <- function(x) pls2(x,alpha, N) - target
        uniroot(f1, lower=lower, upper=upper)$root
    }
    sapply(p, f2)
}

#' Log-series RAD
#' @description Generates a given number of points of the RAD of a
#'     LS, given the total number of species and the parameters of
#'     the distribution.
rad.ls <- function(S, N, alpha, npoints = round(S), ...){
    S.r <- round(S)
    if(missing(alpha))
        alpha <- fishers.alpha(N = N, S = S)
    pp <- rev(ppoints(S.r))
    x <- seq(1,S.r, length=npoints)
    y <- qls2(p = pp[x], N = N, alpha = alpha, ...)
    data.frame(x, y)
    }

#' Continuous approximation for quantile function for TNB distribution
qposnegbin2 <- function(p, size, prob, lower=3e-9, upper=3e9){
    f2 <- function(target){
        f1 <- function(x) pposnegbin(x, size, prob) - target
        uniroot(f1, lower=lower, upper=upper)$root
    }
    sapply(p, f2)
}


#' Zero-truncated negative binomial RAD
#' @details Generates a given number of points of the RAD of a
#'     ZTNB, given the total number of species and the parameters of
#'     the distribution.
#' 
#' @param S total species richness in the RAD.
#' @param size size parameter of the Zero-truncated negative
#'     binomial. 
#' @param prob prob parameter of the Zero-truncated negative binomial.
#' @param npoints number of points along teh RAD to retunr (defaults
#'     to the total number of species, which implies that the expected
#'     abundance will be calculated for each species (can be slow for
#'     a large number of species).
rad.posnegbin <- function(S, size, prob, npoints = round(S),...){
    S.r <- round(S)
    pp <- rev(ppoints(S.r))
    x <- seq(1,S.r, length=npoints)
    y <- qposnegbin2(p = pp[x], size = size, prob = prob, ...)
    data.frame(x, y)
    }


#' Species richness estimate from a zero-truncated beta-binomial
#' distribution.
#' 
#' @param mu parameter mu of the compounding beta distribution (mu =
#'     shape1 / (shape1 + shape2))
#' @param rho parameter rho of the compounding beta distribution (rho
#'     = 1 / (1 + shape1 + shape2))
#' @param size number of trials of the compounding binomial
#'     distribution
#' @param S.obs observed number of species
#' @param shape1 alternative parametrization of the beta distribution,
#'     ignored if mu is provided
#' @param shape2 alternative parametrization of the beta distribution,
#'     ignored if 'rho' is provided
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


#' Fisher logseries predicted species abundance
ls.pred <- function(rank, N, alpha){
    X <- N/(alpha+N)
    alpha*(X^rank)/rank
}

#'Predicted species richness from upsacling a sampled SAD
ls.estS <- function(rad, N){
    y <- sort(rad[rad>0])
    ls.fit <- fitls(y)
    alpha <- coef(ls.fit)[[2]]
    alpha*log(1 + N/alpha)
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

#' Simulates species occupancies in a sample of a RAD
#'
#' @details Given the density of a set of species per plot in a
#'     community, this functions calculates the probability of
#'     recording each species per plot and then simulates the number
#'     of plots each species is recorded in a sample of N plots. The
#'     sample is simulated assuming random distribution of
#'     conspecifics (Poisson sample) or cumpled distribution (negative
#'     binomial sample).
#' 
#' @param mu vector of positive reals, species densities (individuals
#'     per plot unit) in the RAD to be sampled.
#' @param size vector of positive reals, value of parameter 'size'
#'     of the negative binomial for each species in the RAD (see
#'     'dnbinom'). 
#' @param N positive integer, number of plots to be sampled.
#' @param pois.samp logical, if TRUE simulates a Poisson sample;
#'     simulates a negative binomial sample with parameters given by
#'     'size' otherwise. 
#'
#' @return a vector of same length of 'mu' and 'size' with the number
#'     of plots in the simulated sample each species was recorded.
sim.occ <- function(mu, size, N, pois.samp=TRUE){
    if(!pois.samp){
        if(length(mu)!=length(size))
            stop(" 'mu' and 'size' must be of the same length")
        else
            lp <- size*(log(size) - log(mu+size))
    }
    else
        lp <- -mu
    rbinom(n = length(mu), size = N , prob = 1-exp(lp) )
}


#' Simulates persence/absence of species in a sample of a RAD
#'
#' @details Utility function used by 'Pois.samp' and 'NB.samp'
#'     functions, see below.  Given the density of species per plot in
#'     a community, this functions calculates the probability of
#'     recording the species in a sample of N plots and then simulates
#'     an event of recording/not recording sampling from a Bernoulli
#'     distribution. The sample is simulated assuming random
#'     distribution of conspecifics (Poisson sample) or cumpled
#'     distribution (negative binomial sample).
#' 
#' @param mu positive real, species densities (individuals
#'     per plot unit) in the RAD to be sampled.
#' @param size positive real, value of parameter 'size'
#'     of the negative binomial for each species in the RAD (see
#'     'dnbinom'). 
#' @param N positive integer, number of plots to be sampled.
#' @param pois.samp logical, if TRUE simulates a Poisson sample;
#'     simulates a negative binomial sample with parameters given by
#'     'size' otherwise. 
#'
#' @return a value of zero or one, indicating presence/absence the
#'     species in the simulated sample.
sim.pres <- function(mu, size, N, pois.samp=TRUE){
    if(!pois.samp){
        if(length(mu)!=1|length(size)!=1)
            stop("This function is not vectorized, please provide a single value of 'mu' and 'size' ")
        }
    else if(length(mu)!=1)
        stop("This function is not vectorized, please provide a single value of 'mu' ")
    if(pois.samp)
        lp <- -mu
    else
        lp <- size*(log(size) - log(mu+size))
    p0 <- exp(N*lp)
    sample(0:1, size=1, prob=c(p0, 1-p0))
}

#' Simulates samples of population sizes using a Poisson sample
#' 
#' @details This function performs a simulation of which species
#'     abundaces would be included in a Poisson sample of a regional
#'     RAD. The population sizes of the included species are then set
#'     to the population sizes of the 1st, 2nd .. Nth most abundant
#'     species recorded. Two simulated RADs are simulated: (1) a "no estimate
#'     error" vector of abundances, under the assumption that the only
#'     source of uncertainty is which species will be included in the
#'     sample (that is, which species would be detected); and (2) a
#'     "estimate error" vector of abundances, that include the
#'     uncertatinty in the estimation of population sizes (defined by
#'     'lmean.sd' and 'lsd.sd') If nrep > 1 then the simulation is
#'     repeated nrep times and the abundance of the 1st, 2nd ... Nth
#'     species is taken from the mean abundances at each species rank
#'     over repetitions.
#' 
#' @param rad vector of positiev reals, species population sizes in
#'     the RAD to be sampled.
#' @param tot.area positive real, total area of the community to be
#'     sampled. The area unit is one plot.
#' @param n.plots positive integer, number of sampling units
#'     (e.g. plots) to be drawn out of the total number of plots.
#' @param nrep positive integer, number of repetitions of the
#'     simulated sampling
#'
#' @return A dataframe with the vectors of the expected abundances of
#'     the 1nd, 2nd, ... Nth species without and with estimation
#'     errors of population sizes.
Pois.samp <- function(rad, tot.area, n.plots,
                      lmean.sd, lsd.sd, nrep = 1){
    index <- order(rad, decreasing=TRUE)
    rad <- rad[index]
    lmean.sd <- lmean.sd[index]
    m1 <- m2 <- matrix(0,nrow=length(rad), ncol=nrep)
    for(j in 1:nrep){
        y1 <- mapply(sim.pres, mu = rad/tot.area,
                     MoreArgs=list(N = n.plots))
        ## Select species included in the sample
        m1[1:sum(y1),j] <- rad[y1>0]
        ## Sample standard deviations for abundances estimates
        sd1 <- exp( rnorm ( sum(y1), lmean.sd[y1>0], sd = lsd.sd) )
        ## Builds simulated rad with simulated abundance estimates
        ## with Gaussian estimation error. 
        m2[1:sum(y1),j] <- rnorm(sum(y1), mean = rad[y1>0], sd1)
    }
    data.frame( no.est.error = apply(m1,1,mean), with.est.error = apply(m2,1,mean))
}

#' Simulates a Negative Binomial samples of population sizes from a RAD
#'
#' @details This function performs a simulation of which species
#'     abundaces would be included in a Negative Binomial sample of a
#'     regional RAD. The population sizes of the included species are
#'     then set to the population sizes of the 1st, 2nd .. Nth most
#'     abundant species recorded. The expected dispersion parameter of
#'     negative binomial sampling is allowed to vary across species,
#'     and is assumed to have a lognormal error which is added in the
#'     simulations. This is a simulation of how the distribution of
#'     total population sizes would look like, assuming that the there
#'     is a method to estimate the population sizes of the species
#'     recorded. Two simulated RADs are simulated: (1) a "no estimate
#'     error" vector of abundances, under the assumption that the only
#'     source of uncertainty is which species will be included in the
#'     sample (that is, which species would be detected); and (2) a
#'     "estimate error" vector of abundances, that include the
#'     uncertatinty in the estimation of population sizes (defined by
#'     'lmean.sd and 'lsd.sd') If nrep > 1 then the simulation is
#'     repeated nrep times and the abundance of the 1st, 2nd ... Nth
#'     species is taken from the mean abundances at each species rank
#'     over repetitions.
#' 
#' @param rad a vector with the species population sizes in the RAD to be sampled
#' @param tot.area total area of the community to be sampled. The area unit is one plot
#' @param n.plots number of sampling units (e.g. plots) to be drawn out of the total number of plots.
#' @param lmean.k  log of expected value of the dispersion
#'     parameter of the Negative binomial for each species in the
#'     rad. Usually estimated from a linear regression of
#'     log(k)~log(abundance) from a dataset of known values of k and abundances.
#' @param lsd.k log standard deviation of the lmean.k. Can be a single
#'     value or a vector. Usually the standard error from a a linear regression of
#'     log(k)~log(abundance) from a dataset of known values of k and
#'     abundances.
#' @param lmean.sd  log of expected value of the standard deviation of
#'     the estimated population sizes of each species in the
#'     rad. Usually estimated from a linear regression of
#'     log(sd)~log(abundance) from a dataset of known values of
#'     estimated abundances and its standard deviations.
#' @param lsd.sd log standard deviation of the lmean.sd. Can be a
#'     single value or a vector. Usually the standard error from a
#'     linear regression of log(sd)~log(abundance) from a dataset of
#'     known values of estimated abundances and its standard
#'     deviations.
#' @param nrep number of repetitions of the simulated sampling
#' 
#' @return A dataframe with the vectors of the expected abundances of
#'     the 1nd, 2nd, ... Nth species without and with estimation
#'     errors of population sizes.
NB.samp <- function(rad, tot.area, n.plots, lmean.k, lsd.k, lmean.sd,
                    lsd.sd, nrep = 1){
    index <- order(rad, decreasing=TRUE)
    rad <- rad[index]
    lmean.sd <- lmean.sd[index]
    m1 <- m2 <- matrix(0,nrow=length(rad), ncol=nrep)
    for(j in 1:nrep){
        ## Samples aggregation parameters for each species
        k1 <- exp(rnorm(length(rad), mean=lmean.k, sd=lsd.k))
        y1 <- mapply(sim.pres, mu = rad/tot.area, size = k1,
                     MoreArgs=list(N = n.plots, pois.samp=FALSE))
        ## Builds simulated rad with observed abundances
        m1[1:sum(y1),j] <- rad[y1>0]
        ## Sample standard deviations for abundances estimates
        sd1 <- exp( rnorm ( sum(y1), lmean.sd[y1>0], sd = lsd.sd) )
        ## Builds simulated rad with simulated abundance estimates
        ## with Gaussian estimation error. 
        m2[1:sum(y1),j] <- rnorm(sum(y1), mean = rad[y1>0], sd1)
    } 
    data.frame( no.est.error = apply(m1,1,mean), with.est.error = apply(m2,1,mean))
}

#' Generates a community RAD and then simulates Poisson and NB samples
#' from it to define unobserved abundances (to be used in ABC).
#'
#' @details This function generates a logseries, truncated negative
#'     binomial or lognormal species abundance distribution (SAD).  The
#'     expected number of individuals of each species of the community
#'     is calculated from the values of species richness (S) and total
#'     number of individuals (N) provided for logseries, plus
#'     additional parameters for the other two SADs models.  For
#'     truncated negative binomial (tnb), the user should supply a fit
#'     of tnb to an empirical vector of abundances, usually from a
#'     sample of the community to be simulated. The parameters to
#'     simulate the abundances of the theoretical community are
#'     calculated from this object. For lognormal the user should
#'     supply the parameter 'sdlog' of this distribution model.
#' @param S positive integer, total number of species in the community
#'     to be sampled. 
#' @param N positive integer, total number of individuals in the
#'     community to be sampled. 
#' @param sad character, the name of the theoretical distribution
#'     model for the RAD of the community. Currently logseries ("ls"),
#'     truncated negative binomial ("tnb") or lognormal ("lnorm").
#' @param tot.area positive real, total area coverede by the community.
#' @param n.plots positive integer, number of sampling unities (plots)
#' of one unity of area that is drawn from the community to make the
#' sample.
#' @param nb.fit fitsad object, fit of the negative binomial model of
#'     SADs truncated at zero to a vector of species abundances in a
#'     empirical sample.
#' @param ...  further arguments to be passed to the functions called
#'     internally. Should include a named argument 'sdlog', with the
#'     value of the standard deviation of log values of abundances for
#'     the lognormal model of abundance distributions, if sad = lnorm.
#' @return a vector with the abundance of each of the S species
#'     according to the SAD model chosen by argument 'sad'
sim.rad <- function(S, N, sad=c("ls","tnb","lnorm"), nb.fit, ...){
    dots <- list(...)
    if(!is.null(nb.fit)&class(nb.fit)!= "fitsad")
        stop("nb.fit should be an object of class fitsad")
    sad <- match.arg(sad)
    if(sad=="ls"){
        ## Calculate alpha
        alpha <- fishers.alpha(N, S)
        ## Generate rad
        rad <- rad.ls(S, N, alpha, ...)$y
    }
    else
        if(sad=="tnb") {
            S.obs <- length(nb.fit@data$x)
            cf <- coef(nb.fit)
            k <- cf["size"]
            csi.p <- cf["mu"]/sum(cf)
            prob <- tovo.Scsi(S, S.obs, k, csi.p, prob = TRUE)
            rad <- rad.posnegbin(S = S, size = k, prob = prob, ...)$y
        }
    else
        if(sad=="lnorm"){
            if(!"sdlog" %in% names(dots)) stop("please provide the sdlog parameter of the lognormal RAD, as named argument 'sdlog' ")
            sdlog <- dots[["sdlog"]]
            meanlog <- log(N/S) - sdlog^2/2
            rad <- radpred(sad = "lnorm",
                           coef = list(meanlog = meanlog, sdlog = sdlog),
                           S = S, N = N)[,2]
        }
    return(rad)
}

#' simulates Poisson and NB samples from species abundance distribution.
#'
#' @details This function simulates random and clumped samples from a
#'     a vector of expected species abundances in the community. The
#'     log of values of the dispersion parameter over the plots
#'     (argument 'lmean.k' in "NB.samp") are drawn from Gaussian the
#'     estimates from a linear regression (in log scale) of value of
#'     the aggregation parameter as a function of the number of
#'     individuals per sampling unit in real data (argument
#'     'lm.k.fit'). This approach requires a sample of a community
#'     (presumably the same to be simulated) from which the
#'     aggregation parameter of each species has been estimated by
#'     fitting a negative binomial distribution. Usually there is a
#'     positive linear relationship between the standard deviations
#'     and estimated population sizes, and also between dispersion
#'     parameter and the expected abundance of each species in the
#'     sample.
#'
#' @param rad a vector of abundances of species to be sampled
#' @param tot.area positive real, total area coverede by the community.
#' @param n.plots positive integer, number of sampling unities (plots)
#' of one unity of area that is drawn from the community to make the
#' sample.
#' @param lmk.fit lm object, fit of a linear regression of the log of
#'     the aggregation parameter each species over plots (k) as a
#'     function of the log of mean abundance of the species per plot. 
#'     Data from this regression usually comes from an empiriccal
#'     sample of plots from a real community (see details).
#' @param nb.fit fitsad object, fit of the negative binomial model of
#'     SADs truncated at zero to a vector of species abundances in a
#'     empirical sample.
#' @param ...  further arguments to be passed to the functions called
#'     internally. Should include a named argument 'sdlog', with the
#'     value of the standard deviation of log values of abundances fro
#'     the lognormal model of abundance distributions, if sad = lnorm.
#' @param nrep positive integer, number of repetitions of the simulation.
#' @param summary logical, should a summary table of statistics of the simulated samples be returned?
#' @return a list with the simulated abundances of species in the Poisson and Negative binomial samples.
#' }
sim.radsamp<- function(rad,
                    tot.area, n.plots,
                    lmk.fit, nb.fit, ...){
    rad <- sort(rad[rad>0])
    S <- length(rad)
    ## Calculate expected k for each species in rad
    rad.lk <- predict(lmk.fit, newdata=data.frame(dens.ha=rad/tot.area))
    ## standard deviation of k (from regression object)
    rad.lsk <- summary(lmk.fit)$sigma
    ## simulates a value of k for each species
    rad.k <- exp(rnorm(S, mean = rad.lk, sd = rad.lsk))
    ## Poisson sample
    p.samp <- apply(matrix(rpois(S*n.plots, lambda = rad/tot.area),
                           nrow = S), 1, sum)
    ## NB sample
    nb.samp <- apply(matrix(rnbinom(S*n.plots, mu = rad/tot.area,
                                    size = rad.k), nrow = S), 1, sum)
    ## Returns a list with both types of sampling
    list(rnd.samp = p.samp, clump.samp = nb.samp)    
}

#' Generates a community RAD and then simulates Poisson and NB samples
#' from it to define unobserved abundances (to be used in ABC).
#'
#' @details This function simulates random and clumped samples from a
#'     community that follows a theoretical model of SAD (currently
#'     logseries, truncated negative binomial and lognormal).  The
#'     expected number of individuals of each species of the community
#'     is calculated from the values of species richness (S) and total
#'     number of individuals (N) provided for logseries, plus
#'     additional parameters for the other two SADs models.  For
#'     truncated negative binomial (tnb), the user should supply a fit
#'     of tnb to an empirical vector of abundances, usually from a
#'     sample of the community to be simulated. The parameters to
#'     simulate the abundances of the theoretical community are
#'     calculated from this object. For lognormal the user should
#'     supply the parameter 'sdlog' of this distribution model.  The
#'     vector of expected species abundances in the community is then
#'     sampled randomly and with clumping. Once the community
#'     abundance distribution is created, functions 'Pois.samp' and
#'     'NB.samp' are applied to simulate the distribution of estimated
#'     total abundances of species in the community that has been
#'     recorded in sample with random or clumped distribution of
#'     individuals (see help of these functions for further details).
#'     The values of the log of the standard deviations of population
#'     estimates (argument 'lmean.sd' in 'Pois.samp' and 'NB.samp')
#'     are drawn from a Gaussian using the estimates of a linear
#'     regression (in log scale) of standard deviations in function of
#'     estimated population sizes (argument 'lm.sd.fit'). The log of
#'     values of the dispersion parameter over the plots (argument
#'     'lmean.k' in "NB.samp") are also drawn from Gaussian the
#'     estimates from a linear regression (in log scale) of value of
#'     the aggregation parameter as a function of the number of
#'     individuals per sampling unit in real data (argument
#'     'lm.k.fit'). This approach requires a sample of a community
#'     (presumably the same to be simulated) from which the
#'     aggregation parameter of each species has been estimated by
#'     fitting a negative binomial distribution. Usually there is a
#'     positive linear relationship between the standard deviations
#'     and estimated population sizes, and also between dispersion
#'     parameter and the expected abundance of each species in the
#'     sample.
#' @param S positive integer, total number of species in the community
#'     to be sampled. Ignored if 'rad' is provided.
#' @param N positive integer, total number of individuals in the
#'     community to be sampled. Ignored if 'rad' is provided.
#' @param rad a vector of abundances of species to be sampled.
#' @param sad character, the name of the theoretical distribution
#'     model for the RAD of the community. Currently logseries ("ls"),
#'     truncated negative binomial ("tnb") or lognormal ("lnorm").
#' @param tot.area positive real, total area coverede by the community.
#' @param n.plots positive integer, number of sampling unities (plots)
#' of one unity of area that is drawn from the community to make the
#' sample.
#' @param lm.sd.fit lm object, fit of a linear regression of the log
#'     of standard deviation of estimates of population sizes as a
#'     function of the estimated values. Data from this regression
#'     usually comes from an empirical set of total population sizes
#'     of some species and the associated standard deviation of each
#'     estimate (see details).
#' @param lmk.fit lm object, fit of a linear regression of the log of
#'     the aggregation parameter each species over plots (k) as a
#'     function of the log of mean abundance of the species per plot. 
#'     Data from this regression usually comes from an empiriccal
#'     sample of plots from a real community (see details).
#' @param nb.fit fitsad object, fit of the negative binomial model of
#'     SADs truncated at zero to a vector of species abundances in a
#'     empirical sample.
#' @param ...  further arguments to be passed to the functions called
#'     internally. Should include a named argument 'sdlog', with the
#'     value of the standard deviation of log values of abundances fro
#'     the lognormal model of abundance distributions, if sad = lnorm.
#' @param nrep positive integer, number of repetitions of the simulation.
#' @param summary logical, should a summary table of statistics of the simulated samples be returned?
#' @return If summary = TRUE a data frame with the following summary statistics for the
#'     simulated RADs of abundances of species recorded by the random
#'     and clumped sample, without (1) or with (2) simualted
#'     estimation errors of the total population sizes (see help of
#'     Pois.samp and NB.samp):
#' \itemize{
#' \item S1, S2: number of recorded species
#' \item D1, D2: Simpson's species equivalent for the recorded RAD (that is, the
#'     inverse of Simpson index of equitability)
#' \item lmean1, lmean2: mean of the logarithm of the recorded abundances
#' \item lsd1, lsd2: standard deviation of the recorded abundances
#' If summary = FALSE a list with the simulated abundances taken from Poisson and Negative binomial samples.
#' }
sim.abc <- function(S, N, rad, sad=c("ls","tnb","lnorm"),
                    tot.area, n.plots,
                    lm.sd.fit, lmk.fit, nb.fit,
                    nrep = 1, summary = TRUE, ...){
    if(missing(rad))
        rad <- sim.rad(S, N, sad, nb.fit, ...)
    ## Calculate sd of pop estimates for each species in rad
    rad.lmean.sd <- predict(lm.sd.fit, newdata=data.frame(population=rad))
    ## standard deviation of sd (from regression object)
    rad.lsd.sd <- summary(lm.sd.fit)$sigma
    ## Calculate k for each species in rad
    rad.lk <- predict(lmk.fit, newdata=data.frame(dens.ha=rad/tot.area))
    ## standard deviation of k (from regression object)
    rad.lsk <- summary(lmk.fit)$sigma
    ## Poisson sample
    p.samp <- Pois.samp(rad = rad, tot.area = tot.area,
                        n.plots = n.plots,
                        lmean.sd = rad.lmean.sd, lsd.sd = rad.lsd.sd,
                        nrep = nrep)
    ## NB sample
    nb.samp <- NB.samp(rad = rad, tot.area = tot.area,
                       n.plots = n.plots, nrep=nrep,
                       lmean.sd = rad.lmean.sd, lsd.sd = rad.lsd.sd,
                       lmean.k = rad.lk, lsd.k = rad.lsk)
    lista <- list(rnd.samp = p.samp, clump.samp = nb.samp)
    ## Summary statistics
    if(summary){
        results <- data.frame(
            S1 = sapply(lista, function(x) sum(x[,1]>0)),
            D1 = sapply(lista, function(x) D(x[,1])),
            lmean1 = sapply(lista, function(x) mean(log(x[,1][x[,1]>0]))),
            lsd1 = sapply(lista, function(x) sd(log(x[,1][x[,1]>0]))),
            S2 = sapply(lista, function(x) sum(x[,2]>0)),
            D2 = sapply(lista, function(x) D(x[,2])),
            lmean2 = sapply(lista, function(x) mean(log(x[,2][x[,2]>0]))),
            lsd2 = sapply(lista, function(x) sd(log(x[,2][x[,2]>0])))
        )
        rownames(results) <- c("Random", "Clumped")
        return(results)
    }
    else
        return(lista)
        
}

#' utility function: Simpson's Species-equivalent
D <- function(x){
    y <- x/sum(x)
    1/sum(y^2)
    }

#' utility function: mean-square of non-zero log values
MS <- function(x, obs){
    L <- min(length(obs), sum(x>0))
    SS <- (log(x[1:L])-log(obs[1:L]))^2
    mean(SS)
    }

### Shen & He functions ###

#'normalizing factor of Shen & He equation 4 
K <- function(alpha, beta, T){
    A <- (lgamma(alpha)+lgamma(beta))-lgamma(alpha+beta)
    B <- (lgamma(alpha)+lgamma(beta+T))-lgamma(alpha+beta+T)
    1/(exp(A) - exp(B))
    }

#'He & Shen eq. 4, rho for x = 0
rho0 <- function(alpha, beta, t, T, log=FALSE){
        A <- (lgamma(alpha)+lgamma(t+beta))-lgamma(t+alpha+beta)
        B <- (lgamma(alpha)+lgamma(T+beta))-lgamma(T+alpha+beta)
        y <- log(K(alpha, beta, T)) + log(exp(A) - exp(B))
    if(log)
        return(y)
    else
        return(exp(y))
}

#' 'He & Shen eq. 4, rho for x > 0
rho1 <- function(x, alpha, beta, t, T, log=FALSE){
    A <- (lgamma(x+alpha)+lgamma(t+beta-x)) - lgamma(t+alpha+beta)
    y <- lchoose(t,x) + A
    if(log)
        return(y + log(K(alpha, beta, T)))        
    else
        return(K(alpha, beta, T) * exp(y))
}

#' Estimate species richness, unconditional Likelihood
#' @param Y a dataframe with a column with occurrence frequencies (1,
#'     2, ... n) and the other column as the number of species with
#'     each occurrence frequency. 
#' @param t number of plots in the sample
#' @param T total number of plots in the area
shen.S <- function(Y, t, T, ...){
    Y <- Y[Y[,1]>0,]
    D <- sum(Y[,2])
    f1 <- function(lS, lalpha, lbeta){
        S <- exp(lS)
        alpha <- exp(lalpha)
        beta <- exp(lbeta)
        -(
            lfactorial(S) - ( lfactorial(S-D) + sum(lfactorial(Y[,2])) ) +
            (S-D)*rho0(alpha, beta, t=t, T=T, log=TRUE) +
            sum ( Y[,2]*rho1(Y[,1], alpha, beta, t=t, T=T, log=TRUE) )
        )
    }
    mle2(f1, ...)
}


#' Maximum likelihood estimation of Shen & He alpha and beta parameters
shen.ab <- function(Y, t, T, ...){
    D <- sum(Y[,2])
    f1 <- function(lalpha, lbeta){
        alpha <- exp(lalpha)
        beta <- exp(lbeta)
        r0 <- rho0(alpha, beta, t=t, T=T, log=TRUE)
        -(
            (lfactorial(D) - sum(lfactorial(Y[,2]))) +
            sum ( Y[,2] * log(rho1(Y[,1], alpha, beta, t=t, T=T) / (1-r0)) )
        )
    }
     mle2(f1, ...)       
    }

#' Boostrap of She & He estimates
shen.boot <- function(mu, lmean.k, lsd.k, n.samp, N.tot, pois.samp=TRUE, nrep = 100, ...) {
    if(pois.samp)
            y <- sim.occ(mu = rep(mu, nrep), N = n.samp, pois.samp=TRUE)
    else{
        size <- exp(rnorm(length(lmean.k)*nrep, lmean.k, lsd.k))
        y <- sim.occ(mu = rep(mu, nrep), size = size,
                     N = n.samp, pois.samp=FALSE)
    }
    dim(y) <- c(length(mu), nrep)
    Y <- apply(y, 2, function(x) data.frame(table(x)))
    f1 <- function(x){
        x[,1] <- as.integer(as.character(x[,1]))
        return(x)
    }
    Y <- lapply(Y, f1)
    ## Estimate with unconditional likelihood (Eq.3)
    ## restricted to species richness between 1e4 and 2e4
    f2 <- function(Y) {
        z <- try(
            shen.S( Y = Y, t = n.samp, T = N.tot, ...)
        )
        if(class(z)!="try-error")
            coef(z)
        else
            rep(NA, 3)
    }
    t1 <- sapply(Y, f2)
    list(summary = c(mean = mean(exp(t1[1,]), na.rm=TRUE),
         ic.low = quantile(exp(t1[1,]), na.rm=TRUE, probs = 0.025),
         ic.up = quantile(exp(t1[1,]), na.rm=TRUE, probs = 0.975),
         N = sum(!is.na(exp(t1[1,])))),
         boot = t1 )
}


#' Shen & He Estimate species richness, conditional Likelihood
shen.S2 <- function(D, alpha, beta, t, T){
    A <- (lgamma(alpha+beta)-lgamma(beta)) + (lgamma(T+beta)-lgamma(T+alpha+beta))
    B <- (lgamma(alpha+beta)-lgamma(beta)) + (lgamma(t+beta)-lgamma(t+alpha+ beta))
    D * ((1-exp(A))/(1-exp(B)))
    }

#'Ulrich & Ollik estimates of species richness
#' @param x vector of species abundances in the sample
#' @param x.sd vector of standard deviations of species abundances,
#'     needded if 'bootstrap = TRUE'.
#' @param lm.sd.fit lm object, fit of a linear regression of the log
#'     of standard deviation of estimates of population sizes as a
#'     function of log of the estimated values. Needded if 'bootstrap
#'     = TRUE' and 'x.sd' is not provided.
#' @param effort sampling effort, that is, the fraction of the total
#'     area or total number of individuals included in the sample.
#' @param boot logical, should boostrap confidence intervals of
#'     estimated species richness be calculated?
#' @param n.boot number of boostrap samples to calculate CI's.
ulrich <- function(x, x.sd, lm.sd.fit, effort=1, boot=FALSE, n.boot = 100){
    x <- x[x>0]
    ## rad
    x.rad <- rad(x)
    ## Linear regression through central 50% quantiles of the RAD
    p.lm <- lm(log(abund)~rank, data=data.frame(x.rad), 
               subset=rank>max(rank)*.25&rank<max(rank)*.75)
    ## Regression coefficients
    cf.p.lm <- unname(coef(p.lm))
    ## Constant "d"
    d <- log(max(x))-cf.p.lm[1]
    ## Logseries projection (upper bound)
    S.reg1 <- abs((cf.p.lm[1]+log(1/effort))/cf.p.lm[2])
    ## Lognormal projection (lower bound)
    S.reg2 <- abs((2*cf.p.lm[1] +
                   log(max(x)/effort)-2*log(max(x)))/cf.p.lm[2])
    S <- data.frame(estimate = c(S.reg1,S.reg2), boot.mean = NA, boot.CI.low = NA, boot.CI.up = NA)
    rownames(S) <- c("LSE", "LNE")
    if(boot){
        if( missing(x.sd) & missing(lm.sd.fit) )
            stop("To run boostrap please provide 'x.sd' or 'lm.sd.fit'")
        else
            if(missing(x.sd)){
                ## Calculate sd of pop estimates for each species in rad
                lmean.sd <- predict(lm.sd.fit, newdata=data.frame(population=x))
                ## standard deviation of sd (from regression object)
                lsd.sd <- summary(lm.sd.fit)$sigma
                ## Samples abundances sd's from a Gaussian for each bootstrap simulation
                x.sd <- exp(rnorm(length(x)*n.boot, mean = lmean.sd, sd = lsd.sd))
                dim(x.sd) <- c(length(x), n.boot)
                sims <- matrix( nrow=length(x), ncol = n.boot)
                for(i in 1:n.boot)
                    sims[,i] <- rnorm(length(x), mean = x, sd = x.sd[,i])
            }
        else{
                sims <- rnorm(length(x)*n.boot, mean = x, sd = x.sd)
                dim(sims) <- c(length(x), n.boot)
            }
        b1 <- apply(sims, 2, function(x) ulrich(x,  boot=FALSE)$S[,1])
        S[1,2:4] <- c(mean(b1[1,]), quantile(b1[1,], probs=c(0.025, 0.975)))
        S[2,2:4] <- c(mean(b1[2,]), quantile(b1[2,], probs=c(0.025, 0.975)))
    }
    return(
        list(S = S, coefs=c(coef(p.lm), d=d) , lm.fit = p.lm)
        )
}

#' Bootstrap mean and IC of Ulrich & Ollik estimates of species richness
#' @param ... arguments to be passed to Pois.samp (if pois.samp = TRUE) or NB.samp (if pois.samp=FALSE)
#' @param pois.samp logical,  if TRUE simulates a Poisson sample of the regional RAD;
#'     simulates a negative binomial sample with parameters passed to 'NB.samp' otherwise.
#' @param nrep positive integer, number of boostrap simulations to be done.
ulrich.boot <- function(rad, tot.area, n.plots, lm.sd.fit, lm.k.fit,  pois.samp=TRUE, nrep = 100){
    x <-  rad[rad>0]
    y <- vector(mode = "list", length = nrep)
    lmean.sd <- predict(lm.sd.fit, newdata = data.frame(population = x))
    lsd.sd <- summary(lm.sd.fit)$sigma
    if(pois.samp)
    {
        for(i in 1:nrep)
            y[[i]] <- Pois.samp(x, tot.area, n.plots, lmean.sd, lsd.sd)
    }
    else
    {
        lmean.k <- predict(lm.k.fit, newdata = data.frame(dens.ha = x/tot.area))
        lsd.k <- summary(lm.k.fit)$sigma
        for(i in 1:nrep)
            y[[i]] <- NB.samp(x, tot.area, n.plots, lmean.k, lsd.k, lmean.sd, lsd.sd)
    }
    #browser()
    t1 <- sapply(y, function(x) ulrich(x[,2], lm.sd.fit = lm.sd.fit, boot=TRUE)$S[1,2])
    list( summary = c(mean = mean(t1),
                      ic.low = quantile(t1, 0.025),
                      ic.up = quantile(t1, 0.975)),
                      boot = t1 )
}



#' Hui ORC model
#' @param occupancies observed occupancies frequencies in a samples
#' @param effort sampling effort, that is, the fraction of the total
#'     area or total number of individuals included in the sample. 
hui.orc <- function(occupancies, effort=1){
    x <- data.frame(rad(occupancies))
    m1 <- lm(log(abund) ~ rank + log(rank), data = x)
    cf <- unname(coef(m1))
    C1 <- cf[1]-log(effort)
    S.est <- cf[3]*lambertW(cf[2]*exp(-C1/cf[3])/cf[3])/cf[2]
    list(S.est = S.est, model = m1)
    }

#' Parametric boostrap for Hui ORC model
#' @param lmean.k  log of expected value of the dispersion
#'     parameter of the Negative binomial for each species in the
#'     rad. Usually estimated from a linear regression of
#'     log(k)~log(abundance) from a dataset of known values of k and abundances.
#' @param lsd.k log standard deviation of the lmean.k. Can be a single
#'     value or a vector. Usually the standard error from a a linear regression of
#'     log(k)~log(abundance) from a dataset of known values of k and
#'     abundances.
hui.boot <- function(mu, lmean.k, lsd.k, n.samp, N.tot, pois.samp=TRUE, nrep = 100){
    effort <- n.samp/N.tot
    if(pois.samp)
        y <- sim.occ(mu = rep(mu, nrep), N = n.samp, pois.samp=TRUE)
    else{
        size <- exp(rnorm(length(lmean.k)*nrep, lmean.k, lsd.k))
        y <- sim.occ(mu = rep(mu, nrep), size = size,
                     N = n.samp, pois.samp=FALSE)
    }
    dim(y) <- c(length(mu), nrep)
    t1 <- apply(y, 2, function(x) hui.orc(x, effort = effort)$S.est)
    list( summary = c(mean = mean(t1),
                      ic.low = quantile(t1, 0.025),
                      ic.up = quantile(t1, 0.975)),
                      boot = t1 )
}

#' Calulates unbiased men and CIs from the bias simulations
#'
#' @param object list, an object returned from the bias analyses (see scripts in directory 'bias estimation')
#' @param ci.vector, vector of lower and upper limits of original CI of the estimate.
#'
#' @return a data.frame with the mean, and limits of empirical CI corrected from bias between the estimated value and true value,
#' assuming random and clumped sampling.
bias.ci <- function(object, ci.vector){
    ci <- unlist(sort(ci.vector))
    y1 <- with(object$estimates,
               S[S.est.rnd > ci[1] & S.est.rnd < ci[2]])
    y2 <- with(object$estimates,
               S[S.est.clump > ci[1] & S.est.clump < ci[2]])
    lista <- list(y1,y2)
    data.frame(mean = sapply(lista,mean, na.rm =TRUE),
               IC.low = sapply(lista, quantile, 0.025, na.rm=TRUE),
               IC.up = sapply(lista, quantile, 0.975, na.rm=TRUE),
               row.names = c("rnd.samp", "clump.samp"))
    }

#' A function to rule them all
#'
#' Run all the calculations used in Steege et al 2019
#'
#' @details this is just a concatenation of commands that generate all
#'     objects necessary to reproduce the results.  This is not a
#'     generic function and works properly only with the data
#'     structures and functions prepared for the paper, with many
#'     idiosyncrasies.
#' 
#' @param path.to.data character, name of a csv file with data of
#'     plots. Must be a csv with a species per row, and columns named
#'     'species'(species names), 'population' (estimate of total
#'     population sizes), 'N.ind' (number of individuals recorded in
#'     the sample), 'N.plots' (number of plots where the species has
#'     been recorded), 'pop.sd' (standard deviation of estimates of
#'     total population sizes).
#' @param Tot.t positive integer, estimated total number of trees in
#'     the whole region from which the sample of plots came from
#'     (e.g. Whole Amazon)
#' @param Tot.A positive real, total area of the whole region from
#'     which the sample of plots came from.
#' @param N.plots positive integer, Number of plots in the sample
#'     taken from the region.
#' @param Samp.A positive real, area sampled by the set of plots (that
#'     is, sum of plot areas)
#' @param lm.sd.fit lm object, fit of a linear regression of the log
#'     of standard deviation of estimates of population sizes as a
#'     function of log of the estimated values. Ignored if data has a vector
#'     names 'pop.sd'.
#' @return a list with objects needed to reproduce the results in the
#'     paper (some of them are redundant, to be simplified).
atdn.estimates <- function(path.to.data,
                           Tot.t, Tot.A = 6.29e8,
                           N.plots, Samp.A, lm.sd.fit) {
    ## ----Data prep-----------------------------------------------------------
    ## loads data
    data <- read.csv2(path.to.data, as.is=TRUE)
    ## Derived quantities ##
    ## Proportion of total trees in the sample
    p1 <- sum(data$N.ind)/Tot.t
    ## Number of individuals in the sample
    N.ind <- data$N.ind
    ## Number os species in the sample
    Sobs <- length(N.ind)

    cat("\n----------------------------------------------------------------------\n
                         Fitting SADs models to plot data \n
         ----------------------------------------------------------------------\n\n")
    
    ## ----fit pln-------------------------------------------------
    pln <- fitpoilog(N.ind)

    ## ----PLN species richness estimate 2-------------------------------------
    pln.cf <- coef(pln)
    pln.d0 <- dpoilog(0, mu = pln.cf[1], sig=pln.cf[2])

    ## ----fit ls--------------------------------------------------------------
    y.ls <- fitls(N.ind)

    ## ----estimated S ls------------------------------------------------------
    alpha <- coef(y.ls)[[2]]
    S.ls <- alpha*log(1 + Tot.t/alpha)

    ## ----ls and S est for ls-------------------------------------------------
    ls.ci <- confint(y.ls)
    S.ls.ci <- c(ls.ci[1]*log(1 + Tot.t/ls.ci[1]), ls.ci[2]*log(1 + Tot.t/ls.ci[2]))
    ## ----fit NB--------------------------------------------------
    y.nb2 <- fitnbinom(N.ind, 
                       start.value=c(size=0.3, mu=mean(N.ind)))

    ## ----Tovo S estimate-----------------------------------------------------
    cf.nb <- coef(y.nb2)
    csi.p <- unname(cf.nb[2]/(sum(cf.nb)))
    csi <- csi.p/(p1+(1-p1)*csi.p)
    ## Estimated number of species 
    S.nb <- Sobs*(1-(1-csi)^cf.nb[1]) / (1-(1-csi.p)^cf.nb[1])
    S.nb <- unname(S.nb)

    ## ----NB S est CI---------------------------------------------------------
    tovo.S <- tovo(fit = y.nb2, p = p1, CI=TRUE)


    cat("\n----------------------------------------------------------------------\n
                       Estimates from total population sizes\n
         ----------------------------------------------------------------------\n\n")

    ## ---- regression of sd x estimated population  -----------------
    
    if("pop.sd" %in% names(data))
        lm.sd <- lm( log(pop.sd) ~ log(population), data = data)
    else
        lm.sd <- lm.sd.fit
    
    lm.sd.sigma <- summary(lm.sd)$sigma
    
    ## ----Linear extrapolation from RAD  of estimated population sizes (LSE) ------------------------------
    if("pop.sd" %in% names(data))
        S.ulrich <- ulrich(data$population, x.sd = data$pop.sd, boot = TRUE, n.boot = 200)
    else
        S.ulrich <- ulrich(data$population, lm.sd.fit = lm.sd, boot = TRUE, n.boot = 200)
   
    ## ----Amazon S and alpha from LSE ---------------------------------------
    S.r.ls <- S.ulrich$S[1,1]
    alpha.r <- fishers.alpha(N = Tot.t, S = S.r.ls)

    ## ----LSE regional rad-------------------------------------------------------
    reg.ls.rad <- ceiling(
        rad.ls(S = S.r.ls, N = Tot.t, alpha = alpha.r)$y
    )

    ## ----TNB regional RAD-----------------------------------------------------
    reg.nb.rad <- rad.posnegbin(S = S.nb, size = cf.nb[1], 
                                prob = 1-csi)$y

    ## ----k x dens regression-------------------------------------
    ## estimating k parameter of a NB for each species 
    data$dens.ha <- data$N.ind/Samp.A
    data$k <- est.kv(mu=data$dens.ha, 
                     nzeroes=N.plots-data$N.plots, 
                     Nplots=N.plots)
    lm.k <-lm(log(k)~log(dens.ha), 
              data=data, subset=k<1)
    ## Estimated regression standard error
    lm.k.sigma <- summary(lm.k)$sigma

    ## ----sim popsizes LS and NB----------------------------------
    ## Predicted log(sd) values of estimated population sizes for LS rad
    reg.ls.rad.lsd <- predict(lm.sd, 
                             newdata=data.frame(population = reg.ls.rad))
    ## Predicted log(sd) values for TNB rad
    reg.nb.rad.lsd <- predict(lm.sd, 
                             newdata=data.frame(population=reg.nb.rad))
    ## Predicted log(k) values for LS rad
    reg.ls.rad.lk <- predict(lm.k, 
                             newdata=data.frame(dens.ha=reg.ls.rad/Tot.A))
    ## Predicted log(k) values for TNB rad
    reg.nb.rad.lk <- predict(lm.k, 
                             newdata=data.frame(dens.ha=reg.nb.rad/Tot.A))
    
    
    cat("\n----------------------------------------------------------------------\n
               Estimates from occupancy data (Shen & He and Hui estimators)\n
         ----------------------------------------------------------------------\n\n")

    ## ----shen estimate, eval=FALSE-------------------------------------------
    ## table of frequencies of occurrences
    Y <- data.frame(table(data$N.plots))
    Y[,1] <- as.integer(as.character(Y[,1]))
    ## Estimate of alfa and beta (Eq.6)
    ## To be use as starting values for the unconditional estimation below
    ##ab.est <- shen.ab(Y = Y, t = N.plots, T = Tot.A,
    ##                  start=list(lalpha=-2, lbeta=1), method="SANN")
    ## estimated coeficients
    ##cf.st1 <- coef(ab.est)
    ## Estimate with uncoditional likelihood (Eq.3)
    ## restricted to species richness between 1e4 and 2e4
    Shen <- #try(
        shen.S( Y = Y, t = N.plots, T = Tot.A,
               ##start=c(list(lS = log(S.ls)), as.list(cf.st1)),
               start=list(lS = log(1.5e4), lalpha=-2, lbeta=1),
               method="L-BFGS-B",
               upper=c(lS=log(2e4), lalpha=Inf, lbeta=Inf),
               lower=c(lS=log(1e3), lalpha=-Inf, lbeta=-Inf))
    #)
    ##if(class(Shen)!="try-error"){
    cf.st2 <- coef(Shen)
    S.shen <- unname(exp(cf.st2[1]))
        ## ----Shen He profile-----------------------------------------------------
        ## Shen.prf <- profile(Shen, which=1)
    ##    }
    
    
    ## ----Hui ORC estimate----------------------------------------------------
    S.orc <- hui.orc(data$N.plots, effort=Samp.A/Tot.A)
    orc.cf <- coef(S.orc$model)
    

    ##------------------------------------------------------------------------
    ## Save all objects generated by this function in a list
    lista <- mget(ls())
    return(lista[sapply(lista, function(x) class(x)!="function")])                                    
}

