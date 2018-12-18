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

#' Find the value of Csi in Tovo's TNB given observed and total
#' species richness
#' @param S total species richness
#' @param S.obs species richness in the sample
#' @param k estimated size parameter of the Truncated negative binomial 
#' @param csi.p Csi parameter (1 - prob, see Tovo et al) estimated for the sample
tovo.Scsi <- function(S, S.obs, k, csi.p){
    C1 <- S.obs / (1-(1-csi.p)^k)
    -((C1 - S)/C1)^(1/k) + 1
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
    pp <- rev(ppoints(S.r))
    x <- seq(1,S.r, length=npoints)
    y <- qls2(p = pp[x], N = N, alpha = alpha, ...)
    data.frame(x, y)
    }

#' Continuous approximation for quantile function for TNB distribution
qposnegbin2 <- function(p, size, prob, lower=3e-9, upper=3e9){
    f2 <- function(target){
        f1 <- function(x) pposnegbin(x,size, prob) - target
        uniroot(f1, lower=lower, upper=upper)$root
    }
    sapply(p, f2)
}

#' Zero-truncated negative binomial RAD
#' @description Generates a given number of points of the RAD of a
#'     ZTNB, given the total number of species and the parameters of
#'     the distribution.
#' @param S total species richness in the RAD
#' @param size size parameter of the Zero-truncated negative binomial
#' @param prob prob parameter of the Zero-truncated negative binomial
#' @param npoints number of points along teh RAD to retunr (defaults
#'     to the total number of species, which implies that the expected
#'     abundance will be calculated for each species (can be slow for
#'     a large number of species)
rad.posnegbin <- function(S, size, prob, npoints = round(S),...){
    S.r <- round(S)
    pp <- rev(ppoints(S.r))
    x <- seq(1,S.r, length=npoints)
    y <- qposnegbin2(p = pp[x], size = size, prob = prob, ...)
    data.frame(x, y)
    }

#' Find total number of species by sampling probabilities of occurence from a beta distribution
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
## random deviates
rbetabinom.t0 <- function(n, size, mu, rho) rtrunc("betabinom", n, trunc=0, coef=list(size=size, mu=mu, rho=rho))
## Poilog truncated at zero
dpoilog.t0 <- function(x, mu, sig, log=FALSE) dtrunc("poilog", x, trunc=0, coef=list(mu=mu, sig=sig), log=log)
ppoilog.t0 <- function(q, mu, sig, ...) ptrunc("poilog", q, trunc=0, coef=list(mu=mu, sig=sig), ...)
qpoilog.t0 <- function(p, mu, sig, ...) qtrunc("poilog", p, trunc=0, coef=list(mu=mu, sig=sig), ...)
rpoilog.t0 <- function(n, mu, sig) rtrunc("poilog", n, trunc=0, coef=list(mu=mu, sig=sig))

#' testing: estimated species richnedd from zero-truncated negative binomial lognormal
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
## Vectorized version
est.kv <- Vectorize(est.k, c("mu", "nzeroes"))

#' analytical
est.ka <- function(mu, nzeroes, Nplots){
    (mu * nzeroes)/(Nplots - nzeroes)
    }
 
#'Draw N samples from Negative binomial for a pair of parameters of the NB (mu and size) and the sums up these values
rnbinom2 <- function(mu, size, N){
        y <- rnbinom(n = N, mu  = mu, size = size)
        sum(y)
        }

#'Drawn N samples from a Poisson distribution and then sums up these values
rpois2 <- function(lambda, N){
        y <- rpois(n = N, lambda  = lambda)
        sum(y)
    }

#' Simulates the number of non-occupied plots for a Negative Binomial
#' sampling or a Poisson sampling of a RAD
#' @param mu vector of species densities (individuals per plot unit) in the RAD to be sampled. 
#' @param size vector of values of parameter size of the negative
#' binomial for each species in the RAD
#' @param N number of plots to be sampled
#' @param pois.samp logical, if TRUE simulates a Poisson sample,
#' simulates a negative binomial sample with parameters k otherwise.
sim.occ <- function(mu, size, N, pois.samp=TRUE){
    if(pois.samp)
        lp <- -mu
    else
        lp <- size*(log(size)- log(mu+size))
    p0 <- exp(N*lp)
    sample(0:1, size=1, prob=c(p0, 1-p0))
}

#' Simulates samples of population sizes using a Poisson sample
#' @details This function performs a simulation of which species would
#'     be included in a Poisson sample of a regional RAD. The
#'     population sizes of the included species are then set to the
#'     population sizes in the RAD. This is a simulation of how
#'     the distribution of total population sizes would look like,
#'     assuming that the there is a method to estimate the exact
#'     population size and only source of uncertainty is which species
#'     will be included in the sample.
#' @param rad a vector with the species population sizes in the RAD to be sampled
#' @param tot.area total area of the community to be sampled. The area unit is one plot
#' @param n.plots number of sampling units (e.g. plots) to be drawn out of the total number of plots.
#' @param nrep number of repetitions of the simulated sampling
Pois.samp <- function(rad, tot.area, n.plots, nrep){
    rad <- sort(rad, decreasing=TRUE)
    m1 <- matrix(0,nrow=length(rad), ncol=nrep)
    for(j in 1:nrep){
        y1 <- mapply(sim.occ, mu = rad/tot.area,
                     MoreArgs=list(N = n.plots))
        m1[1:sum(y1),j] <- rad[y1>0]
        }
    apply(m1,1,mean)    
}

#' Simulates a Negative Binomial samples of population sizes of a  RAD
#' @param rad a vector with the species population sizes in the RAD to be sampled
#' @param tot.area total area of the community to be sampled. The area unit is one plot
#' @param n.plots number of sampling units (e.g. plots) to be drawn out of the total number of plots.
#' @param lmean.k  log of expected value of the aggregation
#'     parameter of the Negative binomial for each species in the
#'     rad. Usually estimated from a linear regression of
#'     log(k)~log(abundance) from a dataset of known values of k
#' @param lsd.k log standard deviation of the lmean.k. Can be a singel
#'     value or a vector. Usually the standard error from a a linear regression of
#'     log(k)~log(abundance) from a dataset of known values of k
#' @param nrep number of repetitions of the simulated sampling
#' @details This function performs a simulation of which species would
#'     be included in a Poisson sample of a regional RAD. The
#'     population sizes of the included species are then set to the
#'     population sizes in the RAD. This is a simulation of how
#'     the distribution of total population sizes would look like,
#'     assuming that the there is a method to estimate the exact
#'     population size and only source of uncertainty is which species
#'     will be included in the sample.
NB.samp <- function(rad, tot.area, n.plots, lmean.k, lsd.k, nrep){
    rad <- sort(rad, decreasing=TRUE)
    m1 <- matrix(0,nrow=length(rad), ncol=nrep)
    for(j in 1:nrep){
        ## Samples aggregation parameters for each species
        k1 <- exp(rnorm(length(rad), mean=lmean.k, sd=lsd.k))
        y1 <- mapply(sim.occ, mu = rad/tot.area, size = k1,
                     MoreArgs=list(N = n.plots, pois.samp=FALSE))
        m1[1:sum(y1),j] <- rad[y1>0]
        }
    apply(m1,1,mean)    
}

#' Simulates Poisson and NB samples from a LS rad for ABC
#' @param
sim.abc <- function(S, N, tot.area, n.plots, lmk,
                    nb.fit, LS=TRUE, obs.values, nrep = 2, ...){
    if(!is.null(nb.fit)&class(nb.fit)!= "fitsad")
        stop("nb.fit should be an object of class fitsad")
    if(LS){
        ## Calculate alpha
        alpha <- fishers.alpha(N, S)
        ## Generate rad
        rad <- rad.ls(S, N, alpha, ...)$y
    }
    else{
        S.obs <- length(nb.fit@data$x)
        cf <- coef(nb.fit)
        k <- cf["size"]
        csi.p <- cf["mu"]/sum(cf)
        csi <- tovo.Scsi(S, S.obs, k, csi.p)
        rad <- rad.posnegbin(S, k, 1-csi, ...)$y
        }
    ## Calculate k for each species in rad
    rad.lk <- predict(lmk, newdata=data.frame(dens.ha=rad/tot.area))
    ## standard deviation of k (from regression object)
    rad.lsk <- summary(lmk)$sigma
    ## Poisson sample
    p.samp <- Pois.samp(rad = rad, tot.area = tot.area,
                        n.plots = n.plots, nrep = nrep)
    ## NB sample
    nb.samp <- NB.samp(rad = rad, tot.area = tot.area,
                       n.plots = n.plots, nrep=nrep,
                       lmean.k = rad.lk, lsd.k = rad.lsk)
    lista <- list(p.samp, nb.samp)
    ## Summary statistics
    data.frame(
        S = sapply(lista, function(x) sum(x>0)),
        D = sapply(lista, D),
        lmean = sapply(lista, function(x) mean(log(x[x>0]))),
        lsd = sapply(lista, function(x) sd(log(x[x>0])))
    )
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
#' @param Y a dataframe with a column with occurrence frequencies (1, 2, ... n) and the other column as the number of species with each occurrence frequency
#' @param t number of plots in the sample
#' @param T total number of plots in the area
shen.S <- function(Y, t, T, ...){
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


#' likelihood for Shen & He alpha and beta parameters
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

#' Shen & He Estimate species richness, conditional Likelihood
shen.S2 <- function(D, alpha, beta, t, T){
    A <- (lgamma(alpha+beta)-lgamma(beta)) + (lgamma(T+beta)-lgamma(T+alpha+beta))
    B <- (lgamma(alpha+beta)-lgamma(beta)) + (lgamma(t+beta)-lgamma(t+alpha+ beta))
    D * ((1-exp(A))/(1-exp(B)))
    }

#'Ulrich & Ollik estimates
#' @param x vector of species abundances in the sample
#' @param effort sampling effor, that is the fraction of the total area or total number of individuals included in the sample
ulrich <- function(x, effort=1){
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
    ##S.reg2 <- abs((cf.p.lm[1]-d)/cf.p.lm[2])
    S.reg2 <- abs((2*cf.p.lm[1] + log(max(x)/effort)-2*log(max(x)))/cf.p.lm[2])    
  
    return(list(S=c(S.reg1, S.reg2), coefs=c(coef(p.lm), d=d) ))
}

#' Hui ORC model
#' @param occupancies observed occupancies frequencies in a samples
#' @param effort sampling effor, that is the fraction of the total area or total number of individuals included in the sample
hui.orc <- function(occupancies, effort=1){
    x <- data.frame(rad(occupancies))
    m1 <- lm(log(abund) ~ rank + log(rank), data = x)
    cf <- unname(coef(m1))
    C1 <- cf[1]-log(effort)
    S.est <- cf[3]*lambertW(cf[2]*exp(-C1/cf[3])/cf[3])/cf[2]
    list(S.est = S.est, model = m1)
    }
