## ----setup------------------------------------
library(VGAM)
library(untb)
library(sads)
library(abc)
source("functions.R")

#' A function to rule them all
#'
#' Run all the calculations used in Steege et al 2019
#'
#' @details this is just a concatenation of commands that generate all objects necessary to reporduce the results.
#' This is not a generic function and works properly only with the data structures and function prepared
#' for the paper, with many idiosyncrasies.
#' 
#' @param path.to.data character, name of a csv file with data of plots. Must be a csv with a species per row,
#' and columns named 'species'(species names), 'population' (estimate of total population sizes),
#' 'N.ind' (number of individuals recorded in the sample), 'N.plots' (number of plots where the species has been recorded
#' @param path.to.abc character, name of a binary RData file with results of ABC simulation (see file abc_run.R)
#' @param Tot.t positive integer, estimated total number of trees in the whole region from which
#' the sample of plots came from (e.g. Whole Amazon)
#' @param Tot.A positive real, total area of the whole region from which the sample of plots came from.
#' @param N.plots positive integer, Number of plots in the sample taken from the region.
#' @param Samp.A positive real, area sampled by the set of plots (that is, sum of plot areas)
#' 
#' @return All objects needed to reproduce the results in the paper.
atdn.estimates <- function(path.to.data, path.to.abc, Tot.t, Tot.A = 6.29e8, N.plots, Samp.A ) {
    ## ----Data prep-----------------------------------------------------------
    data <- read.csv2(path.to.data, as.is=TRUE)
    ##loads ABC simulation data
    load(path.to.abc)
    ## Derived quantities ##
    ## Proportion of total trees in the sample
    p1 <- sum(data$N.ind)/Tot.t
    ## Number of individuals in the sample
    N.ind <- data$N.ind
    ## Number os species in the sample
    Sobs <- length(N.ind)
    
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


    ## ----Linear extrapolation from regional RAD------------------------------
    S.ulrich <- ulrich(data$population)
    S.r.ls <- S.ulrich$S[1]

    ## ----Amazon alpha--------------------------------------------------------
    alpha.r <- fishers.alpha(N = Tot.t, S = S.r.ls)

    ## ----amazon LS rad-------------------------------------------------------
    reg.ls.rad <- ceiling(
        rad.ls(S = S.r.ls, N = Tot.t, alpha = alpha.r)$y
    )


    ## ----TNB regionl RAD-----------------------------------------------------
    reg.nb.rad <- rad.posnegbin(S = S.nb, size = cf.nb[1], 
                                prob = 1-csi)$y

    ## ----nbinom rad, echo=FALSE----------------------------------------------
    cf.u <- S.ulrich$coefs

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


    ## ----abc model selection-------------------------------------------------
    ## Target: observed number of species, Simpson's 1/D, 
    ## lmean, sdmean             
    target <- c(Sobs, 
                D(data$population), 
                mean(log(data$population)), 
                sd(log(data$population)))
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
                   tol=0.01, method="rejection")
    S.post1.s <- summary(S.post1)


    ## ----clumped LS RAD and CI-----------------------------------
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

    ## ----sim popsizes LS and NB----------------------------------
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


    ## ----shen estimate, eval=FALSE-------------------------------------------
    ## table of frequencies of occurrences
    Y <- data.frame(table(data$N.plots))
    Y[,1] <- as.integer(as.character(Y[,1]))
    ## Estimate of alfa and beta (Eq.6)
    ## To be use as starting values for the unconditional estimation below
    ab.est <- shen.ab(Y = Y, t = N.plots, T = Tot.A,
                      start=list(lalpha=-5, lbeta=3), method="SANN")
    ## estimated coeficients
    cf.st1 <- coef(ab.est)
    ## Estimate with uncoditional likelihood (Eq.3)
    ## restricted to species richness between 1e4 and 2e4
    ShenHe <- shen.S( Y = Y, t = N.plots, T = Tot.A,
                     start=c(list(lS = log(1.5e4)), as.list(cf.st1)),
                     method="L-BFGS-B",
                     upper=c(lS=log(2e4), lalpha=Inf, lbeta=Inf),
                     lower=c(lS=log(1e3), lalpha=-Inf, lbeta=-Inf))
    cf.st2 <- coef(ShenHe)
    S.sh <- unname(exp(cf.st2[1]))

    ## ----Shen He confint-----------------------------------------------------
    ShenHe.prf <- profile(ShenHe, which=1)

    ## ----Hui ORC estimate----------------------------------------------------
    S.orc <- hui.orc(data$N.plots, effort=Samp.A/Tot.A)
    orc.cf <- coef(S.orc$model)

    ##------------------------------------------------------------------------
    ## Save all objects generated by this function in a list
    lista <- mget(ls())
    return(lista[sapply(lista, function(x) class(x)!="function")])                                    
}
