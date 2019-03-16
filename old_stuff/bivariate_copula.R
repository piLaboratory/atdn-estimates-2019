library(copula)
library(VineCopula)
library(sads)
library(VGAM)
library(TruncatedDistributions)
## 1st try with bivariate copulas to fit the joint distribution of
## SAD e SOD
Nplots <- dados$N.plots
Nind <- dados$N.ind
Np.u <- pobs(Nplots)
Ni.u <- pobs(Nind)
## Cheking the best copula model
copsel <- BiCopSelect(Np.u, Ni.u)
(copsel.s <- summary(copsel)) ## t- copula

## Creating Poilog-betabinomial copula with parameters from the fit above
## and fitted marginals (Poilog and betabinom)
## First test with untruncated marginals
myMvd <- mvdc(copula = tCopula(param=copsel.s$par, df=copsel.s$par2),
              margins = c("poilog", "betabinom"),
              paramMargins = list(list(mu = pln.cf[1], sig = pln.cf[2]), list(mu = bb.mle[1], rho = bb.mle[2])))
plot(Nind, Nplots) ## na escala linear eh um problema
plot(Nind, Nplots, log="xy")
## start parameters
a.0 <- sin(cor(Nind, Nplots, method = "kendall") * pi/2)
cfit1 <- fitMvdc(cbind(Nind,Nplots), myMvd,
                 start=c(pln.cf,bb.mle2, a.0, 5)) ## not converging

## Same for a bivariate weibull
wtr <- 0.2
Nweib <- fitsad(dados$N.ind, "weibull", trunc=wtr)
Pweib <- fitsad(dados$N.plots, "weibull", trunc=wtr)
dtw2 <- function(x, shape, scale, log=FALSE){
    y <- dtweibull(x, shape, scale, a=wtr)
    if(log)
        return(log(y))
    else
        return(y)
}
ptw2 <- function(q, shape, scale, log=FALSE){
    y <- ptweibull(q, shape, scale, a=wtr)
    if(log)
        return(log(y))
    else
        return(y)
}
qtw2 <- function(p, shape, scale)
    qtweibull(q, shape, scale, a=wtr)
    
## First test with untruncated marginals
myMvd2 <- mvdc(copula = tCopula(param=copsel.s$par, df=copsel.s$par2),
               margins = c("tw2", "tw2"),
               paramMargins = list(as.list(coef(Nweib)), as.list(coef(Pweib))))
plot(Nind, Nplots) ## na escala linear eh um problema
plot(Nind, Nplots, log="xy")
## start parameters
a.0 <- sin(cor(Nind, Nplots, method = "kendall") * pi/2)
cfit2 <- fitMvdc(cbind(Nind,Nplots), myMvd2,
                 start=c(coef(Nweib), coef(Pweib), a.0, 5))
summary(cfit2)
(cfit2.cf <- coef(cfit2))
    
## Fitted copula
myMvd.f <- mvdc(copula = tCopula(param=cfit2.cf[5], df=cfit2.cf[5]),
               margins = c("tw2", "tw2"),
               paramMargins = list(list(shape=cfit2.cf[1],scale=cfit2.cf[2]), list(shape=cfit2.cf[3],scale=cfit2.cf[4])))
