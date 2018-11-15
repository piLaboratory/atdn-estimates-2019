library(copula)
library(VineCopula)
library(sads)
library(VGAM)
## 1st try with bivariate copulas to fit the joint distribution of
## SAD e SOD
Nplots <- dados$N.plots
Nind <- dados$N.ind
Np.u <- pobs(Nplots)
Ni.u <- pobs(Nind)
## Cheking the best copula model
copsel <- BiCopSelect(Np.u, Ni.u)
summary(copsel) ## t- copula

## Creating Poilog-betabinomial copula with parameters from the fit above
## and fitted marginals (Poilog and betabinom)
## Fisrt test with untruncated marginals
myMvd <- mvdc(copula = tCopula(0.77, df=5.87),
              margins = c("poilog", "betabinom"),
              paramMargins = list(list(mu = pln.cf[1], sig = pln.cf[2]), list(mu = bb.mle[1], rho = bb.mle[2])))
plot(Nind, Nplots) ## na escala linear eh um problema
plot(Nind, Nplots, log="xy")
## start parameters
a.0 <- sin(cor(Nind, Nplots, method = "kendall") * pi/2)
cfit1 <- fitMvdc(cbind(Nind,Nplots), myMvd,
                 start=c(pln.cf,bb.mle2, a.0, 5)) ## not converging
