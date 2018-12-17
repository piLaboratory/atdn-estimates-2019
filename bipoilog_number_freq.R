library(poilog)
library(sads)
pln <- fitpoilog(dados$N.ind)
fpln <- fitpoilog(dados$N.plots)
pln.cf <- coef(pln)
fpln.cf <- coef(fpln)
bipoi <- bipoilogMLE(dados$N.ind, dados$N.plots,
                     startVals=c(mu1=pln.cf[1], sig1=pln.cf[2],
                                 mu2=fpln.cf[1], sig2=fpln.cf[2],
                                 rho=cor(dados$N.ind,dados$N.plots)))
