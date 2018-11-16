library(VGAM)
library(sads)
source("functions.R")

y <- dados$N.ind
Sobs <- length(y)
## Fit of NB
## With VGAM
y.nb <- vglm(y ~ 1, posnegbinomial)
## With sads
y.nb2 <- fitnbinom(y, start.value=c(size=0.3, mu=mean(y)))
cf.nb <- coef(y.nb2)

## Parte misteriosa do codigo
## Fit of the empirical RSA at the local scale with a negative binomial distribution
## normalized from 0 to infinity
y<-as.numeric(y)
pdata <- data.frame(munb = cf.nb[2], size = cf.nb[1])
pdata <- transform(pdata, y = y)
