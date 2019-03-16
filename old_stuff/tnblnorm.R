################################################################################
## Following instructions at https://github.com/stan-dev/rstan/wiki/Installing-RStan-on-Mac-or-Linux
## To set compiler to install rstan in linux, check  if needed in oter OS (runned once)
dotR <- file.path(Sys.getenv("HOME"), ".R")
if (!file.exists(dotR)) dir.create(dotR)
M <- file.path(dotR, "Makevars")
if (!file.exists(M)) file.create(M)
cat("\nCXXFLAGS=-O3 -mtune=native -march=native -Wno-unused-variable -Wno-unused-function", 
    file = M, sep = "\n", append = TRUE)
cat("\nCXXFLAGS+=-flto -ffat-lto-objects  -Wno-unused-local-typedefs", 
    file = M, sep = "\n", append = TRUE)
################################################################################

library(rstan)
library(parallel)
library(sads)
source(/home/paulo/work/orientacoes/recentes/lima/posdoc/artigo_sads_mata_atlantica/functions.R)
#library(magritr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## 1. Truncated Negative binomial Logormal
## Data 
Y <- dados$N.ind
## List with all data need for the model
m1.data <- list(s = length(Y), y=Y, a=2.048e3/5.79e8, lN=log(567*5.5e8))
## To set a priori distributions: fit a poilog to the data
Y.pl <- fitpoilog(Y)
cf2 <- coef(Y.pl)
## The parameter sigma
cf2[2]
## The parameter mu is the estimated parameter - log(fraction sampled)
cf2[1] - log(1.9e3/5.5e8)


## Fit the model (ca. 60 min in my computer)
## Three chains in parallel, first 2000 for warming, thining = 100 steps
## Parameters tunned because default values had some convergence issues and high autocorrelation
## With these (very) conservative settings diagnostics looks ok
fit1 <- stan(file="tnblnorm.stan", data=m1.data, iter=12000,
             warmup=2000, refresh=1000, chains=3, cores=3, thin=100)
fit.nbln <- stan(file="tnblnorm.stan", data=m1.data, iter=200, chains=3, cores=3)
## Diagnostics
traceplot(fit.nbln, c("mu", "sig", "S1", "phi")) ## no convergence, weird estimates
pairs(fit1, pars=c("mu", "sig", "S1", "phi"))
summary(fit.nbln, c("mu", "sig", "S1", "phi"))$summary ## No Rhat > 1.1, effective sample size close to the true size, everything looks ok
print(fit1, c("mu", "sig", "S1", "phi"))
sampler_params <- get_sampler_params(fit1, inc_warmup = TRUE)
lapply(sampler_params, summary, digits = 2)
##  A dataframe with a posteriori distributions
fit1.df <- as.data.frame(fit1, c("mu", "sig", "S1", "phi"))
## To check the fit (see plot #0 below) ##
## Simulates samples of the same size as the observed from the parameter values of the posterior
## A list with one simulated sample sad for each parameter combination from the posterior distribution
indices <- 1:nrow(fit1.df)
lista.fit <- vector(mode="list", length=length(indices))
for(i in 1:length(lista.fit))
    lista.fit[[i]] <- sim.nbln(indices[i], df=fit1.df, N=fit1.df$S2[indices[i]], a=sum(Y)/2.27e11)

## Species richness estimators for the current and original whole community x sampled sads ##

## A conservative estimator of species richness:
## Uses the mean of the lognormal sad behind the BN-lognormal.
## This mean is the mean abundance/species, which allows to estimate the total number
## of species if you know the total number of individuals (227 millions in the original Mata Atlantica)
## Credibility interval of this estimate
quantile(fit1.df$S1, na.rm=TRUE, c(0.025, .975)) ## problems: wide interval that includes values smaller thab the known number of species



####################################################################################################
## Old codes from mata Atlantica analyses, to be checked and adapt for Amazon data
####################################################################################################
## A second estimator of total species richness, using simulated samples with the estimated parameters fo the NB-Lognormal
## To simulate samples from the whole system and estimate the number os species with zero abundance in the sample
## First step : Simulates the samples and tally the proportion of species with abundance = zero
fit1.df$P0 <- sapply(1:nrow(fit1.df), P0, df=fit1.df, N=1e5, a=m1.data$a) ## Estimates the Proportion of species with abundance =0 from the NB-lognormal
## Estimates the total richness as the observed richness times 1+proportion of species with zero abundance in the sample
fit1.df$S2 <- length(Y)/(1-fit1.df$P0) ## Estimated richness = number of  known spp + estimated proportion of unrecorded species
## Empirical credibility interval
quantile(fit1.df$S2, na.rm=TRUE, c(0.025, .975)) ## More reasonable and precise than previous estimator


## The following code estimates the posterior mean and credibility intervals
## for the number of species in each octave in the sample, and in the whole community (current and original area of
## the Mata Atlantica)
## A list with one simulated regional sad (42.6 million trees) for each parameter combination from the posterior distribution
indices <- 1:nrow(fit1.df)
lista <- vector(mode="list", length=length(indices))
for(i in 1:length(lista))
    lista[[i]] <- sim.nbln(indices[i], df=fit1.df, N=fit1.df$S2[indices[i]])
## Octav plots of relative abundances
## range of octaves to be used in all octav tables
oc1 <- -36:0
## List with the octave tables from the simulated samples of the size as the observed sample
oc2 <- 0:20
df.fit2 <- sapply(lista.fit, function(x) octav(x, oct=oc2)$Freq)
## Predicted octavplot for the current area area of Amazon
Mean.fit <- apply(df.fit2,1,mean) # Expected number of especies in each octave
Low.fit <- apply(df.fit2,1,quantile,0.025)
Up.fit <- apply(df.fit2,1,quantile,0.975)

## List with the octave tables from the simulated samples from the current area of AF (42.6 million trees)
lista2 <- lapply(lista, function(x) octav(x/sum(x), oct=oc1))
df2 <- sapply(lista2, function(x) x$Freq)
## Predicted octavplot for the current area of the Mata Atlantica (42.6 million of trees)
Mean <- apply(df2,1,mean) # Expected number of especies in each octave
Low <- apply(df2,1,quantile,0.025)
Up <- apply(df2,1,quantile,0.975)

## Predicted Octavplot for the original area of the Mata Atlantica, from the underlying lognormal distribution
## with the estimated number of species
df3 <- sapply(1:nrow(fit1.df), op.lnorm, df=fit1.df, oct=oc1)
Mean3 <- apply(df3,1,mean) # Expected number of especies in each octave
Low3 <- apply(df3,1,quantile,0.025)
Up3 <- apply(df3,1,quantile,0.975)

## Estimates do evaluate Hyperdominance ##
## Estimates of hyperdominance: number of species that accumulate 75% and 95% of the total number of individuals
## In the current area of Mata Atlantica (which is taken as binomial negative sample of a lognormal community)
hypd1 <- sapply(lista, cum.quant, c(0.75, 0.95))
## Mean values
apply(hypd1,1,mean)## about 400 spp accumulates 75% of all individuals and about 1200 accumulates 95% of the indivuals
## Credibility intervals
apply(hypd1,1,quantile, c(0.025,0.975))
## In the sample
cum.quant(obs.sad, c(0.75, 0.95)) ## Very close to the estimated for the remaining part of the MA.
## In the original area (the lognormal sad assumed for the whole area)
## Simulates 1000 sads from the lognormal with the parameters taken from the posterior
lista4 <- lapply(1:nrow(fit1.df), function(i) with(fit1.df, rlnorm(n=ceiling(S2[i]), meanlog=mu[i], sdlog=sig[i])))
hypd2 <- sapply(lista4, cum.quant, c(0.75, 0.95))
## Mean values: in the hypothetic whole communtiy the dominance was about one third. 
apply(hypd2,1,mean)
## Credibility intervals
apply(hypd2,1,quantile, c(0.025,0.975))


## The plots ##

## 0. Observed x fitted sad sampled: maybe suppl material to show that the fit is nice and better than zero-truncated poilog
plot(octav(obs.sad, oct=oc2), ylim=c(0,400))
points(oc2-0.5, Mean.fit, col="blue", pch=19, type="b")
segments(oc2-0.5, Low.fit, oc2-0.5, Up.fit, col="blue")
## To compate with the Fit to zero-truncated poilog
Y.plt <- fitpoilog(obs.sad)
lines(octavpred(Y.plt), col="red")

## 1. observed x estimated sad for current area
## Grey bars: observed sad from the pooled plots, relative abundances
## Blue: estimated sad for current area of the MA (NB-lognormal sample of 42.6/277 Million of trees)
plot(octav(obs.sad/sum(obs.sad), oct=oc1), ylim=c(0,400), x.oct=TRUE, xlab="Relative abundance class (Log2)")
lines(oc1-0.5, Mean, pch=19, col="blue", type="b")
segments(x0=oc1-0.5, y0=Low, x1=oc1-0.5, y1=Up, col="blue")
## Alternative commented: confidence intervals as a ribbon
## polygon(c(rev(oc1-0.5), oc1-0.5), c(rev(Low),Up), col=rgb(0,0,1, alpha=0.5), border=NA)
## lines(oc1-0.5, Mean, col="darkblue")

# 2. Adds estimated lognormal sad for original area of the Mata Atlantica (not sure if helps)
## Red: estimated  sad fro original area of MA (Lognormal simulated for 227 million trees)
#pdf("nblnorm.pdf", width=9.5, height=6)
plot(octav(obs.sad/sum(obs.sad), oct=oc1), ylim=c(0,1250), x.oct=TRUE, xlab="Relative abundance class (Log2)")
lines(oc1-0.5, Mean, col="darkblue", type="b")
segments(x0=oc1-0.5, y0=Low, x1=oc1-0.5, y1=Up, col="blue")
lines(oc1-0.5, Mean3, col="red", type="b")
segments(x0=oc1-0.5, y0=Low3, x1=oc1-0.5, y1=Up3, col="red")
legend("topleft",
       legend=c(paste("Sample (",round(sum(Y)/1e6,1)," millions of trees)", sep=""),"Current (42.6 millions of trees)", "Original (227 millions of trees)"),
       col=c("grey","blue","red"), lty=c(NA,1,1), pch=c(15,1,2), bty="n")
#dev.off()
## polygon(c(rev(oc1-0.5), oc1-0.5), c(rev(Low),Up), col=rgb(0,0,1, alpha=0.5), border=NA)
## lines(oc1-0.5, Mean, col="darkblue")
## polygon(c(rev(oc1-0.5), oc1-0.5), c(rev(Low3),Up3), col=rgb(1,0,0, alpha=0.5), border=NA)
## lines(oc1-0.5, Mean3, col="red")

## 3. Rad plots of proportion of relative abundances in the sample an in the community (current area of the Mata Atlantica)
## I think octavplots are doing better, but just in case ...
plot(rad(lista[[1]]), type="n", prop=TRUE, xlim=c(1,max(sapply(lista,length))), ylim=c(1e-11,0.1))
sapply(lista, function(x)lines(rad(x), col="grey", prop=TRUE))
lines(rad(obs.sad), prop=TRUE, lwd=2)
