source("functions.R")
library(abc)
library(parallel)
## Basic quantities and original (biased estimates), from script 'dataprep.R'
load("lists_with_all_objects.RData")


################################################################################
## 2013 original data##
################################################################################
## results of simulations for ABC, from script 'simulations_abc/2013/abc_run2013[a,b].R'
## and then applying 'simulations_abc/2013/join_simulations.R'
load("abcFinal2013.RData") 
## Use only the summary statistics of the simulations with noise in
## estimated total population sizes (see abc2019run.R) 
abc2013$sims <- abc2013$sims[,5:8]

## Model selection
## Target: observed number of species, Simpson's species-equivalents,
## mean and sd of log abundances
target <- c( atdn.13$Sobs,
            D(atdn.13$data$population),
            mean(log(atdn.13$data$population)),
            sd(log(atdn.13$data$population))
            )

## Diagnostics plots
## Box plots of each target variable
par(mfrow=c(2,2))
for(i in 1:ncol(abc2013$sims)){
    boxplot(abc2013$sims[,i]~abc2013$labels,
            main=colnames(abc2013$sims)[i], log="y")
    abline(h=target[i], lty=2, col="blue")
    }
par(mfrow=c(1,1))


## Model selection ##

## Cross-validation
cv.modsel <- cv4postpr(abc2013$labels, abc2013$sims,
                       nval=100, tol= c(0.05, 0.025, 0.01),
                       method="rejection")
summary(cv.modsel)

## Model selection
model.sel <- postpr(target = target,
                    index=abc2013$labels,
                    sumstat = abc2013$sims,
                    tol=0.01, method="rejection",
                    corr=TRUE)
summary(model.sel)

## Selected model(s): check from the command above and write the selected model(s) in the following command

index <- abc2013$labels=="LSclump" # write here the code of the selected model(s)

## Goodness of fit of the selected model
nrep <- 200
gof <- gfit(target = target,
            sumstat = abc2013$sims[index,],
            nb.replicate = nrep, tol = 0.025)
sgof <- summary(gof)
p.legend <- ifelse(sgof$pvalue==0,
                   paste("p < 1/",nrep,sep=""),
                   paste("p =",sgof$pvalue))
plot(gof)
mtext(p.legend)


## Parameter estimation ##
## Cross validation
cv.rej <- cv4abc(param = data.frame(S=abc2013$params[index]),
       sumstat = abc2013$sims[index,],
       tols=c(0.05, 0.025, 0.01), nval=100, method="rejection")
cv.neural <- cv4abc(param = data.frame(S=abc2019$params[index]),
       sumstat = abc2019$sims[index,],
       tols=c(0.05, 0.025, 0.01), nval=100, method="neuralnet")
cv.ll <- cv4abc(param = data.frame(S=abc2019$params[index]),
       sumstat = abc2019$sims[index,],
       tols=c(0.05, 0.025, 0.01), nval=100, method="loclinear")
cv.ridge <- cv4abc(param = data.frame(S=abc2019$params[index]),
       sumstat = abc2019$sims[index,],
       tols=c(0.05, 0.025, 0.01), nval=100, method="ridge")
## Checks the method and tolerance with smaller error
summary(cv.rej)
summary(cv.neural)
summary(cv.ll)
summary(cv.ridge)
## Chosen loclinear as ridge with tol=0.01 has a slightly smaller error
## but returns posteriors with negative values
## Chosen tol=0.05 as the difference in loclinear among tolerance levels
## is less than 1e-5, and tol=0.05 enables a posterior drawn from more simulated values

## Posterior distribution of Species richness from the selected model
S.post1 <- abc(target = target, param=data.frame(S=abc2013$params[index]),
              sumstat = abc2013$sims[index,],
              tol=0.05, method="loclinear")
summary(S.post1)
hist(S.post1)

## stores all relevant objects in a list
abc2013.summ <- list(target=target,
                   cv.modsel=cv.modsel, model.sel=model.sel,
                   index=index, gof=gof,
                   sgof=sgof, p.legend=p.legend,
                   cv.rej=cv.rej, cv.neural=cv.neural, cv.ll=cv.ll, cv.ridge=cv.ridge,
                   S.post1=S.post1)


################################################################################
## 2013 revised taxonomy ##
################################################################################
load("abcFinal2013tax2019.RData")
## Use only the summary statistics of the simulations with noise in
## estimated total population sizes (see abc2019run.R) 
abc2013t$sims <- abc2013t$sims[,5:8]

## Model selection
## Target: observed number of species, Simpson's species equivalent,
## mean and sd of log abundances

target <- c( atdn.13.tax$Sobs,
            D(atdn.13.tax$data$population),
            mean(log(atdn.13.tax$data$population)),
            sd(log(atdn.13.tax$data$population))
            )

## Diagnostics plots
## Box plots of each target variable
par(mfrow=c(2,2))
for(i in 1:ncol(abc2013t$sims)){
    boxplot(abc2013t$sims[,i]~abc2013t$labels,
            main=colnames(abc2013t$sims)[i], log="y")
    abline(h=target[i], lty=2, col="blue")
    }
par(mfrow=c(1,1))


## Model selection ##

## Cross-validation
cv.modsel <- cv4postpr(abc2013t$labels, abc2013t$sims,
                       nval=100, tol= c(0.05, 0.025, 0.01),
                       method="rejection")
summary(cv.modsel)

## Model selection
model.sel <- postpr(target = target,
                    index=abc2013t$labels,
                    sumstat = abc2013t$sims,
                    tol=0.01, method="rejection",
                    corr=TRUE)
summary(model.sel)
## Selected model(s)
index <- abc2013t$labels=="LSclump"


## Goodness of fit the models
nrep <- 200
gof <- gfit(target = target,
            sumstat = abc2013t$sims[index,],
            nb.replicate = nrep, tol = 0.025)
sgof <- summary(gof)
p.legend <- ifelse(sgof$pvalue==0,
                   paste("p < 1/",nrep,sep=""),
                   paste("p =",sgof$pvalue))
plot(gof)
mtext(p.legend)


## Parameter estimation ##
## Cross validation
cv.rej <- cv4abc(param = data.frame(S=abc2013t$params[index]),
       sumstat = abc2013t$sims[index,],
       tols=c(0.05, 0.025, 0.01), nval=100, method="rejection")
cv.neural <- cv4abc(param = data.frame(S=abc2019$params[index]),
       sumstat = abc2019$sims[index,],
       tols=c(0.05, 0.025, 0.01), nval=100, method="neuralnet")
cv.ll <- cv4abc(param = data.frame(S=abc2019$params[index]),
       sumstat = abc2019$sims[index,],
       tols=c(0.05, 0.025, 0.01), nval=100, method="loclinear")
cv.ridge <- cv4abc(param = data.frame(S=abc2019$params[index]),
       sumstat = abc2019$sims[index,],
       tols=c(0.05, 0.025, 0.01), nval=100, method="ridge")
## Checks the method and tolerance with smaller error
summary(cv.rej)
summary(cv.neural)
summary(cv.ll)
summary(cv.ridge)
## Loc-linear with tol=0.01

## Posterior distribution of Species richness from the selected model
S.post1 <- abc(target = target, param=data.frame(S=abc2013t$params[index]),
              sumstat = abc2013t$sims[index,],
              tol=0.01, method="loclinear")

summary(S.post1)
hist(S.post1)

## stores all relevant objects in a list
abc2013t.summ <- list(target=target,
                   cv.modsel=cv.modsel, model.sel=model.sel,
                   index=index, gof=gof,
                   sgof=sgof, p.legend=p.legend,
                   cv.rej=cv.rej, cv.neural=cv.neural, cv.ll=cv.ll, cv.ridge=cv.ridge,
                   S.post1=S.post1)

################################################################################
## 2019 ##
################################################################################
load("abcFinal2019.RData")
## Use only the summary statistics of the simulations with noise in
## estimated total population sizes (see abc2019run.R) 
abc2019$sims <- abc2019$sims[,5:8]

## Model selection
## Target: observed number of species, lmean, sdmean of log abundances 
target <- c( atdn.19$Sobs,
           D(atdn.19$data$population), 
           mean(log(atdn.19$data$population)),
           sd(log(atdn.19$data$population))
           )

## Quick diagnostics plots
## Box plots of each target variable
## 
par(mfrow=c(2,2))
for(i in 1:ncol(abc2019$sims)){
    boxplot(abc2019$sims[,i]~abc2019$labels,
            main=colnames(abc2019$sims)[i], log="y")
    abline(h=target[i], lty=2, col="blue")
    }
par(mfrow=c(1,1))


## Model selection ##

## Cross-validation
cv.modsel <- cv4postpr(abc2019$labels, abc2019$sims,
                       nval=100, tol= c(0.05, 0.025, 0.01),
                       method="rejection")
summary(cv.modsel)

## Model selection
model.sel <- postpr(target = target,
                    index=abc2019$labels,
                    sumstat = abc2019$sims,
                    tol=0.01, method="rejection",
                    corr=TRUE)
summary(model.sel)
## Selected model(s)
index <- abc2019$labels=="LSclump"

## Goodness of fit the models
nrep <- 200
gof <- gfit(target = target,
            sumstat = abc2019$sims[index,],
            nb.replicate = nrep, tol = 0.01)
sgof <- summary(gof)
p.legend <- ifelse(sgof$pvalue==0,
                   paste("p < 1/",nrep,sep=""),
                   paste("p =",sgof$pvalue))
plot(gof)
mtext(p.legend)


## Parameter estimation ##
## Cross validation
cv.rej <- cv4abc(param = data.frame(S=abc2019$params[index]),
       sumstat = abc2019$sims[index,],
       tols=c(0.05, 0.025, 0.01), nval=100, method="rejection")
cv.neural <- cv4abc(param = data.frame(S=abc2019$params[index]),
       sumstat = abc2019$sims[index,],
       tols=c(0.05, 0.025, 0.01), nval=100, method="neuralnet")
cv.ll <- cv4abc(param = data.frame(S=abc2019$params[index]),
       sumstat = abc2019$sims[index,],
       tols=c(0.05, 0.025, 0.01), nval=100, method="loclinear")
cv.ridge <- cv4abc(param = data.frame(S=abc2019$params[index]),
       sumstat = abc2019$sims[index,],
       tols=c(0.05, 0.025, 0.01), nval=100, method="ridge")
## Loclinear with tol=0.05

## Posterior distribution of Species richness from the selected model
S.post1 <- abc(target = target, param=data.frame(S=abc2019$params[index]),
              sumstat = abc2019$sims[index,],
              tol=0.05, method="loclinear", numnet=100)

summary(S.post1)
hist(S.post1)


## stores all relevant objects in a list
abc2019.summ <- list(target=target,
                   cv.modsel=cv.modsel, model.sel=model.sel,
                   index=index, gof=gof,
                   sgof=sgof, p.legend=p.legend,
                   cv.rej=cv.rej, cv.neural=cv.neural, cv.ll=cv.ll, cv.ridge=cv.ridge,
                   S.post1=S.post1)

save(abc2013.summ, abc2013t.summ, abc2019.summ, file="abcSummaries.RData")
