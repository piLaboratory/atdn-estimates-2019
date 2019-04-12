source("../functions.R")
library(abc)
library(parallel)
load("../lists_with_all_objects.RData")


## Load simulation ressults (ran in a computer cluster)
load("abc_simulations/abcFinal2013tax2019.RData")
## Uses only the summary statistics of the simulations with noise in estimated total population sizes (see abc2013trun.R)
abc2013t$sims <- abc2013t$sims[,5:8]

## Model selection
## Target: observed number of species, lmean, sdmean and zero of Mean_square with obs values             
target <- c(atdn.13.tax$Sobs, D(atdn.13.tax$data$population),
            mean(log(atdn.13.tax$data$population)), sd(log(atdn.13.tax$data$population)))

## Quick diagnostics plots
## Box plots of each target variable
## 
par(mfrow=c(2,2))
for(i in 1:ncol(abc2013t$sims)){
    boxplot(abc2013t$sims[,i]~abc2013t$labels, main=colnames(abc2013t$sims)[i], log="y")
    abline(h=target[i], lty=2, col="blue")
    }
par(mfrow=c(1,1))

## S in the sample x S total for each model
par(mfrow=c(2,2))
cores <- c(LSrnd="black", LSclump="blue", NBrnd="red", NBclump="green", LNrnd="orange", LNclump="grey")
plot(abc2013t$sims[,"S2"]~abc2013t$params, type="n")
for(n in unique(abc2013t$labels)){
    points(abc2013t$sims[abc2013t$labels==n,"S2"]~abc2013t$params[abc2013t$labels==n], col=cores[n])
}
abline(h=nrow(atdn.13.tax$data), lty=2)
## logmean x S
plot(abc2013t$sims[,"lmean2"]~abc2013t$params, type="n")
for(n in unique(abc2013t$labels)){
    points(abc2013t$sims[abc2013t$labels==n,"lmean2"]~abc2013t$params[abc2013t$labels==n], col=cores[n])
}
legend("topright", unique(abc2013t$labels), pch=1,
       col=c("black", "blue", "red", "green", "orange", "grey"), bty="n")
abline(h=mean(log(atdn.13.tax$data$population)), lty=2)
## logsd x S
plot(abc2013t$sims[,"lsd2"]~abc2013t$params, type="n")
for(n in unique(abc2013t$labels)){
    points(abc2013t$sims[abc2013t$labels==n,"lsd2"]~abc2013t$params[abc2013t$labels==n], col=cores[n])
}
abline(h=sd(log(atdn.13.tax$data$population)), lty=2)
## D x S
plot(abc2013t$sims[,"D2"]~abc2013t$params, type="n")
for(n in unique(abc2013t$labels)){
    points(abc2013t$sims[abc2013t$labels==n,"D2"]~abc2013t$params[abc2013t$labels==n], col=cores[n])
}
abline(h=D(atdn.13.tax$data$population), lty=2)
par(mfrow=c(1,1))

## Simulated S x logmean
par(mfrow=c(2,2))
plot(lmean2~S2, data=abc2013t$sims, type="n")
for(n in unique(abc2013t$labels))
    points(lmean2~S2, data=abc2013t$sims[abc2013t$labels==n,], cex=0.5, col=cores[n])
points(nrow(atdn.13.tax$data), mean(log(atdn.13.tax$data$population)), pch=13, cex=2, col="orange")
legend("topright", unique(abc2013t$labels), pch=1,
       col=c("black", "blue", "red", "green", "orange", "grey"), bty="n")
## Simulated S x logsd
plot(lsd2~S2, data=abc2013t$sims, type="n")
for(n in unique(abc2013t$labels))
    points(lsd2~S2, data=abc2013t$sims[abc2013t$labels==n,], cex=0.5, col=cores[n])
points(nrow(atdn.13.tax$data), sd(log(atdn.13.tax$data$population)), pch=13, cex=2, col="orange")
## lmean x lsd2
plot(lsd2~lmean2, data=abc2013t$sims, type="n")
for(n in unique(abc2013t$labels))
    points(lsd2~lmean2, data=abc2013t$sims[abc2013t$labels==n,], cex=0.5, col=cores[n])
points(mean(log(atdn.13.tax$data$population)), sd(log(atdn.13.tax$data$population)), pch=13, cex=2, col="orange")
## S x D
plot(D2~S2, data=abc2013t$sims, type="n", ylim=range(c(abc2013t$sims$D2, D(atdn.13.tax$data$population))))
for(n in unique(abc2013t$labels))
    points(D2~S2, data=abc2013t$sims[abc2013t$labels==n,], cex=0.5, col=cores[n])
points(atdn.13.tax$Sobs, D(atdn.13.tax$data$population), pch=13, cex=2, col="orange")
par(mfrow=c(1,1))

## Model selection
## Cross-validation
cv.modsel <- cv4postpr(abc2013t$labels, abc2013t$sims, nval=100, tol=0.025, method="rejection")
summary(cv.modsel)
## Model selection
model.sel <- postpr(target = target,
                    index=abc2013t$labels,
                    sumstat = abc2013t$sims,
                    tol=0.025, method="rejection",
                    corr=TRUE)
summary(model.sel)
## Plausible models
index <- abc2013t$labels=="NBclump"|abc2013t$labels=="NBrnd"|abc2013t$labels=="LSrnd"|abc2013t$labels=="LSclump"

## Goodness of fit the models
summary(
    gfit(target = target,
                    sumstat = abc2013t$sims[index,],
                    nb.replicate = 200, tol = 0.05)
    )


## Parameter estimation ##
## Cross validation
cv.nn <- cv4abc(param=abc2013t$params[index],
       sumstat = abc2013t$sims[index,],
       tols=c(0.05, 0.025, 0.01), nval=30, method="neuralnet")
cv.rej <- cv4abc(param=abc2013t$params[index],
       sumstat = abc2013t$sims[index,],
       tols=c(0.05, 0.025, 0.01), nval=30, method="rejection")
cv.rr <- cv4abc(param=data.frame(S=abc2013t$params[index]),
       sumstat = abc2013t$sims[index,],
       tols=c(0.05, 0.025, 0.01), nval=30, method="ridge")
par(mfrow=c(1,3))
plot(cv.nn, caption="Neural net")
plot(cv.rej, caption="Rejection")
plot(cv.rr, caption="Ridge regression")
par(mfrow=c(1,1))

## Posterior distribution of Species richness from the selected model

t1 <- 0.025
S.post1 <- abc(target = target, param=data.frame(S=abc2013t$params[index]),
              sumstat = abc2013t$sims[index,],
              tol=t1, method="rejection")
S.post2 <- abc(target = target, param=data.frame(S=abc2013t$params[index]),
              sumstat = abc2013t$sims[index,],
              tol=t1, method="neuralnet", numnet=100)
S.post3 <- abc(target = target[c(1,3)], param=data.frame(S=abc2013t$params[index]),
              sumstat = abc2013t$sims[index,c(1,3)],tol=t1, method="ridge")
summary(S.post1)
summary(S.post2)
summary(S.post3)
par(mfrow=c(1,3))
hist(S.post1)
hist(S.post2)
hist(S.post3)
par(mfrow=c(1,1))

## Posterior predictive check
## Function to parallelize simulations, that are used in the checks
f1 <- function(x, ...){
    y <- with(atdn.13.tax,
              try(sim.abc(S = x,
                          N = Tot.t,
                          n.plots = N.plots,
                          tot.area= Tot.A,
                          nb.fit = y.nb2,
                          lm.sd.fit = lm.sd,
                          lmk.fit = lm.k,
                          nrep = 1, ...))
              )
    if(class(y)=="try-error")
        return(matrix(NA, nrow=2, ncol=8))
    else
        return(y)
}

### Simulations and predcitive check
LSclump.pc <- mclapply(sample(S.post3$unadj.values,100, replace=TRUE),
                  f1, sad="ls", lower=1e-20, upper=1e20, mc.cores=3)
tmp1 <- LSclump.pc[[1]][2,5:8]
for(i in 2:length(LSclump.pc))
    tmp1 <- rbind(tmp1, LSclump.pc[[i]][2,5:8])
par(mfrow=c(2,2))
for(i in 1:4){
    hist(tmp1[,i], xlim=range(c(tmp1[,i],target[i])), main=names(tmp1[i]))
    abline(v=target[i], col="red")
}
par(mfrow=c(1,1))
         
