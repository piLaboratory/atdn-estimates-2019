source("../functions.R")
library(abc)
library(parallel)
load("../lists_with_all_objects.RData")


## Load simulation ressults (ran in a computer cluster)
load("abc_simulations/abcFinal2013.RData")
## Uses only the summary statistics of the simulations with noise in estimated total population sizes (see abc2013run.R)
abc2013$sims <- abc2013$sims[,5:8]

## Model selection
## Target: observed number of species, Simpson's species-equivalent, lmean and sd mean of estimated population sizes             
target <- c(atdn.13$Sobs, D(atdn.13$data$population), mean(log(atdn.13$data$population)), sd(log(atdn.13$data$population)))

## Quick diagnostics plots
## Box plots of each parget variable
## 
par(mfrow=c(2,2))
for(i in 1:ncol(abc2013$sims)){
    boxplot(abc2013$sims[,i]~abc2013$labels, main=colnames(abc2013$sims)[i], log="y")
    abline(h=target[i], lty=2, col="blue")
    }
par(mfrow=c(1,1))

## S in the sample x S total for each model
par(mfrow=c(2,2))
cores <- c(LSrnd="black", LSclump="blue", NBrnd="red", NBclump="green", LNrnd="orange", LNclump="grey")
plot(abc2013$sims[,"S2"]~abc2013$params, type="n")
for(n in unique(abc2013$labels)){
    points(abc2013$sims[abc2013$labels==n,"S2"]~abc2013$params[abc2013$labels==n], col=cores[n])
}
abline(h=nrow(atdn.13$data), lty=2)
## logmean x S
plot(abc2013$sims[,"lmean2"]~abc2013$params, type="n")
for(n in unique(abc2013$labels)){
    points(abc2013$sims[abc2013$labels==n,"lmean2"]~abc2013$params[abc2013$labels==n], col=cores[n])
}
legend("topright", unique(abc2013$labels), pch=1,
       col=c("black", "blue", "red", "green", "orange", "grey"), bty="n")
abline(h=mean(log(atdn.13$data$population)), lty=2)
## logsd x S
plot(abc2013$sims[,"lsd2"]~abc2013$params, type="n")
for(n in unique(abc2013$labels)){
    points(abc2013$sims[abc2013$labels==n,"lsd2"]~abc2013$params[abc2013$labels==n], col=cores[n])
}
abline(h=sd(log(atdn.13$data$population)), lty=2)
## D x S
plot(abc2013$sims[,"D2"]~abc2013$params, type="n")
for(n in unique(abc2013$labels)){
    points(abc2013$sims[abc2013$labels==n,"D2"]~abc2013$params[abc2013$labels==n], col=cores[n])
}
abline(h=D(atdn.13$data$population), lty=2)
par(mfrow=c(1,1))

## Simulated S x logmean
par(mfrow=c(2,2))
plot(lmean2~S2, data=abc2013$sims, type="n")
for(n in unique(abc2013$labels))
    points(lmean2~S2, data=abc2013$sims[abc2013$labels==n,], cex=0.5, col=cores[n])
points(nrow(atdn.13$data), mean(log(atdn.13$data$population)), pch=13, cex=2, col="orange")
legend("topright", unique(abc2013$labels), pch=1,
       col=c("black", "blue", "red", "green", "orange", "grey"), bty="n")
## Simulated S x logsd
plot(lsd2~S2, data=abc2013$sims, type="n")
for(n in unique(abc2013$labels))
    points(lsd2~S2, data=abc2013$sims[abc2013$labels==n,], cex=0.5, col=cores[n])
points(nrow(atdn.13$data), sd(log(atdn.13$data$population)), pch=13, cex=2, col="orange")
## lmean x lsd2
plot(lsd2~lmean2, data=abc2013$sims, type="n")
for(n in unique(abc2013$labels))
    points(lsd2~lmean2, data=abc2013$sims[abc2013$labels==n,], cex=0.5, col=cores[n])
points(mean(log(atdn.13$data$population)), sd(log(atdn.13$data$population)), pch=13, cex=2, col="orange")
## S x D
plot(D2~S2, data=abc2013$sims, type="n", ylim=range(c(abc2013$sims$D2, D(atdn.13$data$population))))
for(n in unique(abc2013$labels))
    points(D2~S2, data=abc2013$sims[abc2013$labels==n,], cex=0.5, col=cores[n])
points(atdn.13$Sobs, D(atdn.13$data$population), pch=13, cex=2, col="orange")
par(mfrow=c(1,1))

## Model selection
## Cross-validation
cv.modsel <- cv4postpr(abc2013$labels, abc2013$sims, nval=100, tol=0.01, method="rejection")
summary(cv.modsel)
## Model selection
model.sel <- postpr(target = target,
                    index=abc2013$labels,
                    sumstat = abc2013$sims,
                    tol=0.025, method="rejection",
                    corr=TRUE)
summary(model.sel)

## Goodness of fit the models
summary(
    gfit(target = target[c(1,3)],
                    sumstat = abc2013$sims[abc2013$labels=="LSclump",c(1,3)],
                    nb.replicate = 200, tol = 0.05)
)

summary(
    gfit(target = target,
                    sumstat = abc2013$sims[abc2013$labels=="LSclump",],
                    nb.replicate = 200, tol = 0.05)
    )


## Parameter estimation ##
## Cross validation
cv.nn <- cv4abc(param=abc2013$params[abc2013$labels=="LSclump"],
       sumstat = abc2013$sims[abc2013$labels=="LSclump",],
       tols=c(0.05, 0.025, 0.01), nval=30, method="neuralnet")
cv.rej <- cv4abc(param=abc2013$params[abc2013$labels=="LSclump"],
       sumstat = abc2013$sims[abc2013$labels=="LSclump",],
       tols=c(0.05, 0.025, 0.01), nval=30, method="rejection")
cv.ll <- cv4abc(param=data.frame(S=abc2013$params[abc2013$labels=="LSclump"]),
       sumstat = abc2013$sims[abc2013$labels=="LSclump",],
       tols=c(0.05, 0.025, 0.01), nval=30, method="ridge")
par(mfrow=c(1,3))
plot(cv.nn, caption="Neural net")
plot(cv.rej, caption="Rejection")
plot(cv.ll, caption="Ridge regression")
par(mfrow=c(1,1))

## Posterior distribution of Species richness from the selected model
t1 <- 0.025
S.post1 <- abc(target = target, param=data.frame(S=abc2013$params[abc2013$labels=="LSclump"]),
              sumstat = abc2013$sims[abc2013$labels=="LSclump",],
              tol=t1, method="rejection")
S.post2 <- abc(target = target, param=data.frame(S=abc2013$params[abc2013$labels=="LSclump"]),
              sumstat = abc2013$sims[abc2013$labels=="LSclump",],
              tol=t1, method="neuralnet", numnet = 100)
S.post3 <- abc(target = target, param=data.frame(S=abc2013$params[abc2013$labels=="LSclump"]),
              sumstat = abc2013$sims[abc2013$labels=="LSclump",],tol=t1, method="ridge")
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
    y <- with(atdn.13,
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

### Simulations and predictive check
LSclump.pc <- mclapply(sample(S.post1$unadj.values,100, replace=TRUE),
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
         
