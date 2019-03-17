source("functions.R")
library(abc)


## Simulations
## Function to parallelize simulations
f1 <- function(x, ...){
    y <- try(sim.abc(S = x, N = Tot.t, n.plots = N.plots, tot.area= Tot.A,
                     nb.fit=y.nb2, lmk.fit = lm.k,
                     nrep = 1, ...))
    if(class(y)=="try-error")
        return(matrix(NA, nrow=2, ncol=4))
    else
        return(y)
}

## Load simulation ressults (ran in a computer cluster)
load("abc_simulations/abc2019.RData")

## Model selection
## Target: observed number of species, lmean, sdmean and zero of Mean_square with obs values             
target <- c(Sobs, D(atdn.2019$population), mean(log(atdn.2019$population)), sd(log(atdn.2019$population)))

## Quick diagnostics plots
## Box plots of each parget variable
## 
par(mfrow=c(2,2))
for(i in 1:ncol(abc2019$sims)){
    boxplot(abc2019$sims[,i]~abc2019$labels, main=colnames(abc2019$sims)[i], log="y")
    abline(h=target[i], lty=2, col="blue")
    }
par(mfrow=c(1,1))

## S in the sample x S total for each model
par(mfrow=c(2,2))
cores <- c(LSrnd="black", LSclump="blue", NBrnd="red", NBclump="green", LNrnd="orange", LNclump="grey")
plot(abc2019$sims[,"S"]~abc2019$params, type="n")
for(n in unique(abc2019$labels)){
    points(abc2019$sims[abc2019$labels==n,"S"]~abc2019$params[abc2019$labels==n], col=cores[n])
}
abline(h=nrow(atdn.2019), lty=2)
## logmean x S
plot(abc2019$sims[,"lmean"]~abc2019$params, type="n")
for(n in unique(abc2019$labels)){
    points(abc2019$sims[abc2019$labels==n,"lmean"]~abc2019$params[abc2019$labels==n], col=cores[n])
}
legend("topright", unique(abc2019$labels), pch=1,
       col=c("black", "blue", "red", "green", "orange", "grey"), bty="n")
abline(h=mean(log(atdn.2019$population)), lty=2)
## logsd x S
plot(abc2019$sims[,"lsd"]~abc2019$params, type="n")
for(n in unique(abc2019$labels)){
    points(abc2019$sims[abc2019$labels==n,"lsd"]~abc2019$params[abc2019$labels==n], col=cores[n])
}
abline(h=sd(log(atdn.2019$population)), lty=2)
## D x S
plot(abc2019$sims[,"D"]~abc2019$params, type="n")
for(n in unique(abc2019$labels)){
    points(abc2019$sims[abc2019$labels==n,"D"]~abc2019$params[abc2019$labels==n], col=cores[n])
}
abline(h=D(atdn.2019$population), lty=2)
par(mfrow=c(1,1))

## Simulated S x logmean
par(mfrow=c(2,2))
plot(lmean~S, data=abc2019$sims, type="n")
for(n in unique(abc2019$labels))
    points(lmean~S, data=abc2019$sims[abc2019$labels==n,], cex=0.5, col=cores[n])
points(nrow(atdn.2019), mean(log(atdn.2019$population)), pch=19, cex=2, col="orange")
legend("topright", unique(abc2019$labels), pch=1,
       col=c("black", "blue", "red", "green", "orange", "grey"), bty="n")
## Simulated S x logsd
plot(lsd~S, data=abc2019$sims, type="n")
for(n in unique(abc2019$labels))
    points(lsd~S, data=abc2019$sims[abc2019$labels==n,], cex=0.5, col=cores[n])
points(nrow(atdn.2019), sd(log(atdn.2019$population)), pch=19, cex=2, col="orange")
## lmean x lsd
plot(lsd~lmean, data=abc2019$sims, type="n")
for(n in unique(abc2019$labels))
    points(lsd~lmean, data=abc2019$sims[abc2019$labels==n,], cex=0.5, col=cores[n])
points(mean(log(atdn.2019$population)), sd(log(atdn.2019$population)), pch=19, cex=2, col="orange")
## S x D
plot(D~S, data=abc2019$sims, type="n", ylim=range(c(abc2019$sims$D, D(atdn.2019$population))))
for(n in unique(abc2019$labels))
    points(D~S, data=abc2019$sims[abc2019$labels==n,], cex=0.5, col=cores[n])
points(Sobs, D(atdn.2019$population), pch=19, cex=2, col="orange")
par(mfrow=c(1,1))

## Model selection
## Cross-validation
cv.modsel <- cv4postpr(abc2019$labels, abc2019$sims, nval=100, tol=0.01, method="rejection")
summary(cv.modsel)
## Model selection
model.sel <- postpr(target = target,
                    index=abc2019$labels,
                    sumstat = abc2019$sims,
                    tol=0.025, method="rejection",
                    corr=TRUE)
summary(model.sel)

## Goodness of fit the models
summary(
    gfit(target = target[c(1,3)],
                    sumstat = abc2019$sims[abc2019$labels=="LSclump",c(1,3)],
                    nb.replicate = 200, tol = 0.05)
)

summary(
    gfit(target = target[c(1,2,4)],
                    sumstat = abc2019$sims[abc2019$labels=="LSclump",c(1,2,4)],
                    nb.replicate = 200, tol = 0.05)
    )


## Parameter estimation ##
## Cross validation
cv.nn <- cv4abc(param=abc2019$params[abc2019$labels=="LSclump"],
       sumstat = abc2019$sims[abc2019$labels=="LSclump",],
       tols=c(0.05, 0.025, 0.01), nval=30, method="neuralnet")
cv.rej <- cv4abc(param=abc2019$params[abc2019$labels=="LSclump"],
       sumstat = abc2019$sims[abc2019$labels=="LSclump",],
       tols=c(0.05, 0.025, 0.01), nval=30, method="rejection")
cv.ll <- cv4abc(param=data.frame(S=abc2019$params[abc2019$labels=="LSclump"]),
       sumstat = abc2019$sims[abc2019$labels=="LSclump",],
       tols=c(0.05, 0.025, 0.01), nval=30, method="ridge")
par(mfrow=c(1,3))
plot(cv.nn, caption="Neural net")
plot(cv.rej, caption="Rejection")
plot(cv.ll, caption="Ridge regression")
par(mfrow=c(1,1))

## Posterior distribution of Species richness from the selected model
t1 <- 0.05
S.post1 <- abc(target = target, param=data.frame(S=abc2019$params[abc2019$labels=="LSclump"]),
              sumstat = abc2019$sims[abc2019$labels=="LSclump",],
              tol=t1, method="rejection")
S.post2 <- abc(target = target, param=data.frame(S=abc2019$params[abc2019$labels=="LSclump"]),
              sumstat = abc2019$sims[abc2019$labels=="LSclump",],
              tol=t1, method="neuralnet")
S.post3 <- abc(target = target, param=data.frame(S=abc2019$params[abc2019$labels=="LSclump"]),
              sumstat = abc2019$sims[abc2019$labels=="LSclump",],
              tol=t1, method="loclinear")
summary(S.post1)
summary(S.post2)
summary(S.post3)
par(mfrow=c(1,3))
hist(S.post1)
hist(S.post2)
hist(S.post3)

## Posterior predictive check
LSclump.pc <- mclapply(sample(S.post1$unadj.values,100, replace=TRUE),
                  f1, sad="ls", lower=1e-20, upper=1e20, mc.cores=4)
tmp1 <- LSclump.pc[[1]][2,]
for(i in 2:length(LSclump.pc))
    tmp1 <- rbind(tmp1, LSclump.pc[[i]][2,])
par(mfrow=c(2,2))
for(i in 1:4){
    hist(tmp1[,i], xlim=range(c(tmp1[,i],target[i])), main=names(LSclump.pc[[1]])[i])
    abline(v=target[i], col="red")
}
par(mfrow=c(1,1))
         
