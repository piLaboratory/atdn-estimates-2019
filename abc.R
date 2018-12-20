source("functions.R")
library(abc)
library(parallel)
library(untb)
library(sads)

## Simulations
## Function to parallelize simulations
f1 <- function(x, ...){
    y <- try(sim.abc(S = x, N = Tot.t, n.plots = N.plots, tot.area= Tot.A,
                     nb.fit=y.nb2, lmk = lm.k, obs.values=Dec2018$population,
                     nrep = 1, ...))
    if(class(y)=="try-error")
        return(matrix(NA, nrow=2, ncol=4))
    else
        return(y)
}


## Number of Simulated values
nsims <- 6135
##simulated.vals <- runif(nsims, 1e4, S.orc$S.est)
simulated.vals <- runif(nsims, 1e4, 2e4)
## Runs the simulations with mclapply (a bit slow, ca 400 simulations/hour with 4 cores)
## Random sample of LS and TNB
LS.sims <- mclapply(simulated.vals, f1, mc.cores=4, lower=1e-20, upper=1e20)
save.image()
## Clumped sample of LS and TNB
NB.sims <- mclapply(simulated.vals, f1, lower=1e-20, upper=1e20, LS = FALSE, mc.cores=4)
save.image()

## Assembles all simulation results in a matrix
## Simulations from LS rads ##
## exclude failed simulations with random sampling
LS.index <- sapply(LS.sims, function(x) !any(is.na(x)))
j1 <- (1:length(LS.sims))[LS.index]
## Assembles a matrix, each line the result of a simulation
all.sims <- LS.sims[[min(j1)]]
for(i in j1[-min(j1)])
    all.sims <- rbind(all.sims, LS.sims[[i]])
## Simulations from TNB rad
## Excluding failed simulations
NB.index <- sapply(NB.sims, function(x) !any(is.na(x)))
j2 <- (1:length(NB.sims))[NB.index]
## rbind the results of the simulations to the same results matrix
for(i in j2)
    all.sims <- rbind(all.sims, NB.sims[[i]])
## Vector with labels for each simulation
sim.ids <- c(rep(c("LSrnd", "LSclump"), sum(LS.index)),
             rep(c("NBrnd", "NBclump"), sum(NB.index)))## Simulated values for each simulation
## Vector with the parameters used in each simulation
## (which is the total species richness in the regional RADs)
sim.y <- c( rep(simulated.vals[j1],each=2),
           rep(simulated.vals[j2],each=2))             
## ABC ##
## Model selection
## Target: observed number of species, lmean, sdmean and zero of Mean_square with obs values             
target <- c(Sobs, D(Dec2018$population), mean(log(Dec2018$population)), sd(log(Dec2018$population)))

## Quick diagnostics plots
## Box plots of each parget variable
## 
par(mfrow=c(2,2))
for(i in 1:ncol(all.sims)){
    boxplot(all.sims[,i]~sim.ids, main=colnames(all.sims)[i], log="y")
    abline(h=target[i], lty=2, col="blue")
    }
par(mfrow=c(1,1))

## S in the sample x S total for each model
par(mfrow=c(2,2))
cores <- c(LSrnd="black", LSclump="blue", NBrnd="red", NBclump="green")
plot(all.sims[,"S"]~sim.y, type="n")
for(n in unique(sim.ids)){
    points(all.sims[sim.ids==n,"S"]~sim.y[sim.ids==n], col=cores[n])
}
abline(h=nrow(Dec2018), lty=2)
legend("topright", unique(sim.ids), pch=1,
       col=c("black", "blue", "red", "green"), bty="n")
## logmean x S
plot(all.sims[,"lmean"]~sim.y, type="n")
for(n in unique(sim.ids)){
    points(all.sims[sim.ids==n,"lmean"]~sim.y[sim.ids==n], col=cores[n])
}
abline(h=mean(log(Dec2018$population)), lty=2)
## logsd x S
plot(all.sims[,"lsd"]~sim.y, type="n")
for(n in unique(sim.ids)){
    points(all.sims[sim.ids==n,"lsd"]~sim.y[sim.ids==n], col=cores[n])
}
abline(h=sd(log(Dec2018$population)), lty=2)
## D x S
plot(all.sims[,"D"]~sim.y, type="n")
for(n in unique(sim.ids)){
    points(all.sims[sim.ids==n,"D"]~sim.y[sim.ids==n], col=cores[n])
}
abline(h=D(Dec2018$population), lty=2)
par(mfrow=c(1,1))

## Simulated S x logmean
par(mfrow=c(2,2))
plot(lmean~S, data=all.sims, type="n")
for(n in unique(sim.ids))
    points(lmean~S, data=all.sims[sim.ids==n,], cex=0.5, col=cores[n])
points(nrow(Dec2018), mean(log(Dec2018$population)), pch=19, cex=2, col="orange")
legend("topright", unique(sim.ids), pch=1,
       col=c("black", "blue", "red", "green"), bty="n")
## Simulated S x logsd
plot(lsd~S, data=all.sims, type="n")
for(n in unique(sim.ids))
    points(lsd~S, data=all.sims[sim.ids==n,], cex=0.5, col=cores[n])
points(nrow(Dec2018), sd(log(Dec2018$population)), pch=19, cex=2, col="orange")
## lmean x lsd
plot(lsd~lmean, data=all.sims, type="n")
for(n in unique(sim.ids))
    points(lsd~lmean, data=all.sims[sim.ids==n,], cex=0.5, col=cores[n])
points(mean(log(Dec2018$population)), sd(log(Dec2018$population)), pch=19, cex=2, col="orange")
## S x D
plot(D~S, data=all.sims, type="n", ylim=range(c(all.sims$D, D(Dec2018$population))))
for(n in unique(sim.ids))
    points(D~S, data=all.sims[sim.ids==n,], cex=0.5, col=cores[n])
points(Sobs, D(Dec2018$population), pch=19, cex=2, col="orange")
par(mfrow=c(1,1))

## Model selection
## Cross-validation
cv.modsel <- cv4postpr(sim.ids, all.sims, nval=100, tol=0.01, method="rejection")
summary(cv.modsel)
## Model selection
model.sel <- postpr(target = target,
                    index=sim.ids,
                    sumstat = all.sims,
                    tol=0.05, method="rejection",
                    corr=TRUE)
summary(model.sel)

## Goodness of fit the models
summary(
    gfit(target = target[c(1,3)],
                    sumstat = all.sims[sim.ids=="LSclump",c(1,3)],
                    nb.replicate = 200, tol = 0.05)
)

summary(
    gfit(target = target,
                    sumstat = all.sims[sim.ids=="LSclump",],
                    nb.replicate = 200, tol = 0.05)
    )


## Parameter estimation ##
## Cross validation
cv.nn <- cv4abc(param=sim.y[sim.ids=="LSclump"],
       sumstat = all.sims[sim.ids=="LSclump",],
       tols=c(0.05, 0.025, 0.01), nval=30, method="neuralnet")
cv.rej <- cv4abc(param=sim.y[sim.ids=="LSclump"],
       sumstat = all.sims[sim.ids=="LSclump",],
       tols=c(0.05, 0.025, 0.01), nval=30, method="rejection")
cv.ll <- cv4abc(param=data.frame(S=sim.y[sim.ids=="LSclump"]),
       sumstat = all.sims[sim.ids=="LSclump",],
       tols=c(0.05, 0.025, 0.01), nval=30, method="ridge")
par(mfrow=c(1,3))
plot(cv.nn, caption="Neural net")
plot(cv.rej, caption="Rejection")
plot(cv.ll, caption="Ridge regression")
par(mfrow=c(1,1))

## Posterior distribution of Species richness from the selected model
S.post1 <- abc(target = target, param=data.frame(S=sim.y[sim.ids=="LSclump"]),
              sumstat = all.sims[sim.ids=="LSclump",],
              tol=0.025, method="rejection")
summary(S.post1)
hist(S.post1)

## Posterior predictive check
Y <- mclapply(sample(S.post1$unadj.values,100, replace=TRUE),
                  f1, lower=1e-20, upper=1e20, mc.cores=4)
tmp1 <- Y[[1]][2,]
for(i in 2:length(Y))
    tmp1 <- rbind(tmp1, Y[[i]][2,])
par(mfrow=c(2,2))
for(i in 1:4){
    hist(tmp1[,i], xlim=range(c(tmp1[,i],target[i])), main=names(Y[[1]])[i])
    abline(v=target[i], col="red")
}
par(mfrow=c(1,1))
         
