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


## Simulated values
nsims <- 4000
simulated.vals <- runif(nsims, 1e4, S.orc$S.est)
## Runs the simulations with mclapply
LS.sims <- mclapply(simulated.vals, f1, mc.cores=4)
save.image()
NB.sims <- mclapply(simulated.vals, f1, upper=1e20, LS = FALSE, mc.cores=4)
save.image()


## Assembles all simulation results in a matrix
## Simulations from LS rads
## check for failed simulations
LS.index <- sapply(LS.sims, function(x) !any(is.na(x)))
j1 <- (1:length(LS.sims))[LS.index]
all.sims <- LS.sims[[min(j1)]]
for(i in j1[-min(j1)])
    all.sims <- rbind(all.sims, LS.sims[[i]])
## Simulations from TNB rad
## Excluding NAs
NB.index <- sapply(NB.sims, function(x) !any(is.na(x)))
j2 <- (1:length(NB.sims))[NB.index]
for(i in j2)
    all.sims <- rbind(all.sims, NB.sims[[i]])
## Labels for each simulation
sim.ids <- c(rep(c("LSrnd", "LSclump"), sum(LS.index)),
             rep(c("NBrnd", "NBclump"), sum(NB.index)))
## Simulated values for each simulation
sim.y <- c( rep(simulated.vals[j1],each=2),
           rep(simulated.vals[j2],each=2))
           
             
             
## ABC ##
## Model selection
## Target: observed number of species, lmean, sdmean and zero of Mean_square with obs values             
target <- c(Sobs, D(Dec2018$population), mean(log(Dec2018$population)), sd(log(Dec2018$population)))

## Quick diagnostics plots
par(mfrow=c(2,2))
for(i in 1:ncol(all.sims)){
    boxplot(all.sims[,i]~sim.ids, main=colnames(all.sims)[i])
    abline(h=target[i], lty=2, col="blue")
    }
par(mfrow=c(1,1))
## S in the sample x S total for each model
plot(all.sims[sim.ids=="LSrnd","S"]~simulated.vals[j1],
     xlim=range(simulated.vals), ylim=range(all.sims[,"S"]), cex=0.5)
points(all.sims[sim.ids=="LSclump","S"]~simulated.vals[j1], col="blue", cex=0.5)
points(all.sims[sim.ids=="NBrnd","S"]~simulated.vals[j2], col="red", cex=0.5)
points(all.sims[sim.ids=="NBclump","S"]~simulated.vals[j2], col="green", cex=0.5)
abline(h=nrow(Dec2018), lty=2)
legend("topright", unique(sim.ids), pch=1,
       col=c("black", "blue", "red", "green"), bty="n")
## logmean x S
plot(all.sims[sim.ids=="LSrnd","lmean"]~simulated.vals[j1],
     xlim=range(simulated.vals), ylim=range(all.sims[,"lmean"]), cex=0.5)
points(all.sims[sim.ids=="LSclump","lmean"]~simulated.vals[j1], col="blue", cex=0.5)
points(all.sims[sim.ids=="NBrnd","lmean"]~simulated.vals[j2], col="red", cex=0.5)
points(all.sims[sim.ids=="NBclump","lmean"]~simulated.vals[j2], col="green", cex=0.5)
abline(h=mean(log(Dec2018$population)), lty=2)
legend("topright", unique(sim.ids), pch=1,
       col=c("black", "blue", "red", "green"), bty="n")
## logsd x S
plot(all.sims[sim.ids=="LSrnd","lsd"]~simulated.vals[j1],
     xlim=range(simulated.vals), ylim=range(all.sims[,"lsd"]), cex=0.5)
points(all.sims[sim.ids=="LSclump","lsd"]~simulated.vals[j1], col="blue", cex=0.5)
points(all.sims[sim.ids=="NBrnd","lsd"]~simulated.vals[j2], col="red", cex=0.5)
points(all.sims[sim.ids=="NBclump","lsd"]~simulated.vals[j2], col="green", cex=0.5)
abline(h=sd(log(Dec2018$population)), lty=2)
legend("topright", unique(sim.ids), pch=1,
       col=c("black", "blue", "red", "green"), bty="n")
## Simulated S x logmean
plot(lmean~S, data=all.sims, type="n")
points(lmean~S, data=all.sims[sim.ids=="LSrnd",], cex=0.5)
points(lmean~S, data=all.sims[sim.ids=="LSclump",], cex=0.5, col="blue")
points(lmean~S, data=all.sims[sim.ids=="NBrnd",], cex=0.5, col="red")
points(lmean~S, data=all.sims[sim.ids=="NBclump",], cex=0.5, col="green")
points(nrow(Dec2018), mean(log(Dec2018$population)), pch=3, cex=2)
legend("topright", unique(sim.ids), pch=1,
       col=c("black", "blue", "red", "green"), bty="n")
## Simulated S x logsd
plot(lsd~S, data=all.sims, type="n")
points(lsd~S, data=all.sims[sim.ids=="LSrnd",], cex=0.5)
points(lsd~S, data=all.sims[sim.ids=="LSclump",], cex=0.5, col="blue")
points(lsd~S, data=all.sims[sim.ids=="NBrnd",], cex=0.5, col="red")
points(lsd~S, data=all.sims[sim.ids=="NBclump",], cex=0.5, col="green")
points(nrow(Dec2018), sd(log(Dec2018$population)), pch=19, cex=2, col="orange")
legend("topright", unique(sim.ids), pch=1,
       col=c("black", "blue", "red", "green"), bty="n")
## lmean x lsd
plot(lsd~lmean, data=all.sims, type="n")
points(lsd~lmean, data=all.sims[sim.ids=="LSrnd",], cex=0.5)
points(lsd~lmean, data=all.sims[sim.ids=="LSclump",], cex=0.5, col="blue")
points(lsd~lmean, data=all.sims[sim.ids=="NBrnd",], cex=0.5, col="red")
points(lsd~lmean, data=all.sims[sim.ids=="NBclump",], cex=0.5, col="green")
points(mean(log(Dec2018$population)), sd(log(Dec2018$population)), pch=19, cex=2, col="orange")
legend("topleft", unique(sim.ids), pch=1,
       col=c("black", "blue", "red", "green"), bty="n")
## S x D
plot(D~S, data=all.sims, type="n", ylim=range(c(all.sims$D, D(Dec2018$population))))
points(D~S, data=all.sims[sim.ids=="LSrnd",], cex=0.5)
points(D~S, data=all.sims[sim.ids=="LSclump",], cex=0.5, col="blue")
points(D~S, data=all.sims[sim.ids=="NBrnd",], cex=0.5, col="red")
points(D~S, data=all.sims[sim.ids=="NBclump",], cex=0.5, col="green")
points(Sobs, D(Dec2018$population), pch=19, cex=2, col="orange")
legend("topleft", unique(sim.ids), pch=1,
       col=c("black", "blue", "red", "green"), bty="n")

## Model selection
## Cross-validation
cv.modsel <- cv4postpr(sim.ids, all.sims, nval=100, tol=0.01, method="rejection")
summary(cv.modsel)
## Model selection
model.sel <- postpr(target = target[c(1,3,4)],
                    index=sim.ids,
                    sumstat = all.sims[,c(1,3,4)],
                    tol=0.02, method="rejection",
                    corr=TRUE)
summary(model.sel)

## Goodness of fit the models
NBclump.gof <- gfit(target = target[c(1,3,4)],
                    sumstat = all.sims[sim.ids=="NBclump",c(1,3,4)],
                    nb.replicate = 200,
                    tol = 0.02)
summary(NBclump.gof)

## Parameter estimation ##
## Cross validation
cv.nn <- cv4abc(param=simulated.vals[j1],
       sumstat = all.sims[sim.ids=="LSclump",],
       tols=c(0.05, 0.025, 0.01), nval=30, method="neuralnet")
cv.rej <- cv4abc(param=simulated.vals[j1],
       sumstat = all.sims[sim.ids=="LSclump",],
       tols=c(0.05, 0.025, 0.01), nval=30, method="rejection")
cv.ll <- cv4abc(param=data.frame(S=simulated.vals[j1]),
       sumstat = all.sims[sim.ids=="LSclump",],
       tols=c(0.05, 0.025, 0.01), nval=30, method="ridge")
par(mfrow=c(1,3))
plot(cv.nn, caption="Neural net")
plot(cv.rej, caption="Rejection")
plot(cv.ll, caption="Ridge regression")
par(mfrow=c(1,1))

## Posterior distribution of Species richness from the selected model
S.post1 <- abc(target = target[c(1,3,4)], param=data.frame(S=simulated.vals[j2]),
              sumstat = all.sims[sim.ids=="NBclump",c(1,3,4)],
              tol=0.1, method="rejection")
summary(S.post1)
hist(S.post1)

## Posterior predictive check
LS.clump.pc <- mclapply(sample(S.post1$unadj.values,100, replace=TRUE), f1, mc.cores=4)
tmp1 <- LS.clump.pc[[1]]
for(i in 2:length(LS.clump.pc))
    tmp1 <- rbind(tmp1, LS.clump.pc[[i]])
par(mfrow=c(2,2))
for(i in 1:ncol(tmp1)){
    hist(tmp1[,i], xlim=range(c(tmp1[,i],target[i])))
    abline(v=target[i], col="red")
}
         
