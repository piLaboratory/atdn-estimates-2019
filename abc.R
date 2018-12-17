source("functions.R")
library(abc)
library(parallel)

## Simulations
## Function to parallelize simulations
f1 <- function(x, ...){
    y <- try(sim.abc(S = x, N = Tot.t, n.plots = N.plots, tot.area= Tot.A,
                     nb.fit=y.nb2, lmk = lm.k, obs.values=Dec2018$population,
                     nrep = 2, ...))
    if(class(y)=="try-error")
        return(matrix(NA, nrow=2, ncol=3))
    else
        return(y)
}
## Simulated values
nsims <- 3000
simulated.vals <- runif(nsims, 1e4, S.orc$S.est)
## Runs the simulations
LS.sims <- mclapply(simulated.vals, f1, mc.cores=3)
save.image()
NB.sims <- mclapply(simulated.vals, f1, upper=1e20, LS = FALSE, mc.cores=3)
save.image()

## Assembles all simulation results in a matrix
## Simulations from LS rads
all.sims <- LS.sims[[1]]
for(i in 2:length(LS.sims))
    all.sims <- rbind(all.sims, LS.sims[[i]])
## Simulations from TNB rad
## Excluding NAs
oksim <- sapply(NB.sims, function(x) !any(is.na(x)))
for(i in (1:length(NB.sims))[oksim])
    all.sims <- rbind(all.sims, NB.sims[[i]])
## Labels for each simulation
sim.ids <- c(rep(c("LSrnd", "LSclump"), length(LS.sims)),
             rep(c("NBrnd", "NBclump"), sum(oksim)))
             
## ABC ##
## Model selection
## Target: observed number of species, lmean, sdmean and zero of Mean_square with obs values             
target <- c(nrow(Dec2018), mean(log(Dec2018$population)), sd(log(Dec2018$population)), 0)

## Quick diagnostics plots
par(mfrow=c(2,2))
for(i in 1:ncol(all.sims)){
    boxplot(all.sims[,i]~sim.ids, main=colnames(all.sims)[i])
    abline(h=target[i], lty=2, col="blue")
    }
par(mfrow=c(1,1))
## S in the sample x S total for each model
plot(all.sims[sim.ids=="LSrnd","S"]~simulated.vals,
     xlim=range(simulated.vals), ylim=range(all.sims[,"S"]), cex=0.5)
points(all.sims[sim.ids=="LSclump","S"]~simulated.vals, col="blue", cex=0.5)
points(all.sims[sim.ids=="NBrnd","S"]~simulated.vals[oksim], col="red", cex=0.5)
points(all.sims[sim.ids=="NBclump","S"]~simulated.vals[oksim], col="green", cex=0.5)
abline(h=nrow(Dec2018), lty=2)
legend("topright", unique(sim.ids), pch=1,
       col=c("black", "blue", "red", "green"), bty="n")
## logmean x S
plot(all.sims[sim.ids=="LSrnd","lmean"]~simulated.vals,
     xlim=range(simulated.vals), ylim=range(all.sims[,"lmean"]), cex=0.5)
points(all.sims[sim.ids=="LSclump","lmean"]~simulated.vals, col="blue", cex=0.5)
points(all.sims[sim.ids=="NBrnd","lmean"]~simulated.vals[oksim], col="red", cex=0.5)
points(all.sims[sim.ids=="NBclump","lmean"]~simulated.vals[oksim], col="green", cex=0.5)
abline(h=mean(log(Dec2018$population)), lty=2)
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

## Model selection
## Cross-validation
cv.modsel <- cv4postpr(sim.ids, all.sims, nval=100, tol=0.2, method="rejection")
summary(cv.modsel)
## Model selection
model.sel <- postpr(target = target,
                    index=sim.ids,
                    sumstat = all.sims,
                    tol=0.2, method="rejection")
summary(model.sel)

## Posterior distribution of Species richness from the selected model
S.post <- abc(target = target, param=simulated.vals,
              sumstat = all.sims[sim.ids=="LSclump",],
              trasnf= c("log", "log", "none", "none"),
              tol=0.25, method="neuralnet")


tmp <- all.sims[sim.ids=="LSclump",]
tmp2 <- sweep(tmp, 2, apply(tmp,2,mean), "/")
summary(abc(target = target, param=simulated.vals,
              sumstat = tmp2,
              tol=0.1, method="neuralnet"))
