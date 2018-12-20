load("abc_4000_mais_1000_S_D_lmean_lsd.RData")
LS.sims2 <- LS.sims
NB.sims2 <- NB.sims
simulated.vals2 <- simulated.vals
LS.index2 <- sapply(LS.sims2, function(x) !any(is.na(x)))
j12 <- (1:length(LS.sims2))[LS.index2]
all.sims2 <- LS.sims2[[min(j12)]]
for(i in j12[-min(j12)])
    all.sims2 <- rbind(all.sims2, LS.sims2[[i]])
## Simulations from TNB rad
## Excluding NAs
NB.index2 <- sapply(NB.sims2, function(x) !any(is.na(x)))
j22 <- (1:length(NB.sims2))[NB.index2]
for(i in j22)
    all.sims2 <- rbind(all.sims2, NB.sims2[[i]])

## Labels for each simulation
sim.ids2 <- c(rep(c("LSrnd", "LSclump"), sum(LS.index2)),
             rep(c("NBrnd", "NBclump"), sum(NB.index2)))## Simulated values for each simulation
sim.y2 <- c( rep(simulated.vals[j12],each=2),
           rep(simulated.vals[j22],each=2))
all.sims <- rbind(all.sims,all.sims2)
sim.ids <- c(sim.ids, sim.ids2)
sim.y <- c(sim.y, sim.y2)
