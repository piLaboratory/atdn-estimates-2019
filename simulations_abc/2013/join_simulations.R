## Because computational limitations simulations for ABC of 2013 data were done in 2 pieces of 5000 replicates each.
## Each simulation was ran by a separate script (abc_run2013a.R, abc_run22013b.R) and has a separate output
## (ABC2013a.RData, ABC2013b.RData).
## This code join the resulting objects
## Load results
load("ABC2013b.RData")
all.sims2 <- all.sims
sim.ids2 <- sim.ids
sim.y2 <- sim.y
load("ABC2013a.RData")

## Joining simulation inputs and resultig objects in a list
abc2013 <- list( sims = rbind(all.sims, all.sims2),
                labels = c(sim.ids, sim.ids2),
                params = c(sim.y, sim.y2))
save(abc2013, file="abcFinal2013.RData")

