## Because computational limitations simulations for ABC of 2019 data were done in 2 pieces of 5000 replicates each.
## Each simulation was ran by a separate script (abc_run2019a.R, abc_run22019b.R) and has a separate output
## (ABC2019a.RData, ABC2019b.RData).
## This code join the resulting objects
## Load results
load("ABC2019b.RData")
all.sims2 <- all.sims
sim.ids2 <- sim.ids
sim.y2 <- sim.y
## needed.objs2 <- needed.objs
load("ABC2019a.RData")

## Joining simulation inputs and resultig objects in a list
abc2019 <- list( sims = rbind(all.sims, all.sims2),
                labels = c(sim.ids, sim.ids2),
                params = c(sim.y, sim.y2))

save(abc2019, file="abcFinal2019.RData")

