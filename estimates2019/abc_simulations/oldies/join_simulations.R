## Because computational limitations simulations for ABC of 2019 data were done in 3 pieces (7500, 1500 and 1500 replicates)
## Each simulation was ran by a separate script (abc_run.R, abc_run2.R, abc_run3.R, abc_run4.R) and has a separate output
## (testeABC.RData, testeABC2.RData, testeABC3.RData, testeABC4.RData).
## This code join the resulting objects
## Load results
load("testeABC.RData")
load("testeABC2.RData")
load("testeABC3.RData")
load("testeABC4.RData")

## Joining simulation inputs and resultig objects in a list
abc2019 <- list( sims = rbind(all.sims, all.sims2, all.sims3, all.sims4),
                labels = c(sim.ids, sim.ids2, sim.ids3, sim.ids4),
                params = c(sim.y, sim.y2, sim.y3, sim.y4),
                input.objs = c(input.objs, input.objs2, input.objs3, input.objs4))
save(abc2019, file="abc2019.RData")

