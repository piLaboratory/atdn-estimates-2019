## Because computational limitations simulations for ABC of 2013 with new taxonomy
## data were done in 2 pieces of 5000 replicates each.
## Each simulation was ran by a separate script (abc_run2013ta.R, abc_run22013tb.R) and has a separate output
## (ABC2013ta.RData, ABC2013tb.RData).
## This code join the resulting objects
## Load results
load("ABC2013tb.RData")
all.sims2 <- all.sims
sim.ids2 <- sim.ids
sim.y2 <- sim.y
load("ABC2013ta.RData")

## Joining simulation inputs and resultig objects in a list
abc2013t <- list( sims = rbind(all.sims, all.sims2),
                labels = c(sim.ids, sim.ids2),
                params = c(sim.y, sim.y2))
save(abc2013t, file="abcFinal2013tax2019.RData")

