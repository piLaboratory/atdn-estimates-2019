load("../../lists_with_all_objects.R")
needed.objs <- atdn.13.tax[c("pln.cf", "Tot.t", "N.plots", "Tot.A", "y.nb2", "lm.k", "lm.sd", "data")]
save(needed.objs, file = "needed_objs2013_tax.RData")
rm(list=ls())
## pln.cf <- atdn.13$pln.cf
## Tot.t <- atdn.13$Tot.t
## N.plots <- atdn.13$N.plots
## Tot.A <- atdn.13$Tot.A
## y.nb2 <- atdn.13$y.nb2
## lm.k <- atdn.13$lm.k
## lm.sd <- atdn.13$lm.sd
## data <- atdn.13$data

## save(pln.cf,
##      Tot.t,
##      N.plots,
##      Tot.A,
##      y.nb2,
##      lm.k,
##      lm.sd,
##      data , file="needed_objs2013.RData")

