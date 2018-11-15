library(sads)
y <- read.table("abundVetcor.txt")[,1]
## Fit to Poisson-lognormal
y.pln <- fitpoilog(y)

## Fit to a ero-truncated negative binomial
y.nb <- fitnbi
