##
source("functions.R")
atdn2013 <- read.csv2("science_appendix2013.csv", as.is=TRUE)
## 100 RADs from bootstrap samples of the population sizes, from a Gaussian
plot(rad(atdn2013$population), col="grey")
for(i in 1:100){
    y <- rnorm(length(atdn2013$population),
               mean = atdn2013$population,
               sd = atdn2013$pop.sd)
    lines(rad(y), lwd=0.5)
}

## Function "ulrich" calculates linear expansion of logseries (LSE)
## and also the lower limit based on the lognormal (LNE)
## The current version can return mean and CI's of the boostrap above
with(atdn2013,
     ulrich(x = population, x.mean = population, x.sd = pop.sd,
            boot=TRUE, n.boot = 200)$S
     )

## As expected from the plot above, the mean of LSE estimates from bootstrap sample is
## lower than the value we got from the extimated population sizes.
## Not sure yet, but I guess that this is because species switch ranks in bootstrap samples
