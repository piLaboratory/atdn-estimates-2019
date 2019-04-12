source("functions.R")

## Generates a list with all objects needed for the calculations (and some garbage too)


## Original 2013 data set
atdn.13 <- atdn.estimates(path.to.data = "science_appendix2013.csv",
                          Tot.t = 3.9e11,
                          Tot.A = 6.29e8,
                          N.plots = 1170,
                          Samp.A = 1109.31)
## Fix NA in TNB IC, missing because upper lower limit of size CI is NA
## Susbtitutin by 2 x standard error of the coefficient
atdn.13$tovo.S$CIs[4,1] <-
    with(atdn.13,
         tovo(fit = y.nb2,
              cf = c( size = summary(y.nb2)@coef[1,1] - 2*summary(y.nb2)@coef[1,2], mu = tovo.S$CIs[2,1]),
              p = p1)
         )

## 2013 data set with reviewed taxonomy (as 2019)
atdn.13.tax <- atdn.estimates(path.to.data = "Populations_2013_tax_2019.csv",
                          Tot.t = 3.9e11,
                          Tot.A = 5.79e8,
                          N.plots = 1153,
                          Samp.A = 1097,
                          lm.sd.fit = atdn.13$lm.sd)
## 2019 Data set
atdn.19 <- atdn.estimates(path.to.data = "Populations_2019_V1.csv",
                          Tot.t = 3.06e11,
                          Tot.A = 5.79e8,
                          N.plots = 1946,
                          Samp.A = 2038,
                          lm.sd.fit = atdn.13$lm.sd)

save(atdn.13, atdn.19, atdn.13.tax, file = "lists_with_all_objects.RData")
