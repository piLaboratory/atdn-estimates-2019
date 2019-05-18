source("functions.R")

## reorder columns, keep only need variables, rename variables
format.data(input.file ="AppendixS1-Species_data_2013.csv",
            output.file="2013_original.csv", as.is=TRUE)
format.data(input.file ="2013_taxonomy_2019.csv",
            output.file="2013_tax_2019.csv", as.is=TRUE)
format.data(input.file ="2019.csv",
            output.file="2019_May.csv", as.is=TRUE)

## Generates a list with all objects needed for the calculations (and some garbage too)


## Original 2013 data set
atdn.13 <- atdn.estimates(path.to.data = "2013_original.csv",
                          Tot.t = 3.9e11,
                          Tot.A = 1232100*567,
                          N.plots = 1170,
                          Samp.A = 1080)
## Fix NA in TNB IC, missing because upper lower limit of size CI is NA
## Susbtitutin by 2 x standard error of the coefficient
atdn.13$tovo.S$CIs[4,1] <-
    with(atdn.13,
         tovo(fit = y.nb2,
              cf = c( size = summary(y.nb2)@coef[1,1] - 2*summary(y.nb2)@coef[1,2], mu = tovo.S$CIs[2,1]),
              p = p1)
         )

## 2013 data set with reviewed taxonomy (as 2019)
atdn.13.tax <- atdn.estimates(path.to.data = "2013_tax_2019.csv",
                          Tot.t = 3.06e11,
                          Tot.A = 5.79e8,
                          N.plots = 1162,
                          Samp.A = 1080)
## 2019 Data set
atdn.19 <- atdn.estimates(path.to.data = "2019_May.csv",
                          Tot.t = 3.06e11,
                          Tot.A = 5.79e8,
                          N.plots = 1946,
                          Samp.A = 2042)

save(atdn.13, atdn.19, atdn.13.tax, file = "lists_with_all_objects.RData")
