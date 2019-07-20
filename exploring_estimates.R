source("functions.R")
load("lists_with_all_objects.RData")
library(dplyr)
library(tidyr)
library(VCA)

## Read table with all estimates, original and bias corrected
S.estimates <- read.csv("figs_and_tables/estimates_S_table.csv")

## Minimum and amximum values of estimates assuming clumped and random sampling in bias correction   
S.estimates %>%
    filter(bias.corrected==TRUE&type!="LSE TNB") %>%
    mutate(difIC = (IC.up-IC.low)/IC.low)%>%
    group_by(sampling, type) %>%
    summarise(mini=min(mean), maxi = max(mean), difp = (maxi-mini)/mini, dif.ICm=min(difIC))
## Same for each data set
S.estimates %>%
    filter(bias.corrected==TRUE&type!="LSE TNB") %>%
    mutate(difIC = (IC.up-IC.low)/IC.low)%>%
    group_by(sampling, dataset) %>%
    summarise(mini=min(mean, na.rm=TRUE), maxi = max(mean,na.rm=TRUE), difp = (maxi-mini)/mini, dif.ICm=min(difIC))
## Same, only for sampling assumed
S.estimates %>%
    filter(bias.corrected==TRUE&type!="LSE TNB") %>%
    mutate(difIC = (IC.up-IC.low)/IC.low)%>%
    group_by(sampling) %>%
    summarise(mini=min(mean, na.rm=TRUE), maxi = max(mean,na.rm=TRUE), difp = (maxi-mini)/mini, dif.ICm=min(difIC))

## Final figures : average of estimates corrected for bias assuming clumped sample ##
## CHAO excluded as it had a bias above 250% and an uncorrected value below the know number of species

## A weighted mean of estimates and upper and lower limits of CI's
S.estimates %>%
    filter(bias.corrected==TRUE&type!="LSE TNB"&sampling=="clump"&type!="TNB"&type!="CHAO") %>%
    mutate(se = (IC.up-IC.low)/4)%>%
    group_by(dataset) %>%
    summarise(w.mean = weighted.mean(mean, w=1/se), min=min(IC.low), max=max(IC.up))

S.estimates %>%
    filter(bias.corrected==TRUE&type!="LSE TNB"&sampling=="clump"&type!="CHAO") %>%
    mutate(se = (IC.up-IC.low)/4)%>%
    group_by(dataset) %>%
    summarise(w.mean = weighted.mean(mean, w=1/se), min=min(IC.low), max=max(IC.up))


## Max, min IC limits for each data set
S.estimates %>%
    filter(bias.corrected==TRUE&type!="LSE TNB"&sampling=="clump"&type!="CHAO") %>%
    group_by(dataset) %>%
    summarise(min.low=min(IC.low), max.up=max(IC.up))
## Excluding TNB
S.estimates %>%
    filter(bias.corrected==TRUE&type!="LSE TNB"&sampling=="clump"&type!="TNB"&type!="CHAO") %>%
    group_by(dataset) %>%
    summarise(min.low=min(IC.low), max.up=max(IC.up))

## Variance partition
## Including TNB
S.estimates %>%
    filter(bias.corrected==TRUE&type!="LSE TNB"&sampling=="clump"&type!="CHAO") %>%
    anovaVCA(mean~dataset + type, Data=.)
## Excluding TNB
S.estimates %>%
    filter(bias.corrected==TRUE&type!="LSE TNB"&sampling=="clump"&type!="TNB"&type!="CHAO") %>%
    anovaVCA(mean~dataset + type, Data=.)
