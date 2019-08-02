source("functions.R")
load("lists_with_all_objects.RData")
library(dplyr)
library(tidyr)
library(VCA)

## Read table with all estimates, original and bias corrected
S.estimates.all <- read.csv("figs_and_tables/estimates_S_table.csv")

## Minimum and amximum values of estimates assuming clumped and random sampling in bias correction   
S.estimates.all %>%
    filter(bias.corrected==TRUE&type!="LSE TNB") %>%
    mutate(difIC = (IC.up-IC.low)/IC.low)%>%
    group_by(sampling, type) %>%
    summarise(mini=min(mean), maxi = max(mean), difp = (maxi-mini)/mini, dif.ICm=min(difIC))
## Same for each data set
S.estimates.all %>%
    filter(bias.corrected==TRUE&type!="LSE TNB") %>%
    mutate(difIC = (IC.up-IC.low)/IC.low)%>%
    group_by(sampling, dataset) %>%
    summarise(mini=min(mean, na.rm=TRUE), maxi = max(mean,na.rm=TRUE), difp = (maxi-mini)/mini, dif.ICm=min(difIC))
## Same, only for sampling assumed
S.estimates.all %>%
    filter(bias.corrected==TRUE&type!="LSE TNB") %>%
    mutate(difIC = (IC.up-IC.low)/IC.low)%>%
    group_by(sampling) %>%
    summarise(mini=min(mean, na.rm=TRUE), maxi = max(mean,na.rm=TRUE), difp = (maxi-mini)/mini, dif.ICm=min(difIC))

## Final figures : average of estimates corrected for bias assuming clumped sample ##
## CHAO excluded as it had a bias above 250% and an uncorrected value below the know number of species

## A weighted mean of estimates and upper and lower limits of CI's
S.estimates.all %>%
    filter(bias.corrected==TRUE&type!="LSE TNB"&sampling=="clump"&type!="TNB"&type!="CHAO") %>%
    mutate(se = (IC.up-IC.low)/4)%>%
    group_by(dataset) %>%
    summarise(w.mean = weighted.mean(mean, w=1/se), min=min(IC.low), max=max(IC.up))

S.estimates.all %>%
    filter(bias.corrected==TRUE&type!="LSE TNB"&sampling=="clump"&type!="CHAO") %>%
    mutate(se = (IC.up-IC.low)/4)%>%
    group_by(dataset) %>%
    summarise(w.mean = weighted.mean(mean, w=1/se), min=min(IC.low), max=max(IC.up))

## Max, min IC limits for each data set
S.estimates.all %>%
    filter(bias.corrected==TRUE&type!="LSE TNB"&sampling=="clump"&type!="CHAO") %>%
    group_by(dataset) %>%
    summarise(min.low=min(IC.low), max.up=max(IC.up))
## Excluding TNB
S.estimates.all %>%
    filter(bias.corrected==TRUE&type!="LSE TNB"&sampling=="clump"&type!="TNB"&type!="CHAO") %>%
    group_by(dataset) %>%
    summarise(min.low=min(IC.low), max.up=max(IC.up))

## Variance partition
## Including TNB
S.estimates.all %>%
    filter(bias.corrected==TRUE&type!="LSE TNB"&sampling=="clump"&type!="CHAO") %>%
    anovaVCA(mean~dataset + type, Data=.)
## Excluding TNB
S.estimates.all %>%
    filter(bias.corrected==TRUE&type!="LSE TNB"&sampling=="clump"&type!="TNB"&type!="CHAO") %>%
    anovaVCA(mean~dataset + type, Data=.)
## Including bias as a component (only for LS, LSE, TNB)
S.estimates.all %>%
    filter(type!="LSE TNB"&type!="CHAO"&type!="ABC"&(sampling=="clump"|is.na(sampling))) %>%
    select(dataset, type, mean, bias.corrected) %>%
    mutate(type = ifelse(type=="LSE LS", "LSE", as.character(type)),
           bias.corrected = ifelse(bias.corrected, "corrected", "uncorrected")) %>%
    anovaVCA(mean~dataset + type + bias.corrected, Data=.)

## Contrasting 2013 and 2019 data sets
S.estimates.all %>%
    filter(bias.corrected==TRUE&type!="LSE TNB"&sampling=="clump"&type!="CHAO") %>%
    select(dataset, type, mean) %>%
    spread(dataset, mean) %>%
    mutate(d13.13t = (.[,"2013 updated"] -.[,"2013"])/.[,"2013"],
           d13t.19 = (.[,"2019"] -.[,"2013 updated"])/.[,"2013 updated"])
## Uncorrected
S.estimates.all %>%
    filter(bias.corrected==FALSE&type!="LSE TNB"&type!="CHAO")%>%
    select(dataset, type, mean) %>%
    spread(dataset, mean) %>%
    mutate(d13.13t = (.[,"2013 updated"] -.[,"2013"])/.[,"2013"],
           d13t.19 = (.[,"2019"] -.[,"2013 updated"])/.[,"2013 updated"])

## Checking bias in each data set
S.estimates.all %>%
    filter(type!="LSE TNB"&type!="CHAO"&type!="ABC"&(sampling=="clump"|is.na(sampling))) %>%
    select(dataset, type, mean, bias.corrected) %>%
    mutate(type = ifelse(type=="LSE LS", "LSE", as.character(type)),
           bias.corrected = ifelse(bias.corrected, "corrected", "uncorrected")) %>%
    spread(bias.corrected, mean)%>%
    mutate(bias = round((corrected-uncorrected)/corrected*100)) %>%
    arrange(type, dataset)
