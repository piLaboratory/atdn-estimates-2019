source("functions.R")
load("lists_with_all_objects.RData")
library(dplyr)
library(tidyr)

S.estimates <- read.csv("figs_and_tables/estimates_S_table.csv")

S.estimates %>%
    filter(bias.corrected==TRUE&type!="LSE TNB") %>%
    mutate(difIC = (IC.up-IC.low)/IC.low)%>%
    group_by(sampling, type) %>%
    summarise(mini=min(mean), maxi = max(mean), difp = (maxi-mini)/mini, dif.ICm=min(difIC))

S.estimates %>%
    filter(bias.corrected==TRUE&type!="LSE TNB") %>%
    mutate(difIC = (IC.up-IC.low)/IC.low)%>%
    group_by(sampling, dataset) %>%
    summarise(mini=min(mean, na.rm=TRUE), maxi = max(mean,na.rm=TRUE), difp = (maxi-mini)/mini, dif.ICm=min(difIC))

S.estimates %>%
    filter(bias.corrected==TRUE&type!="LSE TNB") %>%
    mutate(difIC = (IC.up-IC.low)/IC.low)%>%
    group_by(sampling) %>%
    summarise(mini=min(mean, na.rm=TRUE), maxi = max(mean,na.rm=TRUE), difp = (maxi-mini)/mini, dif.ICm=min(difIC))
## A weighted mean of estimates and upper and lower limits of CI's
## Excluding TNB
S.estimates %>%
    filter(bias.corrected==TRUE&type!="LSE TNB"&sampling=="clump"&type!="TNB") %>%
    mutate(se = (IC.up-IC.low)/4)%>%
    group_by(dataset) %>%
    summarise(w.mean = weighted.mean(mean, w=1/se), min=min(IC.low), max=max(IC.up))

## Max, min IC limits for each data set exluding LSE TNB)
S.estimates %>%
    filter(bias.corrected==TRUE&type!="LSE TNB"&sampling=="clump") %>%
    group_by(dataset) %>%
    summarise(min.low=min(IC.low), max.up=max(IC.up))
## Max, min IC limits for each data set exluding TNB)
S.estimates %>%
    filter(bias.corrected==TRUE&type!="LSE TNB"&sampling=="clump"&type!="TNB") %>%
    group_by(dataset) %>%
    summarise(min.low=min(IC.low), max.up=max(IC.up))

## Variance partition
S.estimates %>%
    filter(bias.corrected==TRUE&type!="LSE TNB"&sampling=="clump") %>%
    anovaVCA(mean~dataset + type, Data=.)

S.estimates %>%
    filter(bias.corrected==TRUE&type!="LSE TNB"&sampling=="clump"&type!="TNB") %>%
    anovaVCA(mean~dataset + type, Data=.)
