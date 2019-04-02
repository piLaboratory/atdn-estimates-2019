source("functions.R")
load("lists_with_all_objects.RData")
load("abcSummaries.RData")
library(ggplot2)

## Table with all estimates of species richness and CI's 95%
S.estimates <- expand.grid(
    dataset = c("2013", "2013t", "2019"),
    type = c("tnb", "ls", "LSE", "LSE bc", "ABC"),
    mean = NA,
    IC.low = NA,
    IC.up = NA
)
## Including the values
## Truncated negative binomial
S.estimates[S.estimates$type=="tnb"&S.estimates$dataset=="2013",3:5] <-
    with(atdn.13, c(tovo.S$S.est, tovo.S$CIs[4,2:1]))
S.estimates[S.estimates$type=="tnb"&S.estimates$dataset=="2013t",3:5] <-
    with(atdn.13.tax, c(tovo.S$S.est, tovo.S$CIs[4,2:1]))
S.estimates[S.estimates$type=="tnb"&S.estimates$dataset=="2019",3:5] <-
    with(atdn.19, c(tovo.S$S.est, tovo.S$CIs[4,2:1]))
## Log-series
S.estimates[S.estimates$type=="ls"&S.estimates$dataset=="2013",3:5] <-
    with(atdn.13, c(S.ls, S.ls.ci))
S.estimates[S.estimates$type=="ls"&S.estimates$dataset=="2013t",3:5] <-
    with(atdn.13.tax, c(S.ls, S.ls.ci))
S.estimates[S.estimates$type=="ls"&S.estimates$dataset=="2019",3:5] <-
    with(atdn.19, c(S.ls, S.ls.ci))
## Linear extension of RAD of estimated pop sizes (LSE)
S.estimates[S.estimates$type=="LSE"&S.estimates$dataset=="2013",3:5] <-
    with(atdn.13, S.r.ls.boot["LS rnd",])
S.estimates[S.estimates$type=="LSE"&S.estimates$dataset=="2013t",3:5] <-
    with(atdn.13.tax, S.r.ls.boot["LS rnd",])
S.estimates[S.estimates$type=="LSE"&S.estimates$dataset=="2019",3:5] <-
    with(atdn.19, S.r.ls.boot["LS rnd",])
## Linear extension of RAD of estimated pop sizes (LSE), bias corrected
S.estimates[S.estimates$type=="LSE bc"&S.estimates$dataset=="2013",3:5] <-
    with(atdn.13, predict(ulrich.bias.fit, data.frame(S.est = S.ulrich$S[1,1]), interval = "confidence"))
S.estimates[S.estimates$type=="LSE bc"&S.estimates$dataset=="2013t",3:5] <-
    with(atdn.13.tax, predict(ulrich.bias.fit, data.frame(S.est = S.ulrich$S[1,1]), interval = "confidence"))
S.estimates[S.estimates$type=="LSE bc"&S.estimates$dataset=="2019",3:5] <-
    with(atdn.19, predict(ulrich.bias.fit, data.frame(S.est = S.ulrich$S[1,1]), interval = "confidence"))
## ABC
S.estimates[S.estimates$type=="ABC"&S.estimates$dataset=="2013",3:5] <-
    summary(abc2013.summ$S.post1)[c(4,2,6),]
S.estimates[S.estimates$type=="ABC"&S.estimates$dataset=="2013t",3:5] <-
    summary(abc2013t.summ$S.post1)[c(4,2,6),]
S.estimates[S.estimates$type=="ABC"&S.estimates$dataset=="2019",3:5] <-
    summary(abc2019.summ$S.post1)[c(4,2,6),]
## ## Shen & He estimates
## S.estimates[S.estimates$type=="ShenHe"&S.estimates$dataset=="2013",3:5] <-
##     with(atdn.13, S.Shen.boot["LS rnd",-4])
## S.estimates[S.estimates$type=="ShenHe"&S.estimates$dataset=="2013t",3:5] <-
##     with(atdn.13.tax, S.Shen.boot["LS rnd",-4])
## S.estimates[S.estimates$type=="ShenHe"&S.estimates$dataset=="2019",3:5] <-
##     with(atdn.19, S.Shen.boot["LS rnd",-4])
## ## Hui estimates
## S.estimates[S.estimates$type=="Hui"&S.estimates$dataset=="2013",3:5] <-
##     with(atdn.13, S.orc.boot["LS rnd",-4])
## S.estimates[S.estimates$type=="Hui"&S.estimates$dataset=="2013t",3:5] <-
##     with(atdn.13.tax, S.orc.boot["LS rnd",-4])
## S.estimates[S.estimates$type=="Hui"&S.estimates$dataset=="2019",3:5] <-
##     with(atdn.19, S.orc.boot["LS rnd",-4])

    
p1 <- S.estimates %>%
    ggplot(aes(dataset, mean, group=type)) +
    geom_point(aes(colour=type),size=5) +
    geom_line(aes(colour=type)) +
    geom_linerange(aes(ymin=IC.low, ymax=IC.up, colour=type), size=3, alpha=0.25) +
    theme_bw()

p1 + facet_wrap(~type)

## Dotplot
with(S.estimates,
     dotchart(mean, labels = dataset, groups = type, color=1:3, pch=19,
              xlim=range(c(IC.low, IC.up),na.rm=TRUE))
     )
Ys <- 23:25 - rep(seq(0,22, by=5), each=3)
cores <- rep(1:3,5)
for(i in 1:nrow(S.estimates))
    segments(x0=S.estimates$IC.low[i], x1=S.estimates$IC.up[i], y0=Ys[i], y1=Ys[i], col=cores[i])

## Population rad expected by ABC
## 2013
S1 <- S.estimates[S.estimates$type=="ABC"&S.estimates$dataset=="2013",3:5]
abc13.sim.pop <- with(atdn.13,
                 lapply(S1, sim.abc, sad = "ls",
                        N = Tot.t,  tot.area = Tot.A, n.plots = N.plots,
                        lm.sd.fit = lm.sd, lmk.fit = lm.k, nb.fit =y.nb2,
                        summary = FALSE, upper=1e16) )
## Rad plots for clumped samples
plot(rad(atdn.13$data$population), xlim=c(1,6e3), col="grey")
lines(rad(abc13.sim.pop$mean$clump.samp[,2]))
lines(rad(abc13.sim.pop$IC.low$clump.samp[,2]), lty=2)
lines(rad(abc13.sim.pop$IC.up$clump.samp[,2]), lty=2)
## 2013 with updated taxonomy
abc13t.sim.pop <- with(atdn.13.tax,
                      lapply(S.estimates[S.estimates$type=="ABC"&S.estimates$dataset=="2013t",3:5],
                             sim.abc, sad = "ls",
                             N = Tot.t,  tot.area = Tot.A, n.plots = N.plots,
                             lm.sd.fit = lm.sd, lmk.fit = lm.k, nb.fit =y.nb2,
                             summary = FALSE, upper=1e16) )
## Rad plots for clumped samples
plot(rad(atdn.13.tax$data$population), xlim=c(1,6e3), col="grey")
lines(rad(abc13t.sim.pop$mean$clump.samp[,2]))
lines(rad(abc13t.sim.pop$IC.low$clump.samp[,2]), lty=2)
lines(rad(abc13t.sim.pop$IC.up$clump.samp[,2]), lty=2)
## 2019 
abc19.sim.pop <- with(atdn.19,
                      lapply(S.estimates[S.estimates$type=="ABC"&S.estimates$dataset=="2019",3:5],
                             sim.abc, sad = "ls",
                             N = Tot.t,  tot.area = Tot.A, n.plots = N.plots,
                             lm.sd.fit = lm.sd, lmk.fit = lm.k, nb.fit =y.nb2,
                             summary = FALSE, upper=1e16) )
## Rad plots for clumped samples
plot(rad(atdn.19$data$population), xlim=c(1,6e3), col="grey")
lines(rad(abc19.sim.pop$mean$clump.samp[,2]))
lines(rad(abc19.sim.pop$IC.low$clump.samp[,2]), lty=2)
lines(rad(abc19.sim.pop$IC.up$clump.samp[,2]), lty=2)

