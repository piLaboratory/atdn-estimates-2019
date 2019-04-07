source("functions.R")
load("lists_with_all_objects.RData")
load("abcSummaries.RData")

library(ggplot2)

################################################################################
## Estimation bias of ls and tnb
################################################################################

nf <- layout(matrix(c(1:3,7:9,4:6,10:12),2,6,byrow = TRUE),
             rep(c(0.25,0.25,3),4), rep(3,12))
## 2013
par(mar = c(5, 0, 4, 0))
with(subset(bias13$ls$estimates, S.est.clump>atdn.13$S.ls.ci[1] & S.est.clump<atdn.13$S.ls.ci[2]),
     boxplot(S, axes=FALSE, ylim = atdn.13$S.ls.ci*c(.9,1.1), border="red"))
with(subset(bias13$ls$estimates, S.est.rnd>atdn.13$S.ls.ci[1] & S.est.rnd<atdn.13$S.ls.ci[2]),
     boxplot(S, axes=FALSE, ylim = atdn.13$S.ls.ci*c(.9,1.1), border="blue"))
par(mar = c(5, 4, 4, 5))
plot(S ~ S.est.rnd, data = bias13$ls$estimates, ylim = atdn.13$S.ls.ci*c(.9,1.1),
     xlim = atdn.13$S.ls.ci*c(.9,1.1), col="blue", ylab = "")
box()
axis(1)
points(S ~ S.est.clump, data = bias13$ls$estimates, col="red")
abline(0,1)
segments(x0 = c(atdn.13$S.ls.ci,0,0),
         y0 = c(0,0,atdn.13$S.ls.ci),
         x1 = c(atdn.13$S.ls.ci,atdn.13$S.ls.ci),
         y1 = c(atdn.13$S.ls.ci,atdn.13$S.ls.ci),
         lty=2)
##
par(mar = c(5, 0, 4, 0))
with(subset(bias13$tnb$estimates, S.est.clump>atdn.13$tovo$CIs[4,2] & S.est.clump<atdn.13$tovo$CIs[4,1]),
     boxplot(S, axes=FALSE, ylim = range(bias13$tnb$estimates$S), border="red"))
with(subset(bias13$tnb$estimates, S.est.rnd>atdn.13$tovo$CIs[4,2] & S.est.clump<atdn.13$tovo$CIs[4,1]),
     boxplot(S, axes=FALSE, ylim = range(bias13$tnb$estimates$S), border="blue"))
par(mar = c(5, 4, 4, 5))
plot(S ~ S.est.rnd, data = bias13$tnb$estimates,
     xlim = range(bias13$tnb$estimates$S)*c(.9,1.1),
     ylim = range(bias13$tnb$estimates$S),
     col="blue", ylab = "")
box()
axis(1)
points(S ~ S.est.clump, data = bias13$tnb$estimates, col="red")
abline(0,1)
segments(x0 = c(atdn.13$tovo$CIs[4,],0,0),
         y0 = c(0,0,atdn.13$tovo$CIs[4,]),
         x1 = c(atdn.13$tovo$CIs[4,],atdn.13$tovo$CIs[4,]),
         y1 = c(atdn.13$tovo$CIs[4,],atdn.13$tovo$CIs[4,]),
         lty=2)
## 2019
par(mar = c(5, 0, 4, 0))
with(subset(bias19$ls$estimates, S.est.clump>atdn.19$S.ls.ci[1] & S.est.clump<atdn.19$S.ls.ci[2]),
     boxplot(S, axes=FALSE, ylim = atdn.19$S.ls.ci*c(.9,1.1), border="red"))
with(subset(bias19$ls$estimates, S.est.rnd>atdn.19$S.ls.ci[1] & S.est.rnd<atdn.19$S.ls.ci[2]),
     boxplot(S, axes=FALSE, ylim = atdn.19$S.ls.ci*c(.9,1.1), border="blue"))
par(mar = c(5, 4, 4, 5))
plot(S ~ S.est.rnd, data = bias19$ls$estimates, ylim = atdn.19$S.ls.ci*c(.9,1.1),
     xlim = atdn.19$S.ls.ci*c(.9,1.1), col="blue", ylab = "")
box()
axis(1)
points(S ~ S.est.clump, data = bias19$ls$estimates, col="red")
abline(0,1)
segments(x0 = c(atdn.19$S.ls.ci,0,0),
         y0 = c(0,0,atdn.19$S.ls.ci),
         x1 = c(atdn.19$S.ls.ci,atdn.19$S.ls.ci),
         y1 = c(atdn.19$S.ls.ci,atdn.19$S.ls.ci),
         lty=2)
par(mar = c(5, 0, 4, 0))
with(subset(bias19$tnb$estimates, S.est.clump>atdn.19$tovo$CIs[4,2] & S.est.clump<atdn.19$tovo$CIs[4,1]),
     boxplot(S, axes=FALSE, ylim = range(bias19$tnb$estimates$S), border="red"))
with(subset(bias19$tnb$estimates, S.est.rnd>atdn.19$tovo$CIs[4,2] & S.est.clump<atdn.19$tovo$CIs[4,1]),
     boxplot(S, axes=FALSE, ylim = range(bias19$tnb$estimates$S), border="blue"))
par(mar = c(5, 4, 4, 5))
plot(S ~ S.est.rnd, data = bias19$tnb$estimates,
     xlim = range(bias19$tnb$estimates$S)*c(.9,1.1),
     ylim = range(bias19$tnb$estimates$S),
     col="blue", ylab = "")
box()
axis(1)
points(S ~ S.est.clump, data = bias19$tnb$estimates, col="red")
abline(0,1)
segments(x0 = c(atdn.19$tovo$CIs[4,],0,0),
         y0 = c(0,0,atdn.19$tovo$CIs[4,]),
         x1 = c(atdn.19$tovo$CIs[4,],atdn.19$tovo$CIs[4,]),
         y1 = c(atdn.19$tovo$CIs[4,],atdn.19$tovo$CIs[4,]),
         lty=2)

################################################################################
## Table with all estimates of species richness and CI's 95%
################################################################################
S.estimates <- expand.grid(
    dataset = c("2013", "2013 rev.", "2019"),
    type = c("TNB", "LS", "LSE bc", "ABC"),
    mean = NA,
    IC.low = NA,
    IC.up = NA
)
## Including the values
## Truncated negative binomial
S.estimates[S.estimates$type=="TNB"&S.estimates$dataset=="2013",3:5] <-
    with(atdn.13, c(tovo.S$S.est, tovo.S$CIs[4,2:1]))
S.estimates[S.estimates$type=="TNB"&S.estimates$dataset=="2013 rev.",3:5] <-
    with(atdn.13.tax, c(tovo.S$S.est, tovo.S$CIs[4,2:1]))
S.estimates[S.estimates$type=="TNB"&S.estimates$dataset=="2019",3:5] <-
    with(atdn.19, c(tovo.S$S.est, tovo.S$CIs[4,2:1]))
## Log-series
S.estimates[S.estimates$type=="LS"&S.estimates$dataset=="2013",3:5] <-
    with(atdn.13, c(S.ls, S.ls.ci))
S.estimates[S.estimates$type=="LS"&S.estimates$dataset=="2013 rev.",3:5] <-
    with(atdn.13.tax, c(S.ls, S.ls.ci))
S.estimates[S.estimates$type=="LS"&S.estimates$dataset=="2019",3:5] <-
    with(atdn.19, c(S.ls, S.ls.ci))
## Linear extension of RAD of estimated pop sizes (LSE)
## S.estimates[S.estimates$type=="LSE"&S.estimates$dataset=="2013",3:5] <-
##     with(atdn.13, S.r.ls.boot)
## S.estimates[S.estimates$type=="LSE"&S.estimates$dataset=="2013 rev.",3:5] <-
##     with(atdn.13.tax, S.r.ls.boot)
## S.estimates[S.estimates$type=="LSE"&S.estimates$dataset=="2019",3:5] <-
##     with(atdn.19, S.r.ls.boot)
## Linear extension of RAD of estimated pop sizes (LSE), bias corrected
S.estimates[S.estimates$type=="LSE bc"&S.estimates$dataset=="2013",3:5] <-
    with(atdn.13, predict(ulrich.bias.fit, data.frame(S.est = S.r.ls.boot)))
S.estimates[S.estimates$type=="LSE bc"&S.estimates$dataset=="2013 rev.",3:5] <-
    with(atdn.13.tax, predict(ulrich.bias.fit, data.frame(S.est = S.r.ls.boot)))
S.estimates[S.estimates$type=="LSE bc"&S.estimates$dataset=="2019",3:5] <-
    with(atdn.19, predict(ulrich.bias.fit, data.frame(S.est = S.r.ls.boot)) )
## ABC
S.estimates[S.estimates$type=="ABC"&S.estimates$dataset=="2013",3:5] <-
    summary(abc2013.summ$S.post1)[c(4,2,6),]
S.estimates[S.estimates$type=="ABC"&S.estimates$dataset=="2013 rev.",3:5] <-
    summary(abc2013t.summ$S.post1)[c(4,2,6),]
S.estimates[S.estimates$type=="ABC"&S.estimates$dataset=="2019",3:5] <-
    summary(abc2019.summ$S.post1)[c(4,2,6),]
## ## Shen & He estimates
## S.estimates[S.estimates$type=="ShenHe"&S.estimates$dataset=="2013",3:5] <-
##     with(atdn.13, S.Shen.boot["LS rnd",-4])
## S.estimates[S.estimates$type=="ShenHe"&S.estimates$dataset=="2013 rev.",3:5] <-
##     with(atdn.13.tax, S.Shen.boot["LS rnd",-4])
## S.estimates[S.estimates$type=="ShenHe"&S.estimates$dataset=="2019",3:5] <-
##     with(atdn.19, S.Shen.boot["LS rnd",-4])
## ## Hui estimates
## S.estimates[S.estimates$type=="Hui"&S.estimates$dataset=="2013",3:5] <-
##     with(atdn.13, S.orc.boot["LS rnd",-4])
## S.estimates[S.estimates$type=="Hui"&S.estimates$dataset=="2013 rev.",3:5] <-
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
     dotchart(mean, labels = dataset, groups = type, color=2:4, pch=1,
              pt.cex = 1.5,
              xlim=range(c(IC.low, IC.up),na.rm=TRUE))
     )
Ys <- 16:18 - rep(seq(0,15, by=5), each=3)
cores <- rep(2:4,5)
for(i in 1:nrow(S.estimates))
    segments(x0=S.estimates$IC.low[i], x1=S.estimates$IC.up[i], y0=Ys[i], y1=Ys[i],
             col=cores[i], lwd =1.5)


################################################################################
## Population rad expected by ABC ##
################################################################################
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

