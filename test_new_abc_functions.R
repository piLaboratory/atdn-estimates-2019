## Regression of log(sd) ~ log(pop sizes)
plot(log(pop.sd) ~ log(population), data = atdn.13$data)
lm.sd <- lm(log(pop.sd) ~ log(population), data = atdn.13$data)
abline(lm.sd, col="blue")
lm.sd.sigma <- summary(lm.sd)$sigma
lmean.sd <- predict(lm.sd, newdata = data.frame(population=atdn.13$reg.ls.rad))

## Test sample functions
teste.nb <- with(atdn.13,
              NB.samp( rad = reg.ls.rad, tot.area = Tot.A, n.plots = N.plots,
                      lmean.k = reg.ls.rad.lk, lsd.k = lm.k.sigma,
                      lmean.sd = lmean.sd, lsd.sd = lm.sd.sigma , nrep =1) )
teste.p <- with(atdn.13,
              Pois.samp( rad = reg.ls.rad, tot.area = Tot.A, n.plots = N.plots,
                        lmean.sd = lmean.sd, lsd.sd = lm.sd.sigma , nrep =1) )
plot(rad(atdn.13$data$population), col="grey")
lines(rad(teste.p[,1]))
lines(rad(teste.p[,2]), lty=2)
lines(rad(teste.nb[,1]), col="red")
lines(rad(teste.nb[,2]), col="red", lty=2)
plot(with.est.error ~ no.est.error, data=teste.p, log="xy")
abline(0,1)
plot(with.est.error ~ no.est.error, data=teste.nb, log="xy")
abline(0,1)

## Test abc.function

teste.abc.ls <- sim.abc(S = 15000, N = atdn.13$Tot.t, sad = "ls", tot.area = atdn.13$Tot.A, n.plots = atdn.13$N.plots,
                     lm.sd.fit = lm.sd, lmk.fit = atdn.13$lm.k, nb.fit = atdn.13$y.nb2, nrep = 3)
