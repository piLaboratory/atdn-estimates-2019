## Part of old code of function "atdn.estimates' to build rads from LS using
## posterior estimates of the LS clump model. Onlye valid if this is the model selected by ABC.
## TBI: generelaize this code for any model (check if it is necessary).

## ----clumped LS RAD and CI-----------------------------------
        ## Predicted log(k) values for LS rad
        reg.ls.rad.lk <- predict(lm.k, 
                                 newdata=data.frame(dens.ha=reg.ls.rad/Tot.A))

        ## LS RAD with the lower CI
        abc.ls.c.rad.l <- rad.ls(S = S.post1.s[2,],
                                 N = Tot.t, 
                                 alpha = fishers.alpha(Tot.t, S.post1.s[2,]))$y
        ## LS RAD with the upper CI
        abc.ls.c.rad.u <- rad.ls(S = S.post1.s[6,],
                                 N = Tot.t, 
                                 alpha = fishers.alpha(Tot.t, S.post1.s[6,]))$y
        ## Simulated clumped samples
        ## from LS with species richness estimated from linear extrapolation
        ls.clump <- NB.samp(rad = reg.ls.rad, tot.area = Tot.A, 
                            n.plots = N.plots,
                            lmean.sd = reg.ls.rad.lsd,
                            lsd.sd = lm.sd.sigma,
                            lmean.k = reg.ls.rad.lk, 
                            lsd.k = lm.k.sigma, 
                            nrep=100)
        ## From LS with richness from posterior credible intervals
        ls.s2 <- NB.samp(rad = abc.ls.c.rad.l, 
                         tot.area = Tot.A, 
                         n.plots = N.plots,
                         lmean.sd =
                             predict(lm.sd, 
                                     newdata=data.frame(population=abc.ls.c.rad.l)),
                         lsd.sd = lm.sd.sigma,
                         lmean.k =  
                             predict(lm.k, 
                                     newdata=data.frame(dens.ha=abc.ls.c.rad.l/Tot.A)),
                         lsd.k = lm.k.sigma, 
                         nrep=100)
        ls.s3 <- NB.samp(rad = abc.ls.c.rad.u, 
                         tot.area = Tot.A, 
                         n.plots = N.plots,
                         lmean.sd =
                             predict(lm.sd, 
                                     newdata=data.frame(population=abc.ls.u.rad.l)),
                         lsd.sd = lm.sd.sigma,
                         lmean.k =  
                             predict(lm.k, 
                                     newdata=data.frame(dens.ha=abc.ls.c.rad.u/Tot.A)),
                         lsd.k = lm.k.sigma, 
                         nrep=100)
