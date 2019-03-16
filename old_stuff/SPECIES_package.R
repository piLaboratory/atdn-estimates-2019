## Using all estimators in SPECIES package
library(SPECIES)
Ab <- as.data.frame(table(dados$N.ind))
names(Ab) <- c("j", "n_j")
Ab$j <- as.integer(Ab$j)
##jackknife method
(Sj <- jackknife(butterfly,k=5))
##using only'ACE coverage method
(SChao92 <- ChaoLee1992(Ab,t=10, method="all"))
##using chao1984 lower bound estimator
(SChao84 <- chao1984(Ab))
##using Chao and Bunge coverage-duplication method
(SChaoB <- ChaoBunge(Ab,t=10))
##penalized NPMLE method
(SNPMLE <- pnpmle(Ab,t=15, C=1))
##unconditonal NPMLE method
(SUNPMLE <-unpmle(Ab,t=10, C=1))
##Poisson-compound Gamma method
(SPG <- pcg(Ab,t=20, C=1))
