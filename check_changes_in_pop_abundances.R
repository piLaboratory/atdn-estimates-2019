pop.2013 <- read.csv2("2013_original.csv", as.is=TRUE)
pop.2013$dens <- pop.2013$N.ind/1048
pop.2013t <- read.csv2("2013_tax_2019.csv", as.is=TRUE)
pop.2013t$dens <- pop.2013t$N.ind/1048
pop.2019 <- read.csv2("2019_May.csv", as.is=TRUE)
pop.2019$dens <- pop.2019$N.ind/2042

pop.all <- merge(
    merge(pop.2013,pop.2013t, by="species", suffixes = c(".13", ".13t")),
    pop.2019, by= "species")
names(pop.all)[12:16] <- paste(names(pop.all)[12:16],"19", sep=".")


my_line <- function(x,y,...){
    points(x,y,...)
    abline(a = 0,b = 1, col="blue")
}

pdf("check_may_data%1d.pdf", onefile=FALSE)
pairs(pop.all[,c(2,7,12)], lower.panel=my_line, upper.panel=NULL, log="xy", col="lightgrey", cex=0.25)
pairs(pop.all[c(6,11,16)], lower.panel=my_line, upper.panel=NULL, log="xy", col="lightgrey", cex=0.25)
pairs(pop.all[,c(5,10,15)], lower.panel=my_line, upper.panel=NULL, log="xy", col="lightgrey", cex=0.25)
dev.off()
