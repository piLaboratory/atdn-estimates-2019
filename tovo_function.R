############################################################################################################
###################     UPSCALING SPECIES RICHNESS AND ABUNDANCES IN TROPICAL FOREST     ###################
############################################################################################################
#
# INPUT:  abundance_vec -> numeric vector containing species' abundances
#         p -> fraction between local area and global area
#         r.estimate -> estimation of the r parameter at the local scale
#         csi.estimate -> estimation of the csi parameter at the local scale

UpscalingNB<-function(abundance_vec,p,r.estimate,csi.estimate)
  
{ print("***********************************************************",quote=FALSE)
  print("************ UPSCALING SPECIES RICHNESS AND ABUNDANCES IN TROPICAL FOREST *************",quote=FALSE)
  print("***********************************************************",quote=FALSE)
  
  require(VGAM)
  require(boot)
  require(TeachingDemos)
  
  # check validity of parameters first estimates
  if(r.estimate<0) stop("r parameter must be positive")
  if(csi.estimate<0||csi.estimate>=1) stop("csi parameter must be non negative and less than one")
  
  print(" ",quote=FALSE)
  print("INITIAL DATA:",quote=FALSE)
  print(" ",quote=FALSE)
 
  # Total number of species 
  S.local<-length(which(abundance_vec!=0))
  print(c("Total number of found species:",S.local),quote=FALSE)
  
  # Total number of individuals
  print(c("Total number of found individuals:",sum(abundance_vec)),quote=FALSE)
  print(" ",quote=FALSE)
  print("***********************************************************",quote=FALSE)
  print(" ",quote=FALSE)
  
  print("FITTING RSA...",quote=FALSE)
  print(" ",quote=FALSE)
  
  # Fit of the empirical RSA at the local scale with a negative binomial distribution
  # normalized from 0 to infinity
  y<-as.numeric(abundance_vec)
  pdata <- data.frame(munb = r.estimate*csi.estimate/(1-csi.estimate), size = r.estimate)
  pdata <- transform(pdata, y = y)

  options(warn=-1)
  fit <- tryCatch(vglm(y ~ 1, posnegbinomial, pdata, trace = FALSE,eps.trig = 1e-8,max.support = Inf),error = function(e) {})
  if(is.null(fit)) stop("RSA fitting failed. Try with different r and csi initial estimates")
  options(warn=0)
  
  coefficients<-coef(fit, matrix = TRUE)
  
  mu.fit<-exp(coefficients[1])
  r.fit<-exp(coefficients[2])
  csi.fit<-mu.fit/(mu.fit+r.fit)
  
  print("RESULTS:",quote=FALSE)
  print(c("fitted r parameter:",signif(r.fit,4)),quote=FALSE)
  print(c("fitted csi parameter:",signif(csi.fit,4)),quote=FALSE)
  print(" ",quote=FALSE)
   
  # Standard errors for the fitted parameters
  options(warn=-1)
  std.error.logemu<-coef(summary(fit))[3]
  std.error.loger<-coef(summary(fit))[4]
  options(warn=0)
  
  std.error.mu<-mu.fit*std.error.logemu
  std.error.r<-r.fit*std.error.loger
  std.error.csi<-1/((r.fit+mu.fit)^2)*sqrt(r.fit^2*std.error.mu^2+mu.fit^2*std.error.r^2)
  
  print("ERRORS:",quote=FALSE)
  print(c("error on r:",signif(std.error.r,4)),quote=FALSE)
  print(c("error on csi:",signif(std.error.csi,4)),quote=FALSE)  
  print(" ",quote=FALSE)
  # Extrapolation of the csi parameter at the global scale
  csi.global<-csi.fit/(p+csi.fit*(1-p))
    
  # Estimation of total biodiversity at the global scale
  S.estimate<-S.local*(1-(1-csi.global)^r.fit)/(1-(1-csi.fit)^r.fit)
  
  # Computing the error on S.estimate
  upscaling<-expression( S.local*(1-(1-csi.fit/(p+csi.fit*(1-p)))^r.fit)/(1-(1-csi.fit)^r.fit) )
  
  der.r<-deriv(upscaling,"r.p")
  der.csi<-deriv(upscaling,"csi.p")
  
  std.error.S<-sqrt(as.numeric(eval(der.csi))^2*std.error.csi^2+as.numeric(eval(der.r))^2*std.error.r^2)
  
  print("***********************************************************",quote=FALSE)
  print(" ",quote=FALSE)
  print("UPSCALING RESULTS:",quote=FALSE)
  print(" ",quote=FALSE)
  print(c("Estimated total biodiversity S:",round(S.estimate)),quote=FALSE)
  print(c("Error on S:",round(std.error.S)),quote=FALSE)
  print(" ",quote=FALSE)
  print("***********************************************************",quote=FALSE)
  print(" ",quote=FALSE)
  print("Printing local and predicted SAD..",quote=FALSE)
  
  
  # Comparison beqtween local and global SAD
  
  # Empirical SAD at the local scale
  occupation_lognumbers<-log(abundance_vec,base=2)
  
  Preston_class<-0:ceiling(max(occupation_lognumbers))
  Preston_counts<-rep(0,length(Preston_class))
  Preston_counts[1]<-length(occupation_lognumbers[occupation_lognumbers==0])/2
  
  for (i in 2:length(Preston_class))
  {low<-length(occupation_lognumbers[occupation_lognumbers==Preston_class[i-1]])/2
  middle<-length(occupation_lognumbers[occupation_lognumbers<Preston_class[i] & occupation_lognumbers>Preston_class[i-1]])
  high<-length(occupation_lognumbers[occupation_lognumbers==Preston_class[i]])/2
  Preston_counts[i]<-low+middle+high
  }
  
  
  # Theoretical SAD at the local scale
  xl<-c(1:2^ceiling(log(max(abundance_vec),base=2)))
  yl<-rep(0,length(xl))
  for(i in 1:length(xl))
    {yl[i]<-S.local*(1-csi.fit)^r.fit*csi.fit^xl[i]*choose(xl[i]+r.fit-1,xl[i])/(1-(1-csi.fit)^r.fit)}
  
  occupation_lognumbers.local<-log(xl,base=2)
  
  Preston_class.local<-0:ceiling(max(occupation_lognumbers.local))
  Preston_counts.local<-rep(0,length(Preston_class.local))
  Preston_counts.local[1]<-yl[1]/2
  Preston_counts.local[2]<-yl[1]/2+yl[2]/2
  
  for (i in 3:length(Preston_class.local))
  {low<-yl[2^(Preston_class.local[i-1])]/2
  middle<-sum(yl[(2^(Preston_class.local[i-1])+1):(2^(Preston_class.local[i])-1)])
  high<-yl[2^Preston_class.local[i]]/2
  Preston_counts.local[i]<-low+middle+high
  }
  
  
  # Theoretical SAD at the global scale
  xg<-c(1:2^20)
  yg<-rep(0,length(xg))
  for(i in 1:length(xg))
  {yg[i]<-S.estimate*(1-csi.global)^r.fit*csi.global^xg[i]*choose(xg[i]+r.fit-1,xg[i])/(1-(1-csi.global)^r.fit)}
  
  occupation_lognumbers.global<-log(xg,base=2)
  
  Preston_class.global<-0:ceiling(max(occupation_lognumbers.global))
  Preston_counts.global<-rep(0,length(Preston_class.global))
  Preston_counts.global[1]<-yg[1]/2
  Preston_counts.global[2]<-yg[1]/2+yg[2]/2
  
  for (i in 3:length(Preston_class.global))
  {low<-yg[2^(Preston_class.global[i-1])]/2
  middle<-sum(yg[(2^(Preston_class.global[i-1])+1):(2^(Preston_class.global[i])-1)])
  high<-yg[2^Preston_class.global[i]]/2
  Preston_counts.global[i]<-low+middle+high
  }
  
  # Graphics of the SAD histograms
  xbin<-1.2
  xorigin<-0.7
  xlast<-length(Preston_class.global)*xbin
  xpos<-seq(xorigin,xlast,by=xbin)
  
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  histSAD<-barplot(c(Preston_counts,rep(0,length(Preston_class.global)-length(Preston_class)-1)),xlim=c(0,xlast),ylim=c(0,max(max(Preston_counts),max(Preston_counts.local),max(Preston_counts.global))),col=rgb(0.1,0.1,0.1,0.1),main="SAD",xlab=expression("Preston Class"),ylab="Number of species",mgp=c(1.8,0.5,0),mar=c(5.1, 4.1, 4.1, 8.1))
  axis(1,at=xpos,label=1:(length(xpos)),mgp=c(1.8,0.5,0))
  points(xpos,c(Preston_counts.local,rep(0,length(Preston_class.global)-length(Preston_class))),col="red",type="l",lwd=3)
  points(xpos,Preston_counts.global,col="blue",pch=19,type="l",lwd=3)
  legend("topright",inset=c(-0.55,0.2),legend=c("local SAD","global SAD"),lwd=c(3,3),col=c("red","blue"),seg.len=0.6)
}
