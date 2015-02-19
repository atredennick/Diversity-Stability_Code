# Import and format growth parameters
# then define growth function

# growth parameters
Gpars=list(intcpt=rep(NA,Nspp),intcpt.yr=matrix(0,Nyrs,Nspp), intcpt.gr=matrix(0,6,Nspp),
  slope=rep(NA,Nspp),slope.yr=matrix(0,Nyrs,Nspp),
  nb=matrix(0,Nspp,Nspp),alpha=matrix(NA,Nspp,Nspp),
  sigma2.a=rep(NA,Nspp),sigma2.b=rep(NA,Nspp))
for(i in 1:Nspp){
  infile=paste("growth/cache/Growth_params_",sppList[i],".csv",sep="")
  Gdata=read.csv(infile)
  Gpars$intcpt[i]=Gdata$Intercept[1]
  tmp=which(names(Gdata)=="Group")
  if(length(tmp)>0) Gpars$intcpt.gr[,i]=Gdata$Group[!is.na(Gdata$Group)] 
  Gpars$intcpt.yr[,i]=Gdata$Intercept.yr
  Gpars$slope[i]=Gdata$logarea.t0[1]
  # random effects on slope
  tmp=which(names(Gdata)=="logarea.t0.yr")
  if(length(tmp)>0) Gpars$slope.yr[,i]=Gdata[,tmp]
  # get competition coefficients
  tmp=paste("crowd",1:length(sppList),sep="")
  tmp=which(is.element(names(Gdata),tmp))
  if(length(tmp)>0) Gpars$nb[i,]=as.numeric(Gdata[1,tmp])
  # get competition X size interactions coefficients
#   tmp=paste("logarea.t0.W",1:length(sppList),sep="")
#   tmp=which(is.element(names(Gdata),tmp))
#   if(length(tmp)>0) Gpars$slopeXnb[i,]=as.numeric(Gdata[1,tmp])   
  # get yr random effects on competition
#   tmp=paste("W",1:length(sppList),".yr",sep="")
#   tmp=which(is.element(names(Gdata),tmp))
#   if(length(tmp)>0) Gpars$nb.yr[i,,]=as.matrix(Gdata[,tmp])
  Gpars$alpha[i,]=Gdata$alpha[1:length(sppList)]
  Gpars$sigma2.a[i]=Gdata$sigma.a[1]
  Gpars$sigma2.b[i]=Gdata$sigma.b[1]
} # next i
rm(Gdata)

# growth function
G=function(v,u,W,Gpars,doYear,doSpp){
  mu=Gpars$intcpt[doSpp]+Gpars$intcpt.yr[doYear,doSpp]+(Gpars$slope[doSpp]+Gpars$slope.yr[doYear,doSpp])*u+
    W%*%(Gpars$nb[doSpp,])
  sigma2=Gpars$sigma2.a[doSpp]*exp(Gpars$sigma2.b[doSpp]*mu)
  out=dnorm(v,mu,sqrt(sigma2))
  out
}

