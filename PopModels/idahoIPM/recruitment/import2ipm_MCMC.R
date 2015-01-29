# Import and format recruitment parameters
# then define recruitment function

# recruitment parameters
Rpars=list(intcpt.mu=rep(0,Nspp),intcpt.yr=matrix(0,Nyrs,Nspp),intcpt.tau=rep(100,Nspp),
           intcpt.gr=matrix(NA,6,Nspp),g.tau=rep(NA,Nspp),
           dd=matrix(NA,Nspp,Nspp),theta=rep(NA,Nspp),sizeMean=rep(NA,Nspp),sizeVar=rep(NA,Nspp),
           recSizes=list(1))

#Get random row selection from the MCMC chain
mcDraw <- sample(seq(1,2000,1), size = 1)
infile=paste("recruitment/recruitmentParamsMCMC_AllSpp.rds")
Rdata <- readRDS(infile)
Rdata <- as.data.frame(Rdata[mcDraw,]) #get one row from the MCMC chain


Rdata$Coef <- c(rep("dd", 4*Nspp),
                rep("tauGroup", Nspp),
                rep("Group", 6*Nspp),
                rep("Intercept", Nspp),
                rep("Intercept.tau", Nspp),
                rep("Intercept.yr", Nyrs*Nspp),
                rep("theta", Nspp),
                rep("u", Nspp))
Rdata$Species <- c(rep(sppList, each=Nspp),
                   rep(sppList, times=1),
                   rep(sppList, each=6),
                   rep(sppList, times=1),
                   rep(sppList, times=1),
                   rep(sppList, each=Nyrs),
                   rep(sppList, times=1),
                   rep(sppList, times=1))
  
colnames(Rdata)[1] <- "value"

for(i in 1:Nspp){
  infile=paste("../speciesData/",sppList[i],"/recSize.csv",sep="")
  recSize=read.csv(infile)
  Rpars$sizeMean[i]=mean(log(recSize$area))
  Rpars$sizeVar[i]=var(log(recSize$area))
  #Rpars$recSizes[[i]]=recSize$area
}

#needs to be like this: c[i,j] = effect of j on i
for(i in 1:Nspp){
  di <- which(Rdata$Coef=="dd")
  tmpDD <- Rdata[di,]
  di2 <- which(tmpDD$Species==sppList[i])
  Rpars$dd[i,] <- tmpDD$value[di2]
}

Rpars$intcpt.mu <- Rdata$value[which(Rdata$Coef=="Intercept")]

for(i in 1:Nspp){
  tmpInt <- Rdata[which(Rdata$Coef=="Intercept.yr"),]
  Rpars$intcpt.yr[,i] <- tmpInt$value[which(tmpInt$Species==sppList[i])]
}
Rpars$intcpt.tau <- Rdata$value[which(Rdata$Coef=="Intercept.tau")]

for(i in 1:Nspp){
  tmpG <- Rdata[which(Rdata$Coef=="Group"),]
  Rpars$intcpt.gr[,i] <- tmpG$value[which(tmpG$Species==sppList[i])]
}
Rpars$g.tau <- Rdata$value[which(Rdata$Coef=="tauGroup")]
Rpars$theta <- Rdata$value[which(Rdata$Coef=="theta")] 

rm(Rdata)

# define recruitment function
#number of recruits per area produced 
# cover is stored in absolute area (cm^2)
get.rpa=function(Rpars,cover,doYear){
  # cover is in m^2 per m^2; convert to % scale:
  cover2=cover*100
  # calculate recruits
  Nspp=length(cover)
  mu=rep(NA,Nspp)
  for(i in 1:Nspp){
    mu[i]=cover2[i]*exp(Rpars$intcpt.yr[doYear,i]+sqrt(cover2)%*%Rpars$dd[i,]) 
  }
  if(sum(is.na(mu))>0) browser() # stop for errors
  rpa=mu/(cover*A)  # convert from number recruits to recruits per cm^2
  return(rpa)
}

# Fecundity function, expected number of recruits of size y produced by a size x individual
# The size distribution of recruits is on the log scale
f=function(v,u,Rpars,rpa,doSpp) { 
  nRecruits = rpa[doSpp]*exp(u)
  #probability of producing a seedling of size v
  tmp=dnorm(v,Rpars$sizeMean[doSpp],sqrt(Rpars$sizeVar[doSpp]))/(1-pnorm(-1.61,Rpars$sizeMean[doSpp],sqrt(Rpars$sizeVar[doSpp])))
  #number recruits of each size 
  f=nRecruits*tmp
  return(f)
}   