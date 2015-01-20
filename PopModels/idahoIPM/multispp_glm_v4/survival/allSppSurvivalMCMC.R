#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

library(rjags)

sppList=c("ARTR","HECO","POSE","PSSP")
alpha.effect=c(0.004,0.048,0.040,0.017) # for spp in alphabetical order

for(spp in 1:length(sppList)){
  # for(spp in 2:2){
  doSpp=sppList[spp]
  outfile=paste("survParamsMCMC_",doSpp,".rds",sep="")
  
  growDfile=paste("../../speciesData/",doSpp,"/survD.csv",sep="")
  growD=read.csv(growDfile)
  #   growD$Group=as.factor(substr(growD$quad,1,1)) ##add Group information
  
  D=growD  #subset(growD,allEdge==0)
  D$logarea=log(D$area)
  D$quad=as.character(D$quad)
  D$year=as.factor(D$year)
  
  # calculate crowding 
  for(i in 1:length(sppList)){
    distDfile=paste("../../speciesData/",sppList[i],"/",sppList[i],"_genet_xy.csv",sep="")
    if(i==1){
      distD=read.csv(distDfile)
      distD$nbSpp=sppList[i]  
    }else{
      tmp=read.csv(distDfile)
      tmp$nbSpp=sppList[i] 
      distD=rbind(distD,tmp)
    }
  }
  
  distD=distD[,c("quad","year","trackID","area","nbSpp","x","y")]
  W=matrix(NA,dim(D)[1],length(sppList))
  for(i in 1:dim(D)[1]){
    tmpD=subset(distD,year==D$year[i] & quad==D$quad[i])
    focal=which(tmpD$trackID==D$trackID[i] & tmpD$nbSpp==doSpp)
    xx=tmpD$x[focal] ; yy=tmpD$y[focal]
    tmpD$distance=sqrt((xx-tmpD$x)^2+(yy-tmpD$y)^2)
    tmpD=subset(tmpD,distance>0)
    if(dim(tmpD)[1]>0){
      for(k in 1:length(sppList)){
        sppI=which(tmpD$nbSpp==sppList[k])
        if(length(sppI)>0){
          W[i,k]=sum(exp(-1*alpha.effect[k]*tmpD$distance[sppI]^2)*tmpD$area[sppI])         
        }else{
          W[i,k]=0
        }
      }
    }else{
      W[i,]=0
    }   
  }
  
  crowd=W
  
  #Set up data for JAGS
  Nyears = length(unique(D$year))
  Ngroups = length(unique(D$Group))
  Nspp = length(sppList)
  dataJ <- list(yr = as.numeric(as.factor(D$year)),
                Nyears = Nyears,
                grp = as.numeric(D$Group),
                Ngroups = Ngroups,
                Nspp = Nspp,
                crowd = crowd,
                x = log(D$area),
                y = D$survives,
                Nobs = nrow(D))
  params=c("intcpt.mu","intcpt.tau","intcpt.yr","intG",
           "slope.mu","slope.tau","slope.yr",
           "NBbeta.mu","NBbeta.tau","tauGroup",
           "tau","tauSize")
  
  inits=list(1)
  inits[[1]]=list(
    intcpt.mu=0,
    intcpt.tau=1,
    intcpt.yr=rep(0,Nyears),
    intG=rep(0,Ngroups),  
    slope.mu=0.5,
    slope.tau=1,
    slope.yr=rep(0,Nyears),
    NBbeta.mu=rep(0,Nspp),
    tauGroup=1,
    NBbeta.tau=0.2
  )
  
  inits[[2]]=list(
    intcpt.mu=1,
    intcpt.tau=10,
    intcpt.yr=rep(1,Nyears),
    intG=rep(1,Ngroups),  
    slope.mu=1,
    slope.tau=10,
    slope.yr=rep(1,Nyears),
    NBbeta.mu=rep(1,Nspp),
    tauGroup=1,
    NBbeta.tau=0.3
  )
  
  iterations <- 5000
  adapt <- 2000
  mod <- jags.model("survival_JAGS.R", data=dataJ, n.chains=2, inits=inits, n.adapt=adapt)
  update(mod, n.iter = (iterations))
  out <- coda.samples(mod, params, n.iter=iterations, n.thin=10)
  outStats <- summary(out)$stat
  outQuants <- summary(out)$quantile
  
  outMCMC <- rbind(out[[1]][(iterations-999):iterations,], 
                   out[[2]][(iterations-999):iterations,])
  saveRDS(object = outMCMC, file = outfile)
}#end species loop
