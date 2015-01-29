#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))
library(rjags)
sppList=sort(c("PSSP","HECO","POSE","ARTR"))
alpha.effect=c(0.014,0.017,0.026,0.018) # for spp in alphabetical order

for(spp in 1:length(sppList)){
  doSpp=sppList[spp]
  outfile=paste("Growth_paramsMCMCind_",doSpp,".rds",sep="")
  outfile2=paste("growthStatsInd_",doSpp,".csv",sep="")
  outfile3=paste("growthQuantsInd_",doSpp,".csv",sep="")
  growDfile=paste("../../speciesData/",doSpp,"/growDnoNA.csv",sep="")
  growD=read.csv(growDfile)
  
  D=growD  #subset(growD,allEdge==0)
  D$logarea.t0=log(D$area.t0)
  D$logarea.t1=log(D$area.t1)
  D$quad=as.character(D$quad)
  
  ##then we moved some specific points:
  ##remove outliers (large plants that obviously do not turn into tiny plants) for ARTR only
  if(doSpp=="ARTR"){
    tmp=which(D$quad=="Q23" & D$year==45 & D$trackID==67)
    tmp=c(tmp,which(D$quad=="Q12" & D$year==55 & D$trackID==25))
    tmp=c(tmp,which(D$quad=="Q26" & D$year==45 & D$trackID==73))
    D=D[-tmp,]
  }else{
    D=D
  }

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
  
  #crowd=W[,which(sppList==doSpp)] #for single species
  crowd = W #for multispecies
  crowd[crowd<1e-99]=0 ###this is very key...
  
  D$plantID <- paste(D$quad,D$trackID,sep = "")
  plant <- as.numeric(as.factor(D$plantID))
  #Set up data for JAGS
  Nyears = length(unique(D$year))
  Ngroups = length(unique(D$Group))
  Nspp = length(sppList)
  Nplants = length(unique(plant))
  dataJ <- list(yr = as.numeric(as.factor(D$year)),
                Nyears = Nyears,
                grp = as.numeric(D$Group),
                Ngroups = Ngroups,
                Nspp = Nspp,
                crowd = crowd,
                x = log(D$area.t0),
                y = log(D$area.t1),
                Nobs = nrow(D),
                Nplants = Nplants,
                plant = plant)
  params=c("intcpt.mu","intcpt.tau","intcpt.yr","intG",
           "slope.mu","slope.tau","slope.yr",
           "NBbeta.mu","NBbeta.tau","tauGroup",
           "tau","tauSize", "intI")
  
  #Run MCMC through JAGS
  inits=list(1)
  inits[[1]]=list(
    intcpt.mu=0,
    intcpt.tau=1,
    intcpt.yr=rep(0,Nyears),
    intG=rep(0.1,Ngroups),
    slope.mu=0.5,
    slope.tau=1,
    slope.yr=rep(0.1,Nyears),
    NBbeta.mu=rep(-0.2,Nspp),
    NBbeta.tau=0.2,
    tauGroup=2,  
    tau=1,
    tauSize=1
  )
  
  inits[[2]]=list(
    intcpt.mu=-1,
    intcpt.tau=0.1,
    intcpt.yr=rep(-1,Nyears),
    intG=rep(-0.5,Ngroups),
    slope.mu=1,
    slope.tau=0.1,
    slope.yr=rep(0.5,Nyears),
    NBbeta.mu=rep(-1,Nspp),
    NBbeta.tau=0.1,
    tauGroup=1,
    tau=0.1,
    tauSize=0.1
  )
  iterations <- 5000
  adapt <- 2000
  mod <- jags.model("growth_JAGS_ind.R", data=dataJ, n.chains=2, inits=inits, n.adapt=adapt)
  update(mod, n.iter = (iterations))
  out <- coda.samples(mod, params, n.iter=iterations, n.thin=10)
  outStats <- summary(out)$stat
  outQuants <- summary(out)$quantile
  
  outMCMC <- rbind(out[[1]][(iterations-999):iterations,], 
                   out[[2]][(iterations-999):iterations,])
  saveRDS(object = outMCMC, file = outfile)
  write.csv(outStats, outfile2)
  write.csv(outQuants, outfile3)
}#next species for MCMC
