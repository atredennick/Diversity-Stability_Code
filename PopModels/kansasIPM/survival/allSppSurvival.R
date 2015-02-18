#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

sppList=c("BOCU","BOHI","SCSC") # for spp in alphabetical order
alpha.effect=c(0.006,0.012,0.007) # from Chengjin's code

for(spp in 1:length(sppList)){
# for(spp in 2:2){
  doSpp=sppList[spp]
  outfile=paste("Surv_params_",doSpp,".csv",sep="")
  
  growDfile=paste("../../../Data/Kansas/",doSpp,"/survD.csv",sep="")
  growD=read.csv(growDfile)
  ##excluding the records later than 1968, to keep the same random year effect...
  ######################################################################
  ##to remove some points:
  #for q25
  tmp1<-which(growD$quad=="q25" & (growD$year<35 | growD$year>62))
  #for q27
  tmp2<-which(growD$quad=="q27")
  #for q28
  tmp3<-which(growD$quad=="q28")
  #for q30
  tmp4<-which(growD$quad=="q30")
  #for q31
  tmp5<-which(growD$quad=="q31" & (growD$year<35 | growD$year>39))
  #for q32
  tmp6<-which(growD$quad=="q32" & (growD$year<35 | growD$year>41))
  tmp<-c(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)
  growD<-growD[-tmp,]
  growD<-subset(growD,year<68)
#   growD$Group=as.factor(substr(growD$quad,1,1)) ##add Group information
  
  D=growD  #subset(growD,allEdge==0)
  D$logarea=log(D$area)
  D$quad=as.character(D$quad)
  D$year=as.factor(D$year)
  
  
  ##then we moved some specific points:
#   tmp2<-which(D$quad=="A12" & D$year==44)
#   tmp3<-which(D$quad=="B1"  & D$year==44)
#   tmp41<-which(D$quad=="E4" & D$year==33) 
#   tmp42<-which(D$quad=="E4" & D$year==34) 
#   tmp43<-which(D$quad=="E4" & D$year==43)
#   tmp44<-which(D$quad=="E4" & D$year==44)
#   
#   tmpONE<-c(tmp2,tmp3,tmp41,tmp42,tmp43,tmp44)
#   if(length(tmpONE)>0) D<-D[-tmpONE,]

  # calculate crowding 
  for(i in 1:length(sppList)){
    distDfile=paste("../../../Data/Kansas/",sppList[i],"/",sppList[i],"_genet_xy.csv",sep="")
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

  library(INLA)
  #Set up ID variables for INLA random effects
  D$yearID <- D$year #for random year offset on intercept
  D$GroupID <- as.numeric(D$Group)
  ngroup <- length(unique(D$Group))
  
  #Instead of full model, match the structure of the quadrat-based IBM regressions
  formula2 <- survives ~ logarea+crowd+
    f(yearID, model="iid", prior="normal",param=c(1,0.001))+
    f(GroupID, model="iid", prior="normal",param=c(0,0.001))+
    f(year, logarea, model="iid", prior="normal",param=c(0,0.001))
  
  outINLA <- inla(formula2, data=D,
                  family=c("binomial"), verbose=TRUE,
                  control.compute=list(dic=T,mlik=T),
                  control.predictor = list(link = 1),
                  control.inla = list(h = 1e-6),
                  Ntrials=rep(1,nrow(D)))
  summary(outINLA)
      
  #Collect parameters
  #random year and group effects
  params <- as.data.frame(outINLA$summary.random$yearID[,1:2])
  params <- cbind(params, outINLA$summary.random$year[,2])
  params <- cbind(params, c(outINLA$summary.random$GroupID[,2],rep(NA,nrow(params)-ngroup)))
  names(params) <- c("Year", "Intercept.yr", "logarea.yr", "Group")
  #fixed effects
  fixed <- as.data.frame(outINLA$summary.fixed)[,1:2]
  tmp=matrix(NA,dim(params)[1],nrow(fixed))
  colnames(tmp)=rownames(fixed)
  tmp[1,]=fixed[,1]
  params=cbind(params,tmp)
  params$alpha=NA; params$alpha[1:length(sppList)]=alpha.effect
  colnames(params)[5] <- "Intercept"

  # write output
  write.table(params,outfile,row.names=F,sep=",")
}
