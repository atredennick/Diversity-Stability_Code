#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

sppList=sort(c("PSSP","HECO","POSE","ARTR"))
alpha.effect=c(0.014,0.017,0.026,0.018) # for spp in alphabetical order

# fixedEffects <- list()
# SSRs <- matrix(ncol=5, nrow=(length(sppList)))

for(spp in 1:length(sppList)){
  doSpp=sppList[spp]
#   climD=read.csv("../../../weather/Climate.csv")
  
  outfile=paste("Growth_params_",doSpp,".csv",sep="")
  
  growDfile=paste("../../../../Data/Idaho/",doSpp,"/growDnoNA.csv",sep="")
  growD=read.csv(growDfile)
#   growD$Group=as.factor(substr(growD$quad,1,1)) ##add Group information
  
  D=growD  #subset(growD,allEdge==0)
  D$logarea.t0=log(D$area.t0)
  D$logarea.t1=log(D$area.t1)
  D$quad=as.character(D$quad)
#   climD$year=climD$year-1900
#   D=merge(D,climD)
#   D$year=as.factor(D$year)
  
  
  ##then we moved some specific points:

#   # remove outliers (large plants that obviously do not turn into tiny plants) for ARTR only
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
    distDfile=paste("../../../../Data/Idaho/",sppList[i],"/",sppList[i],"_genet_xy.csv",sep="")
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
  
#   crowd=W[,which(sppList==doSpp)] #for single species
crowd = W #for multispecies
crowd[crowd<1e-99]=0 ###this is very key...
# colnames(crowd) <- c("crowd1", "crowd2", "crowd3", "crowd4")

  library(INLA)
  #Set up ID variables for INLA random effects
  D$yearID <- D$year+max(D$year) #for random year offset on intercept
  D$GroupID <- as.numeric(D$Group)
#   D$age2 <- D$age
#   D <- cbind(D, crowd)
  
  #Instead of full model, match the structure of the quadrat-based IBM regressions
  formula2 <- logarea.t1 ~ logarea.t0+crowd+
    f(yearID, model="iid", prior="normal",param=c(0,0.001))+
    f(GroupID, model="iid", prior="normal",param=c(0,0.001))+
    f(year, logarea.t0, model="iid", prior="normal",param=c(0,0.001))
  
  outINLA <- inla(formula2, data=D,
                  family=c("gaussian"), verbose=TRUE,
                  control.predictor = list(link = 1),
                  control.compute=list(dic=T,mlik=T),
                  control.inla = list(h = 1e-6))
  summary(outINLA)

library(lme4)
# # model 7 is best, but behaves poorly, use model 5
out=lmer(logarea.t1~logarea.t0+W+(1|Group)+(logarea.t0|year),data=D)
summary(out)
AIC(out)

  #fit variance
  x = outINLA$summary.fitted.values$mean #fitted values from INLA
  y = (D$logarea.t1-outINLA$summary.fitted.values$mean)^2 #calculates the variance (residuals^2)
  plot(x,y)
  outVar=nls(y~a*exp(b*x),start=list(a=1,b=0)) #model the variance as function of fitted size
  
  #Collect parameters
  #random year and group effects
  params <- as.data.frame(outINLA$summary.random$yearID[,1:2])
  params <- cbind(params, outINLA$summary.random$year[,2])
  names(params) <- c("Year", "Intercept.yr", "logarea.t0.yr")
  tmp <- as.data.frame(outINLA$summary.random$Group[,2])
  names(tmp) <- "Group"
  tmp$GrpName <- outINLA$summary.random$Group[,1]
  tmp[(NROW(tmp)+1):NROW(params),]=NA
  params=cbind(params,tmp)
  #fixed effects
  fixed <- as.data.frame(outINLA$summary.fixed)[,1:2]
  tmp=matrix(NA,dim(params)[1],nrow(fixed))
  colnames(tmp)=rownames(fixed)
  tmp[1,]=fixed[,1]
  params=cbind(params,tmp)
  colnames(params)[6] <- "Intercept"
  params$alpha=NA; params$alpha[1:length(sppList)]=alpha.effect
  #variance 
  params$sigma.a=NA; params$sigma.a[1]=coef(outVar)[1] 
  params$sigma.b=NA; params$sigma.b[1]=coef(outVar)[2]
  
  # write output
  write.table(params,outfile,row.names=F,sep=",")
}
