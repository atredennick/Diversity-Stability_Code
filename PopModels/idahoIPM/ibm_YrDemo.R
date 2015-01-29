################################################################################
# Below is an IBM for simulating long-term population dynamics with demographic
#   and environmental stochasticity. Parameters for vital rate functions come
#   from previously fit statistical models.
#
# email: atredenn@gmail.com

################################################################################
## Clear the workspace:
rm(list=ls(all=TRUE))

################################################################################
## Required libraries:
library(boot) #for inverse logit function (inv.logit)

################################################################################
## Some global variables:
nSims=5
nYears=20
nGrps=6
startYr=32
#qName="Q1"
#doGroup=6
L=100 # dimension of square quadrat (cm)
expand=1  # 1 = 1x1 m^2, 2 = 2x2m^2, etc
sppList=c("ARTR","HECO","POSE","PSSP")
minSize=0.25
maxSize=c(8000,500,500,500)
myCol=c("black","gold1","blue","red")

################################################################################
## Get vital rate parameters:


################################################################################
## Vital rate functions:
survive=function(Spars,doSpp,doYear,sizes,crowding){
  logsizes=log(sizes)
  mu=Spars$intcpt.yr[doYear,doSpp]+Spars$slope.yr[doYear,doSpp]*logsizes+Spars$nb[doSpp,]%*%crowding
  out=inv.logit(mu)
  out=rbinom(length(sizes),1,out)
  out
}

grow=function(Gpars,doSpp,doYear,sizes,crowding){
  # crowding and nb are vectors of dim Nspp
  logsizes=log(sizes)
  mu=Gpars$intcpt.yr[doYear,doSpp]+Gpars$slope.yr[doYear,doSpp]*logsizes+Gpars$nb[doSpp,]%*%crowding
  tmp=which(mu<log(minSize)*1.5)  # we will kill vanishingly small plants...below
  mu[tmp]=log(minSize) # truncate tiny sizes (they cause problems in sigma2)
  sigma2=Gpars$sigma2.a[doSpp]*exp(Gpars$sigma2.b[doSpp]*mu)
  out=exp(rnorm(length(sizes),mu,sqrt(sigma2)))
  if(sum(is.na(out))>0) browser()
  out[tmp]=0   # here's the killing of really small plants
  out[out>maxSize[doSpp]]=maxSize[doSpp] #truncate big plants for maximum observed size
  out
}

recruit=function(Rpars,sizes,spp,doYear,lastID,L,expand){
  # dd is a matrix of dim Nspp X Nspp
  # sizes and spp are vectors of the same length (=N plants)
  # calculate total areas
  totArea=aggregate(sizes,by=list("spp"=spp),FUN=sum)
  # put in missing zeros
  tmp=data.frame("spp"=1:length(sppList))
  totArea=merge(totArea,tmp,all.y=T)
  totArea[is.na(totArea)]=0
  totArea=totArea[,2]/((L*expand)^2)*100  # scale to % cover   
  # calculate recruits
  lambda=rep(NA,Nspp) # seed production
  for(i in 1:Nspp){
    lambda[i]=totArea[i]*exp(Rpars$intcpt.yr[doYear,i]+sqrt(totArea)%*%Rpars$dd[i,])
  }
  # number of draws from distribution depends on size of landscape
  NN=rnbinom(length(lambda)*expand^2,mu=lambda,size=Rpars$theta)  
  NN=rowSums(matrix(NN,length(lambda),expand^2))      
  x=y=spp=id=size=NULL
  for(i in 1:Nspp){
    if(NN[i]>0){
      #get recruit sizes 
      size=c(size,exp(rnorm(NN[i],Rpars$sizeMean[i],sqrt(Rpars$sizeVar[i]))))
      if(sum(is.na(size))>0) stop("Check recruit sizes")
      #assign random coordinates
      x=c(x,expand*L*runif(NN[i])); y=c(y,expand*L*runif(NN[i]))
      spp=c(spp,rep(i,NN[i]))
      #assign genet ID's
      if(length(id)>0) lastID=max(id)
      id=c(id,(1:NN[i]+lastID))
    }      
  } # next i
  # output
  size[size<minSize]=minSize
  out=cbind(spp,size,x,y,id)
  return(out)
}

# crowding function, assumes toroidal landscape
getCrowding=function(plants,alpha,L,expand){
  # plants is a matrix: sizes in column 1; x,y coords in columns 2 and 3
  # d is the distance weighting parameter
  # function returns a vector of length = rows in plants
  if(dim(plants)[1]>1){
    xdiff=abs(outer(plants[,3],plants[,3],FUN="-"))
    ydiff=abs(outer(plants[,4],plants[,4],FUN="-"))
    distMat=sqrt(xdiff^2+ydiff^2) 
    distMat[distMat==0]=NA
    #distMat[distMat<1]=1  # CLUDGE TO PREVENT GROWTH EXPLOSIONS
    distMat=exp(-1*alpha[plants[,1]]*distMat^2)
    sizeMat=matrix(plants[,2],dim(plants)[1],dim(plants)[1])
    out=aggregate(distMat*sizeMat,by=list("spp"=plants[,1]),FUN=sum,na.rm=T)
    # put in missing zeros
    tmp=data.frame("spp"=c(1:length(sppList)))
    out=merge(out,tmp,all.y=T)
    out[is.na(out)]=0
    out=out[order(out$spp),]
    out=as.matrix(out[,c(2:NCOL(out))])  # drop spp column
  }else{
    out=rep(0,Nspp)
  }
  out
}

################################################################################
## Function to simulate one time step:
simOne <- function(nSpp, doGroup, plants){
  doYr=sample(1,yrList) #draw year effects
  nextplants=plants
  ##recruitment
  newplants=recruit(Rpars,sizes=plants[,2],spp=plants[,1],doYear=doYr,L,expand)
  for(ss in 1:nSpp){ #loop through species for survival and growth
    if(N[tt,ss]>0){ #make sure spp ss is not extinct
      ##growth
      W=getCrowding(plants,Gpars$alpha[ss,],L,expand)
      newsizes=grow(Gpars,doSpp=ss,doYear=doYr,sizes=plants[,2],crowding=W)
      if(sum(newsizes==Inf)>0) browser()
      if(is.na(sum(newsizes))) browser()
      ##survival
      #uses same W as growth            
      live=survive(Spars,doSpp=ss,doYear=doYr,sizes=plants[,2],crowding=W)
      #put it all together
      tmp=which(plants[,1]==ss) #only alter plants of focal spp        
      nextplants[tmp,2]=newsizes[tmp]*live[tmp] #update with G*S
    } #end IF no plants
  } #next species (ss) 
  
  nextplants=nextplants[nextplants[,2]>0,] #remove dead plants 
  nextplants=rbind(nextplants,newplants) #add recruits
  return(nextplants)
}

################################################################################
## Get intitial conditions for model
Nspp=length(sppList)
init.plants=list(NULL)
lastID=rep(0,4)
for(i in 1:length(sppList)){
  infile=paste("../speciesData/",sppList[i],"/",sppList[i],"_genet_xy.csv",sep="")
  tmpD=read.csv(infile)
  tmpD=subset(tmpD, tmpD$year==min(tmpD$year))
  spp=rep(i,dim(tmpD)[1])
  tmpD2=data.frame(cbind(spp,tmpD[,c("area","x","y")]))
  names(tmpD2)=c("spp","size","x","y")
  init.plants=rbind(init.plants,tmpD2)
}

# plot initial conditions
par(mgp=c(2,0.5,0),tcl=-0.2)
symbols(x = init.plants[,3], y = init.plants[,4], circles = sqrt(init.plants[,2]/pi), fg=myCol[init.plants[,1]],
        xlim=c(0,L*expand),ylim=c(0,L*expand),main ="Time=1",xlab="x",ylab="y",inches=F,lwd=2)

#turn init.plants to cover
initArea <- aggregate(init.plants[,2],by=list(init.plants[,1]),FUN=sum)

################################################################################
## Run simulation
A <- matrix(0,nYears,Nspp) #storage matrix
A[1,1:4] <- initArea[,2] #initial conditions
doPlants <- init.plants
for(tt in 1:nYears){
  plantsNow <- simOne(nSpp = 4, doGroup = NA, plants = doPlants)
  A[tt,] <- aggregate(plantsNow[,2],by=list(plantsNow[,1]),FUN=sum)
  doPlants <- plantsNow
}

################################################################################
## Save output


