# Multispecies, spatially implicit IPM
# This version makes it possible to assume "no overlap"
# for intraspecific competition only or intra- and interspecific competition

#This version allows you to "grow" each species in isolation, to obtain
#the species' intrinsic growth rate (r in Lotka-Volterra terms)

#We do this in a fixed environment where all year effects are set to zero

####
#### This version is a sensitivity analysis, chaning Gpars, Rpars, and Spars intercepts
####

setwd(dir = "../")

# ATT 8/26/14
outfile1="ipm_cover_intrinsicGrowth.csv"

A=10000 #Area of 100cm x 100cm quadrat
tlimit=200 ## number of years to simulate
burn.in=100    # years to cut before calculations
sppList=c("ARTR","HECO","POSE","PSSP")
bigM=c(75,75,50,50)     #Set matrix dimension for each species
maxSize=c(3000,202,260,225)    # in cm^2: PSSP=225 HECO=202  POSE=260  ARTR=3000  # minSize=0.2  cm^2
Nyrs=22
doGroup=NA  # NA for spatial avg., values 1-6 for a specific group
constant=T  # constant environment
NoOverlap.Inter=T # no overlap of heterospecifics
compScale=F # not well implemented, but for rescaling competition coefficients

changeVec <- c(-0.25, 0, 0.25)
maxR <- matrix(NA, nrow=length(sppList), ncol=length(changeVec))

for(jjjj in 1:length(sppList)){
  #####################################
  fixSpp=sppList[jjjj] ##need to change if for other species
  fixCov=0.0001 # in % cover
  #####################################
  
  #============================================================
  # (I) LOAD VITAL RATE PARAMETERS & FUNCTIONS
  #============================================================
  Nspp=length(sppList)
  
  # set up survival parameters and function
  source("./survival/import2ipm_noOverlap.r")
  # set up growth parameters and function
  source("./growth/import2ipm_noOverlap.r")
  # set up recruitment parameters and function
  source("./recruitment/import2ipm.r")
  
  # get stable size distribution for fixed species
  infile=paste(fixSpp,"_stable_size.csv",sep="")
  sizeD=read.csv(infile)
  
  # model spatial group variation (or not)
  if(!is.na(doGroup)){
    Spars$intcpt=Spars$intcpt+Spars$intcpt.gr[doGroup,]
    Gpars$intcpt=Gpars$intcpt+Gpars$intcpt.gr[doGroup,]
    Rpars$intcpt.yr=Rpars$intcpt.yr+matrix(Rpars$intcpt.gr[doGroup,],Nyrs,Nspp,byrow=T)
  }
  
  # PERTURB PARAMETERS -------------------------------------
  if(constant==T){
    #turn off random year effects  
    Rpars$intcpt.yr=matrix(Rpars$intcpt.mu,Nyrs,Nspp,byrow=T)
    Gpars$intcpt.yr[]=0;Gpars$slope.yr[]=0
    Spars$intcpt.yr[]=0;Spars$slope.yr[]=0    
  }
  
  if(compScale==T){
    tmp <- matrix(rep(-0.1, length(Gpars$nb)), nrow = 4)
    for(i in 1:Nspp){
      tmp[i,i] <- Gpars$nb[i,i] 
    }
    Gpars$nb <- tmp
  }
  
  intcpSpp <- Gpars$intcpt[jjjj]
  
  for(doSens in 1:length(changeVec)){
    Gpars$intcpt[jjjj] <- intcpSpp + intcpSpp*changeVec[doSens]
    
    #============================================================================================#
    # (II) Simulation length, Matrix size and initial vectors
    #============================================================================================#
    
    v=v.r=b.r=expv=Cr=WmatG=WmatS=list(4)
    h=r.L=r.U=Ctot=numeric(4)
    for(i in 1:Nspp){
      
      # minimum (0.9*minimum size from data) and maximum sizes (1.1*maximum size from data)
      L=log(0.2)
      U=log(maxSize[i])*1.1     
      
      # boundary points b and mesh points y. Note: b chops up the size interval (L-U) into bigM-equal-sized portions.
      b = L+c(0:bigM[i])*(U-L)/bigM[i] 
      
      # v calculates the middle of each n-equal-sized portion.
      v[[i]] = 0.5*(b[1:bigM[i]]+b[2:(bigM[i]+1)])
      
      # step size for midpoint rule. (see equations 4 and 5 in Ellner and Rees (2006) Am Nat.)
      h[i] = v[[i]][2]-v[[i]][1]  
      
      # variables for Wr approximation
      b.r[[i]]=sqrt(exp(b)/pi)
      v.r[[i]]=sqrt(exp(v[[i]])/pi)
      expv[[i]]=exp(v[[i]])
      r.L[i] = sqrt(exp(L)/pi)
      r.U[i] = sqrt(exp(U)/pi)
      WmatG[[i]]=matrix(NA,length(v.r[[i]]),Nspp)  # storage of size-specific W values for each focal species
      WmatS[[i]]=matrix(NA,length(v.r[[i]]),Nspp)
      
      
    } # next species
    tmp=range(v.r)
    size.range=seq(tmp[1],tmp[2],length=50) # range across all possible sizes
    
    #============================================================================================#
    # (III) Utility functions
    #============================================================================================#
    
    # load the necessary libraries
    library(boot)
    library(mvtnorm)
    library(msm)
    library(statmod)  
    
    ## combined kernel
    make.K.values=function(v,u,muWG,muWS, #state variables
                           Rpars,rpa,Gpars,Spars,doYear,doSpp){  #growth arguments
      f(v,u,Rpars,rpa,doSpp)+S(u,muWS,Spars,doYear,doSpp)*G(v,u,muWG,Gpars,doYear,doSpp) 
    }
    
    # Function to make iteration matrix based only on mean crowding
    make.K.matrix=function(v,muWG,muWS,Rpars,rpa,Gpars,Spars,doYear,doSpp) {
      muWG=expandW(v,v,muWG)
      muWS=expandW(v,v,muWS)
      
      K.matrix=outer(v,v,make.K.values,muWG,muWS,Rpars,rpa,Gpars,Spars,doYear,doSpp)
      return(h[doSpp]*K.matrix)
    }
    
    # Function to format the W matrix for the outer product
    expandW=function(v,u,W){
      if(dim(W)[1]!=length(u)) stop("Check size of W")
      Nspp=dim(W)[2]
      W=as.vector(W)
      W=matrix(W,length(W),ncol=length(v))
      W=as.vector(t(W))
      W=matrix(W,nrow=length(u)*length(v),ncol=Nspp)
      return(W)
    }
    
    
    # Function to calculate size-dependent crowding, assuming no overlap
    wrijG=function(r,i,j){
      return(2*pi*integrate(function(z) z*exp(-alphaG[i,j]*(z^2))*Cr[[j]](z-r),r,r+r.U[j])$value+
               pi*Ctot[j]*exp(-alphaG[i,j]*((r+r.U[j])^2))/alphaG[i,j]);   
    }
    WrijG=Vectorize(wrijG,vectorize.args="r")
    
    wrijS=function(r,i,j){
      return(2*pi*integrate(function(z) z*exp(-alphaS[i,j]*(z^2))*Cr[[j]](z-r),r,r+r.U[j])$value+
               pi*Ctot[j]*exp(-alphaS[i,j]*((r+r.U[j])^2))/alphaS[i,j]);   
    }
    WrijS=Vectorize(wrijS,vectorize.args="r")
    
    
    # Function to sum total cover of each species
    sumCover=function(v,nt,h,A){
      out=lapply(1:Nspp,function(i,v,nt,h,A) h[i]*sum(nt[[i]]*exp(v[[i]]))/A,v=v,nt=nt,h=h,A=A)
      return(unlist(out))
    } 
    
    # Function to sum total density of each species
    sumN=function(nt,h){
      out=lapply(1:Nspp,function(i,nt,h) h[i]*sum(nt[[i]]),nt=nt,h=h)
      return(unlist(out))
    }
    
    # Function to calculate size variance of each species
    varN=function(v,nt,h,Xbar,N){
      out=lapply(1:Nspp,function(i,v,nt,h,Xbar,N) h[i]*sum((exp(v[[i]]-Xbar[i])^2)*nt[[i]])/N[i],v=v,nt=nt,h=h,Xbar=Xbar,N=N)
      return(unlist(out))
    }  
    
    # Function to do an image plot of a matrix in the usual orientation, A(1,1) at top left  
    matrix.image=function(x,y,A,col=topo.colors(100),...) {
      nx=length(x); ny=length(y); 
      x1=c(1.5*x[1]-0.5*x[2],1.5*x[nx]-0.5*x[nx-1]); 
      y1=c(1.5*y[1]-0.5*y[2],1.5*y[ny]-0.5*y[ny-1]); 
      image(list(x=x,y=y,z=t(A)),xlim=x1,ylim=rev(y1),col=col,bty="u",...);  
    }
    
    #============================================================================================#
    # (IV) Calculate the equilibrium areas.
    #============================================================================================# 
    fixI=which(sppList==fixSpp)
    
    ## initial population density vector
    nt=v
    for(i in 1:Nspp) nt[[i]][]=0
    # set fix spp to stable size distribution
    nt.fix=sizeD$freq
    # initialize at fix cover value
    tmp=fixCov*100/(h[fixI]*sum(nt.fix*exp(v[[fixI]])))
    nt.fix=nt.fix*tmp
    nt[[fixI]]=nt.fix
    new.nt=nt
    
    # set up matrix to record cover
    covSave=matrix(NA,0,(2+2*Nspp))
    colnames(covSave)=c("time","yrParams",paste(sppList,".t0",sep=""),paste(sppList,".t1",sep=""))
    covSave=rbind(covSave,c(1,NA,sumCover(v,nt,h,A),rep(NA,Nspp)) )
    
    # initial densities 
    Nsave=matrix(NA,tlimit,Nspp)
    Nsave[1,]=sumN(nt,h)
    
    yrSave=rep(NA,tlimit)
    for (i in 2:(tlimit)){
      
      #draw from observed year effects
      allYrs=c(1:Nyrs)
      doYear=sample(allYrs,1)
      yrSave[i]=doYear
      
      #get recruits per area
      #   cover=covSave[i-1,]; N=Nsave[i-1,]
      #   rpa=get.rpa(Rpars,cover,doYear)
      cover=covSave[i-1,3:6]; N=Nsave[i-1,]
      rpa=get.rpa(Rpars,cover,doYear)
      
      #calculate size-specific crowding
      alphaG=Gpars$alpha 
      alphaS=Spars$alpha 
      
      
      if(NoOverlap.Inter==F){#T: heterospecific genets cannot overlap; F: overlap allowed
        for(ii in 1:Nspp){ 
          # first do all overlap W's
          Xbar=cover*A/N       # multiply by A to get cover back in cm^2
          varX=varN(v,nt,h,Xbar,N) 
          
          muWG = pi*Xbar*N/(A*alphaG[ii,])
          muWS = pi*Xbar*N/(A*alphaS[ii,])
          
          muWG[is.na(muWG)]=0
          muWS[is.na(muWS)]=0
          
          WmatG[[ii]]=matrix(muWG,nrow=length(v[[ii]]),ncol=Nspp,byrow=T)
          WmatS[[ii]]=matrix(muWS,nrow=length(v[[ii]]),ncol=Nspp,byrow=T)
          
          # now do conspecific no overlap W
          Ctot[ii]=h[ii]*sum(expv[[ii]]*nt[[ii]]) 
          Cr[[ii]]=splinefun(b.r[[ii]],h[ii]*c(0,cumsum(expv[[ii]]*nt[[ii]])),method="natural")
          
          WmatG[[ii]][,ii]=WrijG(v.r[[ii]],ii,ii)/A
          WmatS[[ii]][,ii]=WrijS(v.r[[ii]],ii,ii)/A
        }
      }else{
        for(ii in 1:Nspp){
          Ctot[ii]=h[ii]*sum(expv[[ii]]*nt[[ii]]) 
          Cr[[ii]]=splinefun(b.r[[ii]],h[ii]*c(0,cumsum(expv[[ii]]*nt[[ii]])),method="natural") 
        }
        for(jj in 1:Nspp){ 
          
          WfunG=splinefun(size.range,WrijG(size.range,jj,jj))
          WfunS=splinefun(size.range,WrijS(size.range,jj,jj))
          
          for(ii in 1:Nspp) { 
            WmatG[[ii]][,jj]=WfunG(v.r[[ii]])/A 
            WmatS[[ii]][,jj]=WfunS(v.r[[ii]])/A 
          }
        }
        
      } # end NoOverlap if
      
      for(doSpp in 1:Nspp){  
        if(cover[doSpp]>0){    
          # make kernels and project
          K.matrix=make.K.matrix(v[[doSpp]],WmatG[[doSpp]],WmatS[[doSpp]],Rpars,rpa,Gpars,Spars,doYear,doSpp)  
          new.nt[[doSpp]]=K.matrix%*%nt[[doSpp]] 
          #       sizeSave[[doSpp]][,i]=new.nt[[doSpp]]/sum(new.nt[[doSpp]])  
        }    
      } # next species
      
      tmp=c(i,doYear,sumCover(v,nt,h,A),sumCover(v,new.nt,h,A))
      covSave=rbind(covSave,tmp)  # store the cover as cm^2/cm^2
      Nsave[i,]=sumN(nt,h)
      nt=new.nt
      
      # return focal spp to fix cover value
      tmp=fixCov*100/(h[fixI]*sum(nt[[fixI]]*exp(v[[fixI]])))
      nt[[fixI]]=nt[[fixI]]*tmp
      
      # return all other species to zero cover value
      tmp2 <- which(c(1:4) != fixI)
      nt[[tmp2[1]]] <- nt[[tmp2[1]]]*0
      nt[[tmp2[2]]] <- nt[[tmp2[2]]]*0
      nt[[tmp2[3]]] <- nt[[tmp2[3]]]*0
      
      print(i);flush.console()
      if(sum(is.na(nt))>0) browser()    
    } # next time step
    #low density growth rate of focal species
    tmp1 <- which(colnames(covSave)==paste(fixSpp, ".t0", sep=""))
    tmp2 <- which(colnames(covSave)==paste(fixSpp, ".t1", sep=""))
    pgrMean <- mean(log(covSave[burn.in:tlimit,tmp2]/covSave[burn.in:tlimit,tmp1]), na.rm=TRUE)
    maxR[jjjj,doSens] <- pgrMean
  } #next perturbation
  
}#next species

barplot(t(maxR), beside=T, names.arg=sppList, legend.text = changeVec, 
        xlab="Species", ylab="intrinsic growth rate")
plotD <- as.data.frame(t(maxR))
colnames(plotD) <- sppList
library(reshape2)
mD <- melt(plotD)
mD$change <- rep(changeVec, times = length(sppList))
library(ggplot2)
ggplot(mD, aes(x=variable, y=value, fill=as.factor(change)))+
  geom_bar(stat="identity", position="dodge")+
  theme_bw()+
  scale_fill_manual(values=c("grey25", "grey50", "grey75"),
                    name=("Proportional intercept change"))+
  xlab("species") + ylab("intrinsic growth rate (r)")

####
#### OUTPUT
####

outfile <- "Rmax_ipm.csv"
output <- as.data.frame(maxR)
output$species <- sppList
write.table(output, outfile, row.names=FALSE, sep=",")


