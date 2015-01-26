# Import and format survival parameters
# then define survival function

SparMCMC <- list()
tlimit=2500 ## number of years to simulate
sppList=c("ARTR","HECO","POSE","PSSP")
Nyrs=22

Nspp=length(sppList)
for(jjj in 1:tlimit){
  # survival parameters
  Spars=list(intcpt=rep(NA,Nspp),intcpt.yr=matrix(0,Nyrs,Nspp),intcpt.gr=matrix(0,6,Nspp),
             slope=rep(NA,Nspp),slope.yr=matrix(0,Nyrs,Nspp),
             nb=matrix(0,Nspp,Nspp),alpha=matrix(NA,Nspp,Nspp))
  
  #Get random row selection from the MCMC chain
  mcDraw <- sample(seq(1,2000,1), size = 1)
  
  for(i in 1:Nspp){
    infile=paste("survParamsMCMC_",sppList[i],".rds",sep="")
    Sdata=readRDS(infile)
    Sdata <- as.data.frame(Sdata[mcDraw,]) #get one row from the MCMC chain
    Sdata$Coef <- c(rep("W", 4),
                    rep("W.tau", 1),
                    rep("Group", 6),
                    rep("Intercept", 1),
                    rep("Intercept.tau", 1),
                    rep("Intercept.yr", 22),
                    rep("logarea.t0", 1),
                    rep("logarea.t0.tau", 1),
                    rep("logarea.t0.yr", 22),
                    rep("tauGroup", 1))
    colnames(Sdata)[1] <- "value"
    
    Spars$intcpt[i]=Sdata$value[which(Sdata$Coef=="Intercept")]
    
    tmp=which(Sdata$Coef=="Group")
    if(length(tmp)>0) Spars$intcpt.gr[,i]=Sdata$value[tmp] 
    
    tmp=which(Sdata$Coef=="Intercept.yr")
    if(length(tmp)>0) Spars$intcpt.yr[,i]=Sdata$value[tmp] 
    
    tmp=which(Sdata$Coef=="logarea.t0")
    Spars$slope[i]=Sdata$value[tmp]
    
    # random effects on slope
    tmp=which(Sdata$Coef=="logarea.t0.yr")
    if(length(tmp)>0) Spars$slope.yr[,i]=Sdata$value[tmp]
    
    # get competition coefficients
    tmp=which(Sdata$Coef=="W")
    if(length(tmp)>0) Spars$nb[i,]=Sdata$value[tmp]
    
    alphaS <- read.csv("alphaSurvival.csv")
    Spars$alpha[i,]=alphaS$alpha
  } # next i
  rm(Sdata)
  SparMCMC[[length(SparMCMC)+1]] <- Spars
  print(paste("Retrieved survival parameter set ", jjj, "out of", tlimit, "."))
}
saveRDS(SparMCMC,"survivalParamsMCMC4IPM.rds")

# growth function
## survival function: probability an individual of size u survives  (u is on log scale)
# S=function(u,W,Spars,doYear,doSpp){
#   mu=Spars$intcpt[doSpp]+Spars$intcpt.yr[doYear,doSpp]+(Spars$slope[doSpp]+Spars$slope.yr[doYear,doSpp])*u+
#     W%*%(Spars$nb[doSpp,])
#   return(inv.logit(mu))
# }
