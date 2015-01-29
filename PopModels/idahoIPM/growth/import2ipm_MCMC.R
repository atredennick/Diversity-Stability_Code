# Import and format growth parameters
# then define growth function

GparMCMC <- list()

for(jjj in 1:tlimit){
  # growth parameters
  Gpars=list(intcpt=rep(NA,Nspp),intcpt.yr=matrix(0,Nyrs,Nspp),intcpt.gr=matrix(0,6,Nspp),
             slope=rep(NA,Nspp),slope.yr=matrix(0,Nyrs,Nspp),
             nb=matrix(0,Nspp,Nspp),alpha=matrix(NA,Nspp,Nspp),
             sigma2.a=rep(NA,Nspp),sigma2.b=rep(NA,Nspp))
  
  #Get random row selection from the MCMC chain
  mcDraw <- sample(seq(1,2000,1), size = 1)
  
  for(i in 1:Nspp){
    infile=paste("growth/Growth_paramsMCMC_",sppList[i],".rds",sep="")
    Gdata=readRDS(infile)
    Gdata <- as.data.frame(Gdata[mcDraw,]) #get one row from the MCMC chain
    Gdata$Coef <- c(rep("W", 4),
                    rep("W.tau", 1),
                    rep("Group", 6),
                    rep("Intercept", 1),
                    rep("Intercept.tau", 1),
                    rep("Intercept.yr", 22),
                    rep("logarea.t0", 1),
                    rep("logarea.t0.tau", 1),
                    rep("logarea.t0.yr", 22),
                    rep("sigma2.a", 1),
                    rep("tauGroup", 1),
                    rep("sigma2.b", 1))
    colnames(Gdata)[1] <- "value"
    
    Gpars$intcpt[i]=Gdata$value[which(Gdata$Coef=="Intercept")]
    
    tmp=which(Gdata$Coef=="Group")
    if(length(tmp)>0) Gpars$intcpt.gr[,i]=Gdata$value[tmp] 
    
    tmp=which(Gdata$Coef=="Intercept.yr")
    if(length(tmp)>0) Gpars$intcpt.yr[,i]=Gdata$value[tmp] 
    
    tmp=which(Gdata$Coef=="logarea.t0")
    Gpars$slope[i]=Gdata$value[tmp]
    
    # random effects on slope
    tmp=which(Gdata$Coef=="logarea.t0.yr")
    if(length(tmp)>0) Gpars$slope.yr[,i]=Gdata$value[tmp]
    
    # get competition coefficients
    tmp=which(Gdata$Coef=="W")
    if(length(tmp)>0) Gpars$nb[i,]=Gdata$value[tmp]
    
    alphaG <- read.csv("./growth/alphaGrowth.csv")
    Gpars$alpha[i,]=alphaG$alpha
    Gpars$sigma2.a[i]=Gdata$value[which(Gdata$Coef=="sigma2.a")]
    Gpars$sigma2.b[i]=Gdata$value[which(Gdata$Coef=="sigma2.b")]
  } # next i
  rm(Gdata)
  GparMCMC[[length(GparMCMC)+1]] <- Gpars
  print(paste("Done with parameter set ", jjj, "out of", tlimit, "."))
} #end list loop



# growth function
G=function(v,u,W,Gpars,doYear,doSpp){
  mu=Gpars$intcpt[doSpp]+Gpars$intcpt.yr[doYear,doSpp]+
    (Gpars$slope[doSpp]+Gpars$slope.yr[doYear,doSpp])*u+
    W%*%(Gpars$nb[doSpp,])
  
  sigma2=Gpars$sigma2.a[doSpp]*exp(Gpars$sigma2.b[doSpp]*mu)
  out=dnorm(v,mu,sqrt(sigma2))
  out
}
