rm(list=ls())

#set working directory
# wd = "/Users/atredenn/Documents/Projects/Diversity_Stability/"
# setwd(wd)
# getwd()

#for source file
# "/Users/atredenn/Documents/Projects/Diversity_Stability/Diversity-Stability_Code/LoreauModel_Multi-species.R"

model <- function(Nspecies, time){
  rm = numeric(Nspecies)
  N = matrix(ncol = Nspecies, nrow = time)
  K = numeric(Nspecies)
  K.log = numeric(Nspecies)
  evar = numeric(Nspecies)
  dvar = numeric(Nspecies)
  beta = matrix(nrow = Nspecies, ncol = (Nspecies-1))
  Ntot = numeric(time)
  r = matrix(ncol=Nspecies, nrow=time)
  comp.spp = matrix(nrow=Nspecies, ncol=(Nspecies-1))
  sum.comp = numeric(Nspecies)

  for (i in 1:Nspecies){
    K.log[i] <- rnorm(1, log(10000), 0.7)
    K[i] <- exp(K.log[i])
    
    evar[i] <- runif(1, 0, 0.02)
    dvar[i] <- runif(1, 0, 2)
    
    rm[i] <- runif(1, 0.2, 1.5)
    
    N[1,i] <- rnorm(1, K[i], 500)
    
    for (j in 1:(Nspecies-1)){
      rand.mean <- round(runif(1, 0, 0.05), 1)
      beta[i,j] <- rnorm(1, rand.mean, 0.02)
    }
  
    beta[2:Nspecies,1] <- 0.4
    evar[1] <- runif(1, 0.1, 0.2)
    
  }
  
  Ntot[1] <- sum(N[1,])
  
  for(t in 2:time){
    for(i in 1:Nspecies){
      
#       if (i != Nspecies){
#         spp <- seq(i+1, Nspecies, 1)
#         leftover <- (Nspecies-1) - (length(spp))
#         if (leftover == 0) {
#           spp <- spp
#         }
#         else{
#           spp2 <- seq(1, leftover, 1)
#           spp <- c(spp, spp2)
#         }
#          
#       }
#       else{
#         spp <- seq(1, (Nspecies-1), 1)
#       }
#       
      n.seq = seq(1, Nspecies, 1)
      n.index = n.seq[n.seq!=i]
      
      for(j in 1:(Nspecies-1)){
        comp.spp[i,j] <- ((beta[i,j] * N[t-1,n.index[j]])/K[n.index[j]])
      }
      sum.comp[i] <- sum(comp.spp[i,])
      whitevar <- rnorm(1, mean=0)
      dnoise <- rnorm(1, mean=0)
      
      r[t, i] <- rm[i] * (1-(N[t-1,i]/K[i]) - sum.comp[i]) + evar[i] * whitevar + (dvar[i] * dnoise/sqrt(N[t-1, i]))    
      N[t, i] <- N[t-1, i] + N[t-1, i] * r[t, i]
    }
    
    Ntot[t] = sum(N[t,])
  }
  
  population <- matrix(ncol=Nspecies, nrow=time)
  population <- cbind(N, Ntot)
  return(population)
}



years = 10
N.sim = c(2,4,8,16,32)
sims = 10
cv = matrix(ncol=length(N.sim), nrow=sims)

for(i in 1:length(N.sim)){
  for(j in 1:sims){
    two.species <- model(Nspecies=N.sim[i], time=years)
    Ntot = two.species[1:years,ncol(two.species)]
    cv[j,i] = mean(Ntot)/sd(Ntot)
  }
  
}


cv = cv[complete.cases(cv),]












write.csv(cv, "stability_simulations.csv")


N = matrix(ncol=5, nrow=length(cv[,1]))
N[,1]=2
N[,2]=4
N[,3]=8
N[,4]=16
N[,5]=32
# N=log(N)
plot(N[,1], cv[,1], xlim=c(min(N)-0.4,max(N)+0.4), pch=21, col="#00000050", cex=0.3,
     xlab="Number of Species", ylab=expression(paste("Stability (", mu/sigma, ")")),
     ylim = c(min(cv), max(cv)), axes=FALSE)
axis(1, at=N[1,], labels=c("2", "4", "8", "16", "32"), cex.axis=0.8)
axis(2)
box()
points(N[,2], cv[,2], pch=21, col="#00000050", cex=0.3)
points(N[,3], cv[,3], pch=21, col="#00000050", cex=0.3)
points(N[,4], cv[,4], pch=21, col="#00000050", cex=0.3)
points(N[,5], cv[,5], pch=21, col="#00000050", cex=0.3)


means = c(mean(cv[,1]), mean(cv[,2]), mean(cv[,3]), mean(cv[,4]), mean(cv[,5]))
lines(N[1,], means, col="red")

points(N[1,1], mean(cv[,1]), pch=17, col="white", cex=2)
points(N[1,1], mean(cv[,1]), pch=17, col="red", cex=1)
points(N[1,2], mean(cv[,2]), pch=17, col="white", cex=2)
points(N[1,2], mean(cv[,2]), pch=17, col="red", cex=1)
points(N[1,3], mean(cv[,3]), pch=17, col="white", cex=2)
points(N[1,3], mean(cv[,3]), pch=17, col="red", cex=1)
points(N[1,4], mean(cv[,4]), pch=17, col="white", cex=2)
points(N[1,4], mean(cv[,4]), pch=17, col="red", cex=1)
points(N[1,5], mean(cv[,5]), pch=17, col="white", cex=2)
points(N[1,5], mean(cv[,5]), pch=17, col="red", cex=1)



cv.long = as.vector(cv)
spp = c(rep(N.sim[1], length(cv[,1])),
      rep(N.sim[2], length(cv[,1])),
      rep(N.sim[3], length(cv[,1])),
      rep(N.sim[4], length(cv[,1])),
      rep(N.sim[5], length(cv[,1])))

mod <- lm(cv.long ~ spp)
summary(mod)
pred = predict(mod)
pred = unique(pred)
lines(N[1,], pred)
