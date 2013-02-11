rm(list=ls())

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
      rand.mean <- round(runif(1, 0.1, 0.9), 1)
      beta[i,j] <- rnorm(1, rand.mean, 0.02)
    }
    
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



years = 600
N = 3
two.species <- model(Nspecies=N, time=years)


t = seq(100,years,1)
N1 = two.species[100:years,1] 
N2 = two.species[100:years,2] 
Ntot = two.species[100:years,3]

min1 = min(N1)
min2 = min(N2)
ymin = min(N1, N2)

plot(t, Ntot, type="l", ylim=c(ymin+500, max(Ntot)+500))
lines(t, N1, col="blue")
lines(t, N2, col="red")



