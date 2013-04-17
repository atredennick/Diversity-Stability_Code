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
    
    evar[i] <- runif(1, 0.1, 0.2)
    dvar[i] <- runif(1, 0, 2)
    
    rm[i] <- runif(1, 0.2, 1.5)
    
    N[1,i] <- rnorm(1, K[i], 500)
    
    for (j in 1:(Nspecies-1)){
      rand.mean <- round(runif(1, 0, 0.05), 1)
      beta[i,j] <- rnorm(1, rand.mean, 0.02)
    }
  }

  
  Ntot[1] <- sum(N[1,])
  
  for(t in 2:time){
    for(i in 1:Nspecies){
      
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



years = 2000
N.sim = c(2,4,8,16,32)
# N.sim = c(2)
sims = 500
cv = matrix(nrow=sims, ncol=length(N.sim))


  for(i in 1:length(N.sim)){
    for(j in 1:sims){
      two.species <- model(Nspecies=N.sim[i], time=years)
      Ntot = two.species[1001:years,ncol(two.species)]
      cv[j,i] = mean(Ntot)/sd(Ntot)
    }
  }


write.csv(cv, "stability_simulations_nodomhivar.csv")
