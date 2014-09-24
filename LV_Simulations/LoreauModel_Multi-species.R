rm(list=ls())

model <- function(Nspecies, time, evar.dom, evar.dom.sd){
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
  }
  beta[2:Nspecies,1] <- 0.4
  evar[1] <- runif(1, evar.dom, evar.dom.sd)
  K[1] <- max(K)
  
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
evar.vec = c(0.02, 0.1)
evar.vec.sd = c(0.02, 0.2)
sims = 500
cv = array(dim=c(length(evar.vec), sims, length(N.sim)))

for(k in 1:length(evar.vec)){
  for(i in 1:length(N.sim)){
    for(j in 1:sims){
      two.species <- model(Nspecies=N.sim[i], time=years, 
                           evar.dom=evar.vec[k], evar.dom.sd=evar.vec.sd[k])
      Ntot = two.species[1001:years,ncol(two.species)]
      cv[k,j,i] = mean(Ntot)/sd(Ntot)
    }
  }
}

write.csv(cv[1,,], "stability_simulations_lovar.csv")
write.csv(cv[2,,], "stability_simulations_hivar.csv")
