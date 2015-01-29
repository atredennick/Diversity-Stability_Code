model{
  for(i in 1:Nobs){
    logit(mu[i]) <- intcpt.yr[yr[i]] + intG[grp[i]] + slope.yr[yr[i]]*x[i] + inprod(NBbeta.mu[1:Nspp],crowd[i,])
    mu2[i]<-max(0.000000000000,min(0.999999999999,mu[i]))
    y[i] ~ dbern(mu2[i])
  }
  
  #priors
  intcpt.mu~dnorm(0,0.001)
  intcpt.tau~dgamma(0.5,0.5)
  slope.mu~dnorm(0,0.001)
  slope.tau~dgamma(0.5,0.5)
  tauGroup~dgamma(0.001,0.001)
  NBbeta.tau~dgamma(0.5,0.5)
  for(j in 1:Ngroups){
    intG[j]~dnorm(0,tauGroup)
  }  
  for(k in 1:Nyears){
    slope.yr[k]~dnorm(slope.mu,slope.tau)
    intcpt.yr[k]~dnorm(intcpt.mu,intcpt.tau)
  }  
  for(m in 1:Nspp){
    NBbeta.mu[m]~dnorm(0,NBbeta.tau) ##OK for both
  }
}