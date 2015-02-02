##############################################
#### Vital rate regression for growth.
####
#### Andrew Tredennick: atredenn@gmail.com
#### 1-30-2015

####
#### Write BUGS code for NIMBLE
####
model{
  for(i in 1:Nobs){
    mu[i] <- intcpt.yr[yr[i]] + intG[grp[i]] + slope.yr[yr[i]]*x[i] + inprod(NBbeta.mu[1:Nspp],crowd[i,])
    tau2[i] <- 1/(tau*exp(tauSize*mu[i])) 
    tau3[i] <- max(tau2[i],0.00000001)  
    y[i] ~ dnorm(mu[i], tau3[i])
  }
  
  #priors
  tau~dnorm(0,0.001)
  tauSize~dnorm(0,0.001)
  
  for(j in 1:Ngroups){
    intG[j]~dnorm(0,tauGroup) ##using hyperprior here
  }
  tauGroup~dgamma(0.001,0.001)
  
  intcpt.mu~dnorm(0,0.001)
  slope.mu~dnorm(0,0.001)
  intcpt.tau~dgamma(0.5,0.5)
  slope.tau~dgamma(0.5,0.5)
  for(k in 1:Nyears){
    slope.yr[k]~dnorm(slope.mu,slope.tau)
    intcpt.yr[k]~dnorm(intcpt.mu,intcpt.tau)
  }
  
  for(m in 1:Nspp){
    NBbeta.mu[m]~dnorm(0,NBbeta.tau) T(-2,-0.000001) ##change this to meet BUGS(too wide to fail)
  }
  NBbeta.tau~dgamma(0.5,0.5)
}

