################################################################################
#### Lotka-Volterra model to test the relationship between
#### synchrony in population dynamics and species responses to 
#### the environment in the presence/absence of competition.
####
#### The idea is to see if the correlation between the synchrony measures
#### can tell us something about how much competition structures
#### the community if we know species response to environment.
####
#### Andrew Tredennick: atredenn@gmail.com
####
#### Model adapted from Loreau and de Mazancourt 2013, Ecology Letters
################################################################################

# Clear the workspace
rm(list=ls())

# Load some libraries
library(mvtnorm)
library(ggplot2)
library(synchrony)
library(msm)
library(reshape2)


####
#### 1. Define the environmental response function -----------------------------
####
get_env <- function(sigE, rho, nTime){
  varcov <- matrix(c(sigE, rho*sigE, rho*sigE, sigE), 2, 2)
  e <- rmvnorm(n = nTime, mean = c(0,0), sigma = varcov)
  return(e)
}

####
#### 2. Define the 2 species competition model function ------------------------
####
model <- function(rm1, rm2, 
                  N1.start, N2.start, 
                  K1, K2,
                  evar1, evar2, 
                  beta12, beta21, 
                  run_time, whitevar, harv){
  r1 = numeric(run_time)
  r2 = numeric(run_time)
  N1 = numeric(run_time)
  N2 = numeric(run_time)
  Ntot = numeric(run_time)
  
  N1[1] <- N1.start
  N2[1] <- N2.start
  Ntot[1] <- N1[1] + N2[1]
  
  for(t in 2:run_time){
    sum.non1 = (beta12 * N2[t-1])/K2
    sum.non2 = (beta21 * N1[t-1])/K1
    r1[t] = rm1 * (1-(N1[t-1]/K1) - sum.non1) + evar1 * whitevar[t,1]
    r2[t] = rm2 * (1-(N2[t-1]/K2) - sum.non2) + evar2 * whitevar[t,2] 
    N1[t] = N1[t-1] + N1[t-1]*r1[t] 
    N2[t] = N2[t-1] + N2[t-1]*r2[t] - N2[t-1]*harv[t]
    Ntot[t] = N1[t] + N2[t]
  }
  
  pops <- matrix(ncol=5, nrow=run_time)
  pops[,1] <- N1
  pops[,2] <- N2
  pops[,3] <- Ntot
  pops[,4] <- r1
  pops[,5] <- r2
  return(pops)
}


####
#### 3. Run model simulations across range of competition ------------------- 
####    and species correlations in env response
####

# Set up global variables
rm1 = 0.8 #species 1 intrinsic growth rate
rm2 = 0.8 #species 2 intrinsic growth rate
K1 = 1500 #species 1 carrying capacity
K2 = 1500 #species 2 carrying capacity
evar1 = 0.05 #species 1 environmental variance
evar2 = 0.05 #species 2 environmental variance
N1.start <- K1
N2.start <- K2
# beta12 = 1 #competition coefficient; effect of spp2 on spp1
# beta21 = 1 #competition coefficient; effect of spp1 on spp2
beta12 =  1
beta21 = 1
rho = 0
run_time = 1500
burn = run_time/2

whitevar <- get_env(sigE = 1, rho = rho, nTime = run_time)

model_null <- model(rm1, rm2, 
                    N1.start, N2.start, 
                    K1, K2,
                    evar1, evar2, 
                    beta12, beta21, 
                    run_time, whitevar,
                    harv=rep(0,run_time))

par(mfrow=c(1,1))
matplot(model_null[(run_time-500):run_time,1:2],type="l", lty=1, col=c("goldenrod2","darkslateblue"),
        xlab="time", ylab="biomass")
matplot(model_null[(run_time-500):run_time,1:2],pch=19, add=TRUE, col=c("goldenrod2","darkslateblue"))
# matplot(model_null[(run_time-50):run_time,4:5],type="l", lty=1, col=c("goldenrod2","darkslateblue"),
#         xlab="time", ylab="per capita growth rate")
# matplot(model_null[(run_time-50):run_time,4:5],pch=19, add=TRUE, col=c("goldenrod2","darkslateblue"))


# print(paste("ABUNDANCE:",as.numeric(community.sync(model_null[burn:run_time,1:2])[1])))
print(paste("COMM PGR:",as.numeric(community.sync(model_null[burn:run_time,4:5])[1])))
print(paste("ENV RESP:",as.numeric(community.sync(whitevar[burn:run_time,1:2])[1])))



