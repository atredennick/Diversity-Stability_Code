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
                  run_time, whitevar){
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
    N2[t] = N2[t-1] + N2[t-1]*r2[t]
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
# Assume symmetrical community
beta12 = beta21 = seq(0,0.9,0.05)
rho = seq(-0.9, 0.9, 0.1)
run_time = 1500
burn = run_time/2
n_sims=20

# Set up storage matrices and lists
out_comm <- matrix(NA, ncol=length(beta12), nrow=length(rho))
out_resp <- matrix(NA, ncol=length(beta12), nrow=length(rho))
tmp_comm <- numeric(n_sims)
tmp_resp <- numeric(n_sims)
slopes <- matrix(NA, nrow=n_sims, ncol=2)
rsquares <- matrix(NA, nrow=n_sims, ncol=2)
out_slopes <- array(dim = c(length(rho),length(beta12),2))
out_squares <- array(dim = c(length(rho),length(beta12),2))


# Run model at different levels of competition and correlation in env response
pb <- txtProgressBar(min=2, max=length(beta12), char="+", style=3, width=65)
for(comp in 1:length(beta12)){
  for(corr in 1:length(rho)){
    for(sim in 1:n_sims){
      # get environment time series
      whitevar <- get_env(sigE = 1, rho = rho[corr], nTime = run_time)
      
      model_null <- model(rm1, rm2, 
                          N1.start, N2.start, 
                          K1, K2,
                          evar1, evar2, 
                          beta12[comp], beta21[comp], 
                          run_time, whitevar)
      
      model_null <- as.data.frame(model_null)
      tmp_comm[sim] <- as.numeric(community.sync(model_null[burn:run_time,1:2])[1])
      tmp_resp[sim] <- as.numeric(community.sync(whitevar[burn:run_time,1:2])[1])
      mod1 <- lm(model_null[burn:run_time,4]~whitevar[burn:run_time,1])
      mod2 <- lm(model_null[burn:run_time,5]~whitevar[burn:run_time,2])
      slopes[sim,] <- c(coef(mod1)[2], coef(mod2)[2])
      rsquares[sim,] <- c(summary(mod1)$r.squared, summary(mod2)$r.squared)
    }#end simulation loop
    out_comm[corr, comp] <- mean(tmp_comm)
    out_resp[corr, comp] <- mean(tmp_resp)
    out_slopes[corr, comp, 1] <- mean(slopes[,1])
    out_slopes[corr, comp, 2] <- mean(slopes[,2])
    out_squares[corr, comp, 1] <- mean(rsquares[,1])
    out_squares[corr, comp, 2] <- mean(rsquares[,2])
  }#end correlation loop
  setTxtProgressBar(pb, comp)
}#end competition loop

colnames(out_comm) <- beta12
rownames(out_comm) <- rho
comm_df <- melt(out_comm)
colnames(comm_df) <- c("rho", "beta", "community")
colnames(out_resp) <- beta12
rownames(out_resp) <- rho
resp_df <- melt(out_resp)
comm_df$response <- resp_df$value
comm_df$difference <- with(comm_df, community-response)

slopes_df <- out_slopes[,,1]
colnames(slopes_df) <- beta12
rownames(slopes_df) <- rho
slopes_m <- melt(slopes_df)
colnames(slopes_m) <-  c("rho", "beta", "slope")
slopes_m$species <- rep(1,nrow(slopes_m))
tmp <- out_slopes[,,2]
colnames(tmp) <- beta12
rownames(tmp) <- rho
tmp_m <- melt(tmp)
colnames(tmp_m) <-  c("rho", "beta", "slope")
tmp_m$species <- rep(2,nrow(tmp_m))
slopes_m <- rbind(slopes_m, tmp_m)

squares_df <- out_squares[,,1]
colnames(squares_df) <- beta12
rownames(squares_df) <- rho
squares_m <- melt(squares_df)
colnames(squares_m) <-  c("rho", "beta", "rsquare")
squares_m$species <- rep(1,nrow(squares_m))
tmp <- out_squares[,,2]
colnames(tmp) <- beta12
rownames(tmp) <- rho
tmp_m <- melt(tmp)
colnames(tmp_m) <-  c("rho", "beta", "rsquare")
tmp_m$species <- rep(2,nrow(tmp_m))
squares_m <- rbind(squares_m, tmp_m)

####
#### 4. Make some plots -------------------------------
####
ggplot(comm_df, aes(x=beta, y=difference, color=rho, group=rho))+
  geom_line(size=1.5, alpha=0.75)+
#   geom_point(size=5, color="white")+
#   geom_point(size=3)+
  xlab(expression(paste("Competition coefficient (", beta[ij], ")")))+
  ylab("Community synchrony - Environmental synchrony")+
  scale_color_gradient(limits=c(-1, 1), low="red")+
  theme_bw()

# ggplot(slopes_m, aes(x=beta, y=slope, color=rho, group=rho))+
#   geom_line()+
#   facet_grid(species~.)
# 
# ggplot(squares_m, aes(x=beta, y=rsquare, color=rho, group=rho))+
#   geom_line()+
#   facet_grid(species~.)
