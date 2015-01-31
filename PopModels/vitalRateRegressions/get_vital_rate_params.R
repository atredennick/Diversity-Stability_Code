##################################################
#### Call BUGS models for vital rate regressions
####
#### Andrew Tredennick: atredenn@gmail.com
#### 1-30-2015

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# IMPORTANT: NIMBLE DOES NOT LIKE GLOBAL ENVIRONMENT     #
# TO BE IN THE PANEL OF RSTUDIO. CHANGE TO PACKAGE:STATS #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

####
#### Clear the workspace -----------------------
####
rm(list=ls(all=TRUE))


####
#### Load libraries -----------------------------
####
library(nimble)
library(reshape2)
library(plyr)

####
#### Source some functions
####
source("../../Functions/get_crowd_fxn.R")
source("../../Functions/growth_BUGS.R")


data_dir <- "../../Data/"
site_list <- list.files(data_dir)
site_list <- site_list[which(site_list!="alpha_list.csv")]
alpha_all <- read.csv("../../Data/alpha_list.csv")
colnames(alpha_all) <- tolower(colnames(alpha_all))
n_sites <- length(site_list)

for(site_now in 1:n_sites){
  site <- site_list[site_now]
  site_dir <- paste(data_dir, site, "/", sep="")
  spp_list <- list.files(site_dir)
  n_spp <- length(spp_list)
  alpha_effect <- subset(alpha_all, site==site)$alpha
  
  dist_data <- data.frame("quad"=NA,"year"=NA,"trackID"=NA,
                          "area"=NA,"x"=NA,"y"=NA, "nbSpp"=NA)
  for(spp_now in 1:n_spp){
    dist_file <- paste(site_dir, spp_list[spp_now], "/", spp_list[spp_now], "_genet_xy.csv", sep="")
    tmp_dist <- read.csv(dist_file)
    tmp_dist$nbSpp <- spp_list[spp_now]
    tmp_dist <- tmp_dist[,c("quad","year","trackID","area","x","y","nbSpp")]
    dist_data <- rbind(dist_data, tmp_dist)
  }
  dist_data <- dist_data[2:nrow(dist_data),]

  for(spp_now in 1:n_spp){
    ## Get growth data
    grow_file <- paste(site_dir, spp_list[spp_now], "/growDnoNA.csv", sep="")
    grow_data <- read.csv(grow_file)
    
    ## Add group column ONLY if missing
    if("group" %in% colnames(grow_data)==FALSE) {grow_data$group <- as.factor(substr(grow_data$quad,1,1))}
    
    ## Calculate observed crowding
    crowd <- get_crowding(grow_data, dist_data, n_spp, spp_list[spp_now], 
                          spp_list, alpha_effect)
    # Need to set the very small crowding value equal 
    # to zero (<1e-100, but bigger than 0), to pass BUGS (C. Chu)
    crowd[crowd<1e-99] <- 0
    
    ## Set up data for NIMBLE MCMC
    Nspp <- n_spp
    Ngroups <- length(unique(grow_data$group))
    Nyears <- length(unique(grow_data$year))
    G <- as.numeric(grow_data$group)
    Yr <- as.numeric(as.factor(grow_data$year))
    logsize0 <- round(log(grow_data$area.t0),4)
    logsize1 <- log(grow_data$area.t1)
    totN <- length(logsize1)
    
    # Get initial values from quick lm
    mod <- lm(logsize1~logsize0+crowd)
    
    constants <- list(totN=totN, lagsize=logsize0)
    data <- list(size=logsize1)
    g_code <- nimbleCode({
      for(i in 1:totN){
        size[i] ~ dnorm(mu[i], tau)
        mu[i] <- intcpt + slope*lagsize[i]
      }
      tau ~ dgamma(0.001, 0.001)
      intcpt ~ dnorm(0, 0.001)
      slope ~ dnorm(0, 0.001)
    })
    Rmodel <- nimbleModel(code = g_code, 
                          constants = constants, 
                          data = data)
    mcmcspec <- configureMCMC(Rmodel, print=FALSE, thin=10)
    Rmcmc <- buildMCMC(mcmcspec)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc, project = Cmodel)
    Cmcmc$run(10000)
    Cmcmc$run(50000)
    samples <- as.matrix(Cmcmc$mvSamples)
    
  } # end species loop for vital rate regressions
} # end site loop
