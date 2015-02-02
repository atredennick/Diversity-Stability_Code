##################################################
#### Call BUGS models for vital rate regressions
####
#### Andrew Tredennick: atredenn@gmail.com
#### 1-30-2015

####
#### Clear the workspace -----------------------
####
rm(list=ls(all=TRUE))


####
#### Load libraries -----------------------------
####
library(rjags)
library(reshape2)
library(plyr)


####
#### Source some functions ----------------------
####
source("../../Functions/get_crowd_fxn.R")


####
#### Set global variables -----------------------
####
data_dir <- "../../Data/"
site_list <- list.files(data_dir)
site_list <- site_list[which(site_list!="alpha_list.csv")]
alpha_all <- read.csv("../../Data/alpha_list.csv")
colnames(alpha_all) <- tolower(colnames(alpha_all))
n_sites <- length(site_list)


####
#### Begin looping through sites and species for vital rate regressions ----------------------
####
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
    
    ## Set up data for JAGS MCMC
    Nyears = length(unique(grow_data$year))
    Ngroups = length(unique(grow_data$group))
    Nspp = length(spp_list)
    data_jags <- list(yr = as.numeric(as.factor(grow_data$year)),
                      Nyears = Nyears,
                      grp = as.numeric(grow_data$group),
                      Ngroups = Ngroups,
                      Nspp = Nspp,
                      crowd = crowd,
                      x = log(grow_data$area.t0),
                      y = log(grow_data$area.t1),
                      Nobs = nrow(grow_data))
    params=c("intcpt.mu","intcpt.tau","intcpt.yr","intG",
             "slope.mu","slope.tau","slope.yr",
             "NBbeta.mu","NBbeta.tau","tauGroup",
             "tau","tauSize")
    
    inits=list(1)
    inits[[1]]=list(intcpt.mu=0, intcpt.tau=1, intcpt.yr=rep(0,Nyears),
                    intG=rep(0.1,Ngroups), slope.mu=0, slope.tau=1,
                    slope.yr=rep(0,Nyears), NBbeta.mu=rep(-0.2,Nspp),
                    NBbeta.tau=0.2, tauGroup=2, tau=1, tauSize=1)
    
    inits[[2]]=list(intcpt.mu=-1, intcpt.tau=0.1, intcpt.yr=rep(-1,Nyears),
                    intG=rep(-0.5,Ngroups), slope.mu=-1, slope.tau=0.1,
                    slope.yr=rep(-1,Nyears), NBbeta.mu=rep(-0.5,Nspp),
                    NBbeta.tau=0.1, tauGroup=1, tau=0.1, tauSize=0.1)
    
    iterations <- 50000
    burn_in <- 10000
    adapt <- 1000
    print(paste("Starting on", spp_list[spp_now], "for", site, sep=""))
    mod <- jags.model("../../Functions/growth_BUGS.R", data=data_jags, inits=inits, n.chains=2, n.adapt=adapt)
    update(mod, n.iter = (burn_in))
    out <- coda.samples(mod, params, n.iter=iterations, n.thin=50)
    
    out_stats <- summary(out)$stat
    out_quants <- summary(out)$quantile
    out_gelman <- gelman.diag(out)
    outMCMC <- rbind(out[[1]], out[[2]])
    
    out_dir <- paste("../vitalRateResults/", site, "/", sep="")
    out_file_stats <- paste(out_dir, "growth_stats_", site, "_", spp_list[spp_now], ".csv", sep="")
    out_file_quants <- paste(out_dir, "growth_quants_", site, "_", spp_list[spp_now], ".csv", sep="")
    out_file_gelman <- paste(out_dir, "growth_gelman_", site, "_", spp_list[spp_now], ".csv", sep="")
    out_file_grow <- paste(out_dir,"MCMC_growth_", site, "_", spp_list[spp_now], ".rds", sep="")
    
    write.csv(out_stats, out_file_stats)
    write.csv(out_quants, out_file_quants)
    write.csv(out_gelman, out_file_gelman)
    saveRDS(object = outMCMC, file = out_file_grow)
    print(paste("Done with", spp_list[spp_now], "for", site, sep=""))
  } # end species loop for vital rate regressions
} # end site loop
