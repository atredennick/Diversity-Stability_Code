##########################################################################
#### Script that calls 'GetTimeSeriesMetricsFxn.R' for each site and set
#### of quadrats. Outputs tables of metrics and plots of metrics and
#### time series.
####
#### Andrew Tredennick: atredenn@gmail.com
#### Date: 1-28-2015
####

#Clear the workspace
rm(list=ls(all=TRUE))

####
#### Load necessary libraries ----------------------------------------
####
library(ggplot2); library(plyr); library(reshape2); library(gridExtra)

####
#### Source the function ----------------------------------------------
####
source("GetTimeSeriesMetricsFxn.R")

####
#### Set up data structures to send to the function -------------------
####
sites <- list.files("../Data")
species <- list()
for(s in 1:length(sites)){
  tmp_spp <- list.files(paste("../Data/", sites[s], sep=""))
  tmp_name <- sites[s]
  species <- c(species, list(site_name = tmp_spp))
}
names(species) <- sites

####
#### Loop through sites and pass values to the function ---------------
####
n_sites <- length(species)
quad_ts_metrics <- list()
for(ss in 1:n_sites){
  tmp_metrics <- get_comm_ts_metrics(spp_list = species[[ss]], site = names(species)[ss])
  quad_ts_metrics[[ss]] <- tmp_metrics
}
names(quad_ts_metrics) <- sites

#Save the returned list
saveRDS(quad_ts_metrics, "quad_ts_metrics_AllSites.rds")
