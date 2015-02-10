##########################################################################
#### Script to assess importance of demographic vs. environmental 
#### stochasticity in driving population growth rate variance
####
#### Andrew Tredennick: atredenn@gmail.com
#### Date: 2-10-2015
####

#Clear the workspace
rm(list=ls(all=TRUE))

####
#### Load necessary libraries ----------------------------------------
####
library(ggplot2); library(plyr); library(reshape2); library(gridExtra)


####
#### Set up data structures -------------------
####
sites <- list.files("../Data")
exclude <- which(sites==grep("^alpha*", sites, value = TRUE))
sites <- sites[-exclude]
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
all_data <- data.frame(quad=NA, year=NA, trackID=NA, area=NA, species=NA, site=NA)
for(ss in 1:n_sites){
  spp_list <- species[[ss]]
  num_spp <- length(spp_list)
  for(dospp in 1:num_spp){ #loop through species to read in data
    spp_now <- spp_list[dospp]
    genet_file <- paste("../Data/", sites[ss],"/",spp_now,"/",spp_now,"_genet_xy.csv",sep="")
    spp_data <- read.csv(genet_file)
    spp_data <- spp_data[,c("quad","year","trackID","area")]
    spp_data$species <- spp_now
    spp_data$site <- sites[ss]
    all_data <- rbind(all_data, spp_data)
  } #end species looping for raw data
}

all_data <- all_data[2:nrow(all_data),] #get rid of NA row
all_density <- ddply(all_data, .(quad, year, species), summarise,
                     num_genets = length(area))
lag_density <- all_density
lag_density$year <- lag_density$year+1
lag_density$lag_genets <- lag_density$num_genets
lag_density <- lag_density[,c("quad","year","species","lag_genets")]
full_density <- merge(all_density, lag_density)
pop_growth_rate <- ddply(full_density, .(species), summarise,
                         pop_growth_rate = mean(log(num_genets/lag_genets)),
                         pop_grate_variance = var(log(num_genets/lag_genets)),
                         avg_genets = mean(num_genets),
                         med_genets = median(num_genets),
                         min_genets = min(num_genets),
                         max_genets = max(num_genets))
pop_growth_rate$expected_genets <- with(pop_growth_rate, (1+pop_growth_rate)/pop_grate_variance)

