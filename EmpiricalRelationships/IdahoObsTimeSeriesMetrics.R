####################################################################
#### Script to plot observed Idaho time series for dominate species
#### and to calculate observed community synchrony in dynamics.
####
#### Andrew Tredennick: atredenn@gmail.com
#### Date began: 1-28-2015
####
#### Make sure working directory is set to source file location.

#clear the workspace
rm(list=ls(all=TRUE))

####
#### Load libraries -----------------------------------------------
####
library(ggplot2); library(plyr); library(reshape2); library(gridExtra)


####
#### Bring in the raw data -----------------------------------------
####
spp_list <- c("ARTR","HECO","POSE","PSSP")
num_spp <- length(spp_list)
all_data <- data.frame(quad=NA, year=NA, totCover=NA, species=NA)

for(dospp in 1:num_spp){ #loop through species to read in data
  spp_now <- spp_list[dospp]
  quad_file <- paste("../PopModels/IdahoIPM/speciesData/",spp_now,"/quadratCover.csv",sep="")
  spp_data <- read.csv(quad_file)
  spp_data$species <- spp_now
  all_data <- rbind(all_data, spp_data)
} #end species looping for raw data

all_data <- all_data[2:nrow(all_data),] #remove first NA row


####
#### Plot individual quadrat time series ----------------------------
####
quad_list <- unique(all_data$quad)
plot_list <- list()
for(q in 1:length(quad_list)){
  tmp_quad_data <- subset(all_data, quad==quad_list[q])
  tmp_plot <- ggplot(tmp_quad_data, aes(x=year, y=totCover))+
                geom_line(aes(linetype=species, color=species))+
                geom_point(aes(shape=species, color=species))+
                guides(linetype=FALSE, shape=FALSE, color=FALSE)+
                ggtitle(paste("Quad", quad_list[q]))
  plot_list <- c(plot_list, list(tmp_plot))
}
plots <- list(plots=plot_list)
png(file = "allQuadratsTimeSeries.png", width = 1500, height = 1000)
do.call(grid.arrange, c(plots$plots, nrow=5, ncol=6))
dev.off()


####
#### Calculate time series metrics for individual quadrats ----------------
####
#Calculate community synchrony as variance of community divided by the squared sum
#of population standard deviations (Loreau and de Mazancourt 2008, Ecol Letts).
quad_synchs <- numeric(length(quad_list))
for(q in 1:length(quad_list)){
  tmp_quad_data <- subset(all_data, quad==quad_list[q])
  community_cover <-  ddply(tmp_quad_data, c("year"), .fun = summarise,
                            all_cover = sum(totCover))
  community_variance <- var(community_cover$all_cover)
  population_stdev <- ddply(tmp_quad_data, c("species"), .fun = summarise,
                            stdev_cover = sd(totCover))
  sum_pop_stdev <- sum(population_stdev$stdev, na.rm = TRUE)
  community_synchrony <- community_variance/(sum_pop_stdev^2)
  quad_synchs[q] <- community_synchrony
}
quadrat_synchrony <- as.data.frame(quad_synchs)
quadrat_synchrony$quad <- quad_list
colnames(quadrat_synchrony)[1] <- "community_synchrony"

#Save plot and data frame
png("quadratCommunitySynchrony.png", width = 500, height = 500)
ggplot(quadrat_synchrony)+
  geom_histogram(aes(x=community_synchrony), color="white", binwidth=0.1)+
  ggtitle("Quadrat Community Synchrony")+
  xlab(expression(phi[C]))+
  ylab("Number of quadrats")
dev.off()
write.csv(quadrat_synchrony, file="quadratCommunitySynchrony.csv")



