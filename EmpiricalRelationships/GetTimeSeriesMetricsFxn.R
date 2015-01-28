##########################################################################
#### Function to extract time series for dominate species at a given site
#### and to calculate observed community synchrony and stability of cover.
####
#### Andrew Tredennick: atredenn@gmail.com
#### Date: 1-28-2015
####


get_comm_ts_metrics <- function(spp_list, site){
  num_spp <- length(spp_list)
  all_data <- data.frame(quad=NA, year=NA, totCover=NA, species=NA)
  for(dospp in 1:num_spp){ #loop through species to read in data
    spp_now <- spp_list[dospp]
    quad_file <- paste("../Data/", site,"/",spp_now,"/quadratCover.csv",sep="")
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
  nc <- round(sqrt(length(quad_list)))+1
  png(file = paste(site,"_allQuadratsTimeSeries.png"), width = 1500, height = 1000)
  do.call(grid.arrange, c(plots$plots, nrow=nc, ncol=nc))
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
  png(paste(site,"_quadratCommunitySynchrony.png"), width = 500, height = 500)
  ggplot(quadrat_synchrony)+
    geom_histogram(aes(x=community_synchrony), color="white", binwidth=0.1)+
    ggtitle("Quadrat Community Synchrony")+
    xlab(expression(phi[C]))+
    ylab("Number of quadrats")
  dev.off()
#   write.csv(quadrat_synchrony, file=paste(site,"_quadratCommunitySynchrony.csv"))
  
  #Calculate temporal stability (mean/sd; Tilman references) for each quadrat
  quad_ts <- numeric(length(quad_list))
  for(q in 1:length(quad_list)){
    tmp_quad_data <- subset(all_data, quad==quad_list[q])
    community_cover <-  ddply(tmp_quad_data, c("year"), .fun = summarise,
                              all_cover = sum(totCover))
    community_ts <- mean(community_cover$all_cover)/sd(community_cover$all_cover)
    quad_ts[q] <- community_ts
  }
  #Save plot and dataframe
  quadrat_stability <- as.data.frame(quad_ts)
  quadrat_stability$quad <- quad_list
  colnames(quadrat_stability)[1] <- "temporal_stability"
  
#   write.csv(quadrat_stability, file=paste(site,"_quadratTemporalStability.csv"))
  
  community_metrics <- merge(quadrat_stability, quadrat_synchrony, by="quad")
  png(paste(site,"_synchrony_vs_stability.png"), width = 500, height = 500)
  ggplot(community_metrics, aes(x=community_synchrony, y=log(temporal_stability)))+
    geom_point()+
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), se=FALSE)+
    stat_smooth(method = "lm", formula = y ~ x, color="darkorange", se=FALSE)+
    xlab(expression(phi[C]))+
    ylab("ln(TS)")
  dev.off()
  
  metrics <- list(quad=quad_list,
                  stability=quadrat_stability$temporal_stability,
                  comm_synchrony=quadrat_synchrony$community_synchrony)
  return(metrics)
}#end of function

