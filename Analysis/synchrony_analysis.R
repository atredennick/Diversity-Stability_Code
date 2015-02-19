#########################################################
#### Script to calculate synchrony metrics from observed
#### data and simulated growth rates
####
#### Andrew Tredennick: atredenn@gmail.com
#### 2-17-2015

## Clear the workspace
rm(list=ls(all=TRUE))

####
#### 0.1. Load libraries ----------------------------------------
####
library(reshape2); library(ggplot2)
library(synchrony); library(plyr); library(tidyr)

####
#### 0.2. Set color-blind pallete -------------------------------
####
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

####
#### 1. Observed community synchrony from time series (population cover) --------
####
# Average over quadrats first, ignoring certain one's for Kansas
# due to over-dominance by BOCU. These excluded quadrats are taken from
# Chnegjin's code. Results in 7 quadrats for the data.

## Kansas
# bring in data
site <- "Kansas"
spp_list <- c("BOCU","BOHI","SCSC")
num_spp <- length(spp_list)
ks_data <- data.frame(quad=NA, year=NA, totCover=NA, species=NA)
for(dospp in 1:num_spp){ #loop through species to read in data
  spp_now <- spp_list[dospp]
  quad_file <- paste("../Data/", site,"/",spp_now,"/quadratCover.csv",sep="")
  spp_data <- read.csv(quad_file)
  spp_data$species <- spp_now
  ks_data <- rbind(ks_data, spp_data)
} #end species looping for raw data
ks_data <- ks_data[2:nrow(ks_data),] #remove first NA row

# now exclude some quadrats (from Chengjin's work)
tmp1<-which(ks_data$quad_data=="q25" & (ks_data$year<35 | ks_data$year>62))
tmp2<-which(ks_data$quad_data=="q27")
tmp3<-which(ks_data$quad=="q28")
tmp4<-which(ks_data$quad=="q30")
tmp5<-which(ks_data$quad=="q31" & (ks_data$year<35 | ks_data$year>39))
tmp6<-which(ks_data$quad=="q32" & (ks_data$year<35 | ks_data$year>41))
tmp<-c(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)
ks_data<-ks_data[-tmp,]
##excluding the records later than 1968, to keep the same random year effect...
ks_data<-subset(ks_data,year<68)

# aggregate the quadrat-level data for average observed cover
# divide totCover by 100 to convert from m2 to percent cover in 1m2 plot
ks_agg <- ddply(ks_data, .(year, species), summarise,
                tot_cover = mean(totCover/100))

# caclulate synchrony and stability
# first stability (sd/mean) by population and then community
ts_ks <- numeric(num_spp+1)
for(i in 1:num_spp){
  tmp <- subset(ks_agg, species==spp_list[i])
  ts_ks[i] <- sd(tmp$tot_cover)/mean(tmp$tot_cover)
}
ts_ks[num_spp+1] <- sd(ks_agg$tot_cover)/mean(ks_agg$tot_cover)

# now synchrony of the population fluctuations
ks_mat <- dcast(ks_agg, formula = year~species)
synch_ks <- community.sync(ks_mat[2:4])

# now synchrony of yearly intrinsic growth rates
ks_yrpgr <- readRDS("../PopModels/kansasIPM/simulations/RmaxYearly_ipm.rds")
synch_ks_yr <- community.sync(ks_yrpgr)

# caclulate observed growth rates
transitions <- nrow(ks_mat)-1
ks_obs_gr <- matrix(nrow=transitions, ncol=num_spp)
for(i in 1:transitions){
  ks_obs_gr[i,] <- as.numeric(log(ks_mat[i+1,2:4]/ks_mat[i,2:4]))
}

ks_growthrates <- cbind(melt(ks_yrpgr[1:33,]), melt(ks_obs_gr)$value)
colnames(ks_growthrates) <- c("year", "species", "simyrpgr", "obsyrgr")
growth_plot_ks <- ggplot(ks_growthrates, aes(x=simyrpgr, y=obsyrgr, color=species))+
                    geom_point(size=3)+
                    stat_smooth(method="lm", se=FALSE, size=1)+
                    scale_color_manual(values=cbPalette[1:num_spp], name="")+
                    xlab("Environmental response (simulated)")+
                    ylab("Yearly growth rate (observed)")

# get R2 from linear models for percent of observed transition
# explained by species response to environment
r_squares <- numeric(num_spp)
for(i in 1:num_spp){
  r_squares[i] <- summary(lm(obsyrgr~simyrpgr, 
                             data=subset(ks_growthrates, species==spp_list[i]))
                          )$r.squared
}
explained_var_ks <- data.frame("species" = spp_list,
                            "var_exp" = r_squares*100)


## Idaho
# bring in data
site <- "Idaho"
spp_list <- sort(c("PSSP","HECO","POSE","ARTR"))
num_spp <- length(spp_list)
id_data <- data.frame(quad=NA, year=NA, totCover=NA, species=NA)
for(dospp in 1:num_spp){ #loop through species to read in data
  spp_now <- spp_list[dospp]
  quad_file <- paste("../Data/", site,"/",spp_now,"/quadratCover.csv",sep="")
  spp_data <- read.csv(quad_file)
  spp_data$species <- spp_now
  id_data <- rbind(id_data, spp_data)
} #end species looping for raw data
id_data <- id_data[2:nrow(id_data),] #remove first NA row

# aggregate the quadrat-level data for average observed cover
# divide totCover by 100 to convert from m2 to percent cover in 1m2 plot
id_agg <- ddply(id_data, .(year, species), summarise,
                tot_cover = mean(totCover/100))

# caclulate synchrony and stability
# first stability (sd/mean) by population and then community
ts_id <- numeric(num_spp+1)
for(i in 1:num_spp){
  tmp <- subset(id_agg, species==spp_list[i])
  ts_id[i] <- sd(tmp$tot_cover)/mean(tmp$tot_cover)
}
ts_id[num_spp+1] <- sd(id_agg$tot_cover)/mean(id_agg$tot_cover)

# now synchrony of the population fluctuations
id_mat <- dcast(id_agg, formula = year~species)
synch_id <- community.sync(id_mat[2:4])

# now synchrony of yearly intrinsic growth rates
id_yrpgr <- readRDS("../PopModels/idahoIPM/simulations/RmaxYearly_ipm.rds")
synch_id_yr <- community.sync(id_yrpgr)

# caclulate observed growth rates
transitions <- nrow(id_mat)-1
id_obs_gr <- matrix(nrow=transitions, ncol=num_spp)
for(i in 1:transitions){
  id_obs_gr[i,] <- as.numeric(log(id_mat[i,2:5]/id_mat[i+1,2:5]))
}

id_growthrates <- cbind(melt(id_yrpgr[1:nrow(id_obs_gr),]), melt(id_obs_gr)$value)
colnames(id_growthrates) <- c("year", "species", "simyrpgr", "obsyrgr")
ggplot(id_growthrates, aes(x=simyrpgr, y=obsyrgr, color=species))+
  geom_point()+
  stat_smooth(method="lm", se=FALSE)+
  scale_color_manual(values=cbPalette[1:num_spp])

#Kansas time series
ggplot()+
  geom_line(data=ks_agg, aes(x=year, y=tot_cover, color=species))+
  geom_point(data=ks_agg, aes(x=year, y=tot_cover, color=species), size=3)+
  scale_color_manual(values=cbPalette[1:num_spp])


# 
# #Get quadrat group information
# quad_ids <- read.csv("../../../Data/Idaho/ARTR/growDnoNA.csv")
# quad_groups <- unique(quad_ids[,c("quad","Group")])
# quad_groups$Group <- as.numeric(quad_groups$Group)
# 
# all_data <- readRDS("RmaxYearly_byGroup_ipm.rds")
# ts_metrics <- readRDS("../../../EmpiricalRelationships/quad_ts_metrics_AllSites.rds")
# test <- as.data.frame(all_data)
# synch <- numeric(6)
# for(i in 1:6){
#   tmp <- test[,c(i,i+6,i+12,i+18)]
#   synch[i] <- community.sync(tmp)[[1]]
# }
# synch <- as.data.frame(synch)
# synch$Group <- c(1:6)
# colnames(synch)[1] <- "yr_synch"
# synch <- merge(synch, quad_groups, by="Group")
# 
# # test$year <- seq(1,22)
# # testM <- melt(test, id.vars = "year")
# # testM$species <- c(rep("ARTR", 22*6),
# #                    rep("HECO", 22*6),
# #                    rep("POSE", 22*6),
# #                    rep("PSSP", 22*6))
# # testM$Group <- rep(c(1:6),each=length(unique(testM$year)))
# # test_all <- merge(testM, quad_groups, by="Group")
# 
# ts_df <- as.data.frame(cbind(ts_metrics$Idaho$quad, ts_metrics$Idaho$stability))
# colnames(ts_df) <- c("quad", "stability")
# ts_synch <- merge(synch, ts_df, by="quad")
# ts_out <- ddply(ts_synch, .(Group), summarise,
#                 stability = mean(as.numeric(stability)),
#                 yr_synch = mean(yr_synch))
# ggplot(ts_out, aes(x=yr_synch, y=log(stability)))+
#   geom_point(size=4)
# 
# 
# allD <- read.csv("RmaxYearly_ipm.csv")
# df <- melt(allD, id.vars = "yearParams")
# 
# synch_yr <- community.sync(allD[,1:4])[[1]]
# ggplot(df, aes(x=yearParams, y=value, color=variable))+
# #   geom_line(aes(size=variable))+
#   geom_line()+
#   geom_point(size=4)+
#   geom_hline(aes(yintercept=0), color="grey59")+
#   xlab("Year")+
#   ylab("Yearly intrinsic growth rate")+
#   scale_color_manual(name="Species", values=cbPalette[1:4])
# #   scale_size_manual(values=c(2,1,1,1))
# #   ggtitle(round(synch,2))
# 
# df_sum <- ddply(df, .(yearParams), summarise,
#                 comm = mean(value))
# df_sum2<- ddply(subset(df, variable!="ARTR"), .(yearParams), summarise,
#                 comm = mean(value))
# 
# ggplot()+
# #   geom_line(data=df, aes(x=yearParams, y=value, color=variable), alpha=0.2)+
# #   geom_point(data=df, aes(x=yearParams, y=value, color=variable), size=4, alpha=0.2)+
#   geom_line(data=df_sum, aes(x=yearParams, y=comm), size=3)+
#   geom_line(data=df_sum2, aes(x=yearParams, y=comm), size=2, color="grey50", linetype=2)+
#   xlab("Year")+
#   ylab("Yearly intrinsic growth rate")+
#   scale_color_manual(name="Species", values=cbPalette[1:4])
# 
# ## Cycle through removing species and see which ones have biggest effect
# synch_test <- numeric(4)
# for(i in 1:4){
#   ll <- which(c(1:4)!=i)
#   synch_test[i] <- community.sync(allD[,ll])[[1]]
# }
# 
# 
# 
# synch2 <- community.sync(allD[,2:4])[[1]]
# g2 <- ggplot(subset(df, variable!="ARTR"), aes(x=yearParams, y=value, color=variable))+
#   geom_line()+
#   geom_point(size=4)+
#   geom_hline(aes(yintercept=0), color="grey59")+
#   xlab("Year")+
#   ylab("Yearly intrinsic growth rate")+
#   scale_color_manual(name="Species", values=cbPalette[2:4])+
#   ggtitle(round(synch2,2))
# 
# 
# # allD <- read.csv("RmaxYearly_ipm.csv")
# # library(synchrony)
# # c <- community.sync(allD[,1:4], nrands = 20000)
# 
# # M <- cor(allD[,1:4])
# # corrplot.mixed(M, lower="ellipse", upper="number", col=c("tomato3","skyblue4"),tl.col="grey35")
# # M
# # 
# # cov(allD[,1:4])
# # 
# # 
# # matplot(allD$yearParams, as.matrix(allD[,1:4]), type="l")
# # matplot(allD$yearParams, as.matrix(allD[,1:4]), type="p", pch=1, add=TRUE)
# # synch <- sd(rowSums(allD[,1:4]))/(sd(allD[,1])+sd(allD[,2])+sd(allD[,3])+sd(allD[,4]))^2