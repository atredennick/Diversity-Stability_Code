#Look at species correlation in intrinsic yearly growth rates

rm(list=ls(all=TRUE))

library(corrplot)
library(reshape2)
library(ggplot2)
library(reshape2)
library(synchrony)
library(plyr)

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Get quadrat group information
quad_ids <- read.csv("../../../Data/Idaho/ARTR/growDnoNA.csv")
quad_groups <- unique(quad_ids[,c("quad","Group")])
quad_groups$Group <- as.numeric(quad_groups$Group)

all_data <- readRDS("RmaxYearly_byGroup_ipm.rds")
ts_metrics <- readRDS("../../../EmpiricalRelationships/quad_ts_metrics_AllSites.rds")
test <- as.data.frame(all_data)
synch <- numeric(6)
for(i in 1:6){
  tmp <- test[,c(i,i+6,i+12,i+18)]
  synch[i] <- community.sync(tmp)[[1]]
}
synch <- as.data.frame(synch)
synch$Group <- c(1:6)
colnames(synch)[1] <- "yr_synch"
synch <- merge(synch, quad_groups, by="Group")

# test$year <- seq(1,22)
# testM <- melt(test, id.vars = "year")
# testM$species <- c(rep("ARTR", 22*6),
#                    rep("HECO", 22*6),
#                    rep("POSE", 22*6),
#                    rep("PSSP", 22*6))
# testM$Group <- rep(c(1:6),each=length(unique(testM$year)))
# test_all <- merge(testM, quad_groups, by="Group")

ts_df <- as.data.frame(cbind(ts_metrics$Idaho$quad, ts_metrics$Idaho$stability))
colnames(ts_df) <- c("quad", "stability")
ts_synch <- merge(synch, ts_df, by="quad")
ts_out <- ddply(ts_synch, .(Group), summarise,
                stability = mean(as.numeric(stability)),
                yr_synch = mean(yr_synch))
ggplot(ts_out, aes(x=yr_synch, y=log(stability)))+
  geom_point(size=4)


allD <- read.csv("RmaxYearly_ipm.csv")
df <- melt(allD, id.vars = "yearParams")

synch_yr <- community.sync(allD[,1:4])[[1]]
ggplot(df, aes(x=yearParams, y=value, color=variable))+
#   geom_line(aes(size=variable))+
  geom_line()+
  geom_point(size=4)+
  geom_hline(aes(yintercept=0), color="grey59")+
  xlab("Year")+
  ylab("Yearly intrinsic growth rate")+
  scale_color_manual(name="Species", values=cbPalette[1:4])
#   scale_size_manual(values=c(2,1,1,1))
#   ggtitle(round(synch,2))

df_sum <- ddply(df, .(yearParams), summarise,
                comm = mean(value))
df_sum2<- ddply(subset(df, variable!="ARTR"), .(yearParams), summarise,
                comm = mean(value))

ggplot()+
#   geom_line(data=df, aes(x=yearParams, y=value, color=variable), alpha=0.2)+
#   geom_point(data=df, aes(x=yearParams, y=value, color=variable), size=4, alpha=0.2)+
  geom_line(data=df_sum, aes(x=yearParams, y=comm), size=3)+
  geom_line(data=df_sum2, aes(x=yearParams, y=comm), size=2, color="grey50", linetype=2)+
  xlab("Year")+
  ylab("Yearly intrinsic growth rate")+
  scale_color_manual(name="Species", values=cbPalette[1:4])

## Cycle through removing species and see which ones have biggest effect
synch_test <- numeric(4)
for(i in 1:4){
  ll <- which(c(1:4)!=i)
  synch_test[i] <- community.sync(allD[,ll])[[1]]
}



synch2 <- community.sync(allD[,2:4])[[1]]
g2 <- ggplot(subset(df, variable!="ARTR"), aes(x=yearParams, y=value, color=variable))+
  geom_line()+
  geom_point(size=4)+
  geom_hline(aes(yintercept=0), color="grey59")+
  xlab("Year")+
  ylab("Yearly intrinsic growth rate")+
  scale_color_manual(name="Species", values=cbPalette[2:4])+
  ggtitle(round(synch2,2))


# allD <- read.csv("RmaxYearly_ipm.csv")
# library(synchrony)
# c <- community.sync(allD[,1:4], nrands = 20000)

# M <- cor(allD[,1:4])
# corrplot.mixed(M, lower="ellipse", upper="number", col=c("tomato3","skyblue4"),tl.col="grey35")
# M
# 
# cov(allD[,1:4])
# 
# 
# matplot(allD$yearParams, as.matrix(allD[,1:4]), type="l")
# matplot(allD$yearParams, as.matrix(allD[,1:4]), type="p", pch=1, add=TRUE)
# synch <- sd(rowSums(allD[,1:4]))/(sd(allD[,1])+sd(allD[,2])+sd(allD[,3])+sd(allD[,4]))^2