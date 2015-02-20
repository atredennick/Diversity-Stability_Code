###############################################
#### Climate fluctuation and trend analysis
####
#### Andrew Tredennick: atredenn@gmail.com

## Clear the workspace
rm(list=ls(all=TRUE))

####
#### 0.1. Load libraries ----------------------------------------
####
library(reshape2); library(ggplot2); library(gridExtra)
library(plyr); library(tidyr)

####
#### 1. Look at trend and fluctuations in Idaho climate data ---------------
####
idaho_climate <- read.csv("../Data/Climate/ClimateIdaho.csv")
idaho_melt <- melt(idaho_climate, id.vars = "year")
ggplot(subset(idaho_melt, variable=="ppt1"|variable=="TmeanSpr1"), aes(x=year, y=value))+
  geom_line(aes(color=variable), alpha=0.5)+
  geom_point(aes(color=variable), size=4, alpha=0.5)+
  stat_smooth(method="lm", se=FALSE, aes(color=variable), size=1)+
  facet_grid(variable~., scales = "free_y")+
  scale_color_manual(values=c("dodgerblue", "darkred"))+
  guides(color=FALSE)




