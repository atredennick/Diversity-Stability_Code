rm(list=ls())

#set working directory
wd = "/Users/atredenn/Documents/Projects/Diversity_Stability/"
setwd(wd)
getwd()

library(plyr)
library(ggplot2)
library(grid)

### Idaho data ###

#read in data from point and polygon files
polydata.mt <- read.csv("/Users/atredenn/Documents/Projects/Diversity_Stability/Montana_Data/Montana_allrecords_cover.csv")
polydata.mt <- as.data.frame(polydata.mt)
polydata.mt$Species <- as.character(polydata.mt$Species)

pointdata.mt <- read.csv("/Users/atredenn/Documents/Projects/Diversity_Stability/Montana_Data/Montana_allrecords_density.csv")

data.pasm <- pointdata.mt[pointdata.mt$Species == "Pascopyrum smithii",]
data.pasm$Area <- 0.000025

data.pasm2 <- subset(data.pasm, select=c(quad, year, Species, Area))


polys.mt <- subset(polydata.mt, select=c(quad, year, Species, Area))

data.mt <- polys.mt
data.mt <- data.mt[order(data.mt$quad, data.mt$year),]
names(data.mt)

#Separate out dominant species, skip this if not interested in just dominant species
dom.spp <- c("Bouteloua gracilis",
             "Pascopyrum smithii",
             "Hesperostipa comata",
             "Carex filifolia",
             "Poa secunda",
             "Opuntia polyacantha")

data.bogr <- data.mt[data.mt$Species == dom.spp[1],]
data.pasm <- data.mt[data.mt$Species == dom.spp[2],]
data.heco <- data.mt[data.mt$Species == dom.spp[3],]
data.cafi <- data.mt[data.mt$Species == dom.spp[4],]
data.pose <- data.mt[data.mt$Species == dom.spp[5],]
data.oppo <- data.mt[data.mt$Species == dom.spp[6],]

data.dom <- rbind(data.bogr, data.pasm, data.heco,
                  data.cafi, data.pose, data.oppo, data.pasm2)

data.dom <- data.dom[order(data.dom$quad, data.dom$year),]


#First sum areas over species within each quad for ALL species
df.agg.mt <- ddply(data.mt, .(quad, year, Species), summarise, 
                   sum = sum(Area),
                   numspecies = length(Species))
df.agg.tot <- ddply(data.mt, .(quad, year), summarise, 
                    sum = sum(Area))
totavg.df <- ddply(df.agg.tot, .(year), summarise,
                   avg = mean(sum))

##################################################################################
##Break from aggregation for diversity-stability to plot dominants throught time##
df.agg.dom <- ddply(data.dom, .(quad, year, Species), summarise,
                    sum = sum(Area))
q1.df <- ddply(df.agg.dom, .(year, Species), summarise, 
               avg = mean(sum))


df.agg.domtot <- ddply(data.dom, .(quad, year), summarise,
                       sum = sum(Area))
df.agg.domavg <- ddply(df.agg.domtot, .(year), summarise,
                       tot = mean(sum))

ddd <- data.dom[data.dom$year == 32,]
sum(ddd$area)


# 
# 
# yrs <- seq(32,45,1)
# dat.all <- data.frame(year = NA,
#                       Species = NA,
#                       avg = NA)
# 
# for (i in 1:length(yrs)){
#   dat <- q1.df[q1.df$year == yrs[i],]
#   for (j in 1:7){
#     ifelse (length(grep(dom.spp[j], dat$Species)) == 0,
#             dat <- rbind(dat, data.frame(year = yrs[i],
#                                          Species = dom.spp[j],
#                                          avg = 0)),
#             dat <- dat)
#   }
#   dat.all <- rbind(dat.all, dat)
# }  
# 
# q1.df <- dat.all[2:nrow(dat.all),]
# 

##Plot dominant species average cover through time
ggplot(data=q1.df, aes(x=year, y=(avg*100))) + 
  geom_line(aes(linetype=Species), color="grey45") +
  geom_point(aes(shape=Species), color="white", size=3) +
  geom_point(aes(shape=Species), size=2) +
  geom_line(aes(x=df.agg.domavg$year, y=(df.agg.domavg$tot*100)), color="steelblue", size=1.1) + 
  geom_line(aes(x=totavg.df$year, y=(totavg.df$avg*100)), color="darkorange", size=1.5) + 
  theme_bw() +
  xlab("Year (19xx)") + ylab("Mean Cover (%)") +
  scale_y_continuous(limits=c(0,40)) +
  theme(axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14, angle=90), 
        axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=12), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.75,0.85),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text=element_text(size=8, face="italic"),
        legend.background = element_rect(colour = NA),
        legend.key.height = unit(1, "mm")
  ) 
