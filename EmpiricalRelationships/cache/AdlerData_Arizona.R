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
polydata.az <- read.csv("/Users/atredenn/Documents/Projects/Diversity_Stability/Arizona_Data/allrecords_polygon_features.csv")
polydata.az <- as.data.frame(polydata.az)
names(polydata.az)
polydata.az$Species <- as.character(polydata.az$Species)


species.all <- read.csv("/Users/atredenn/Documents/Projects/Diversity_Stability/Arizona_Data/species_list.csv")
species.sort <- species.all[sort.list(species.all[,3]), ]
species.sort <- species.sort[!is.na(species.sort[,3]),]
index.doms <- seq(nrow(species.sort)-4, nrow(species.sort), 1)
dom.spp <- species.sort[index.doms,1]


polys.az <- subset(polydata.az, select=c(quad, year, Species, Area))

data.az <- polys.az
data.az <- data.az[order(data.az$quad, data.az$year),]
names(data.az)

data.arsp <- data.az[data.az$Species == dom.spp[1],]
data.hibe <- data.az[data.az$Species == dom.spp[2],]
data.boer <- data.az[data.az$Species == dom.spp[3],]
data.bope <- data.az[data.az$Species == dom.spp[4],]
data.boro <- data.az[data.az$Species == dom.spp[5],]


data.dom <- rbind(data.arsp, data.hibe, data.boer,
                  data.bope, data.boro)

data.dom <- data.dom[order(data.dom$quad, data.dom$year),]


#First sum areas over species within each quad for ALL species
df.agg.az <- ddply(data.az, .(quad, year, Species), summarise, 
                   sum = sum(Area),
                   numspecies = length(Species))
df.agg.tot <- ddply(data.az, .(quad, year), summarise, 
                    sum = sum(Area))
totavg.df <- ddply(df.agg.tot, .(year), summarise,
                   avg = mean(sum))

##################################################################################
##Break from aggregation for diversity-stability to plot dominants throught time##
df.agg.dom <- ddply(data.dom, .(quad, year, Species), summarise,
                    sum = sum(Area))

yrs <- sort(unique(df.agg.dom$year))
quads <- sort(unique(df.agg.dom$quad))
dat.all <- data.frame(quad=NA,
                       year=NA,
                       Species=NA,
                       sum=NA)
dat.quad <- data.frame(quad=NA,
                      year=NA,
                      Species=NA,
                      sum=NA)

for(i in 1:length(quads)){
  dat.quad <- df.agg.dom[df.agg.dom$quad == quads[i],]
  for (j in 1:length(yrs)){
    dat.quadyr <- dat.quad[dat.quad$year == yrs[j],]
    for (k in 1:length(dom.spp)){
      ifelse (length(grep(dom.spp[k], dat.quadyr$Species)) == 0,
              dat.quadyr <- rbind(dat.quadyr, data.frame(quad = quads[i],
                                                         year = yrs[j],
                                                         Species = dom.spp[k],
                                                         sum = 0)),
              dat.quadyr <- dat.quadyr)
    }
    dat.quad <- rbind(dat.quad, dat.quadyr)
  }
  dat.all <- rbind(dat.all, dat.quad)
}

df.agg.dom <- dat.all[2:nrow(dat.all),]

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

q1.df <- q1.df[q1.df[,1]!=47,]
df.agg.domavg <- df.agg.domavg[df.agg.domavg$year!=47,]
totavg.df <- totavg.df[totavg.df$year!=47,]

##Plot dominant species average cover through time
ggplot(data=q1.df, aes(x=year, y=(avg*100))) + 
  geom_line(aes(linetype=Species), color="grey45") +
  geom_point(aes(shape=Species), color="white", size=3) +
  geom_point(aes(shape=Species), size=2) +
  geom_line(aes(x=df.agg.domavg$year, y=(df.agg.domavg$tot*100)), color="steelblue", size=1.1) + 
  geom_line(aes(x=totavg.df$year, y=(totavg.df$avg*100)), color="darkorange", size=1.5) + 
  theme_bw() +
  xlab("Year (19xx)") + ylab("Mean Cover (%)") +
  scale_y_continuous(limits=c(0,10)) +
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
