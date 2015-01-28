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
pointdata.id <- read.csv("/Users/atredenn/Documents/Projects/Diversity_Stability/Idaho_Data/Idaho_allrecords_density.csv")
pointdata.id <- as.data.frame(pointdata.id)
polydata.id <- read.csv("/Users/atredenn/Documents/Projects/Diversity_Stability/Idaho_Data/Idaho_allrecords_cover.csv")
polydata.id <- as.data.frame(polydata.id)
polydata.id$species <-as.character(polydata.id$species)
species.id <- read.csv("/Users/atredenn/Documents/Projects/Diversity_Stability/Idaho_Data/species_list.csv")
species.id <- as.data.frame(species.id)
species.id$species <- as.character(species.id$species)
species.id$growthForm <- as.character(species.id$growthForm)

#remove duplicates bewteen point and poly data
pointdata.id <- pointdata.id[pointdata.id$stem=="N",]

#subset the data for columns of interest
points.id <- subset(pointdata.id, select=c(quad, year, species))
polys.id <- subset(polydata.id, select=c(quad, year, species, area))

polys.id$growthForm <- character(length(polys.id[,1]))

for (i in 1:length(species.id[,1])){
 index <- which(polys.id[,3] == species.id[i,1])
 polys.id[index,5] <- species.id[i,4]
}

#This section commented out for now, until better idea of how to deal with point data.
###############################################
#combine polygon and point data
#first add in "area" column for points data
# points.id$area <- 0
# data.id <- rbind(points.id, polys.id)
##############################################

data.id <- polys.id
data.id <- data.id[order(data.id$quad, data.id$year),]
names(data.id)

data.id <- data.id[data.id$year != 73,]

#Separate out dominant species, skip this if not interested in just dominant species
dom.spp <- c("Poa secunda",
             "Artemisia tripartita",
             "Hesperostipa comata",
             "Pseudoroegneria spicata")

data.poa <- data.id[data.id$species == dom.spp[1],]
data.art <- data.id[data.id$species == dom.spp[2],]
data.hes <- data.id[data.id$species == dom.spp[3],]
data.pse <- data.id[data.id$species == dom.spp[4],]

data.dom <- rbind(data.poa, data.art, data.hes, data.pse)

data.dom <- data.dom[order(data.dom$quad, data.dom$year),]


#First sum areas over species within each quad for ALL species
df.agg.id <- ddply(data.id, .(quad, year, species), summarise, 
                   sum = sum(area),
                   numspecies = length(species))
df.agg.tot <- ddply(data.id, .(quad, year), summarise, 
                   sum = sum(area))
totavg.df <- ddply(df.agg.tot, .(year), summarise,
                   avg = (mean(sum)*100))

##################################################################################
##Break from aggregation for diversity-stability to plot dominants throught time##
df.agg.dom <- ddply(data.dom, .(quad, year, species), summarise,
                    sum = sum(area))
yrs <- sort(unique(df.agg.dom$year))
quads <- sort(unique(df.agg.dom$quad))
dat.all <- data.frame(quad=NA,
                      year=NA,
                      species=NA,
                      sum=NA)
dat.quad <- data.frame(quad=NA,
                       year=NA,
                       species=NA,
                       sum=NA)

for(i in 1:length(quads)){
  dat.quad <- df.agg.dom[df.agg.dom$quad == quads[i],]
  for (j in 1:length(yrs)){
    dat.quadyr <- dat.quad[dat.quad$year == yrs[j],]
    for (k in 1:length(dom.spp)){
      ifelse (length(grep(dom.spp[k], dat.quadyr$species)) == 0,
              dat.quadyr <- rbind(dat.quadyr, data.frame(quad = quads[i],
                                                         year = yrs[j],
                                                         species = dom.spp[k],
                                                         sum = 0)),
              dat.quadyr <- dat.quadyr)
    }
    dat.quad <- rbind(dat.quad, dat.quadyr)
  }
  dat.all <- rbind(dat.all, dat.quad)
}

df.agg.dom <- dat.all[2:nrow(dat.all),]
q1.df <- ddply(df.agg.dom, .(year, species), summarise, 
               avg = mean(sum))


df.agg.domtot <- ddply(data.dom, .(quad, year), summarise,
                    sum = sum(area))
df.agg.domavg <- ddply(df.agg.domtot, .(year), summarise,
                       tot = (mean(sum)*100))

##Plot dominant species average cover through time
ggplot(data=q1.df, aes(x=year, y=(avg*100))) + 
  geom_line(aes(linetype=species), color="grey45") +
#   geom_point(aes(shape=species), color="white", size=3) +
  geom_point(aes(shape=species), size=2) +
#   geom_point(aes(color=species), size=2) +
#   stat_smooth(aes(fill=species), color="white", size=1) +
  geom_line(aes(x=df.agg.domavg$year, y=df.agg.domavg$tot), color="steelblue", size=1.1) + 
  geom_line(aes(x=totavg.df$year, y=totavg.df$avg), color="darkorange", size=1.5) + 
  theme_bw() +
  xlab("Year (19xx)") + ylab("Mean Cover (%)") +
  theme(axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14, angle=90), 
        axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=12), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.25,0.85),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text=element_text(size=8, face="italic"),
        legend.background = element_rect(colour = NA),
        legend.key.height = unit(1, "mm")
  ) 

##################################################################################

#Second, sum for total community area each year in each quad, and calculate species richness
df.spp.id <- ddply(df.agg.id, .(quad, year), summarise,
                   sum = sum(sum),
                   spp = length(species),
                   quadspp = sum(numspecies)
)

#Subsection of code for setting up calculation of Simpson's diversity
quad.names <- unique(df.spp.id$quad)
quadyears.simpson <- matrix(ncol=length(quad.names), nrow=30)
for (i in 1:length(quad.names)){
  index <- which(df.agg.id$quad == quad.names[i])
  data.now <- df.agg.id[index,]
  years <- unique(data.now$year)
  
  for (j in 1:length(years)){
    index.year <- which(data.now$year == years[j])
    data.year <- data.now[index.year,]
    
    simps.start <- data.year$numspecies/sum(data.year$numspecies) 
    quadyears.simpson[j,i] <- sum(simps.start^2)
  }
}

quadyears.simpson <- c(quadyears.simpson)
quadyears.simpson <- quadyears.simpson[!is.na(quadyears.simpson)]
df.spp.id$simpson <- quadyears.simpson


#Lastly, caclulate temporal stability for each quad and average richness through time
df.cv.id <- ddply(df.spp.id, .(quad), summarise,
                  cv = sd(sum)/mean(sum),
                  spp = mean(spp),
                  simp = mean(1-simpson)
)

ggplot(data=df.cv.id) +
  geom_histogram(aes(x=cv), binwidth=0.07, fill="grey", color="white") +
  geom_density(aes(x=cv), fill="black", alpha=0.4) +
  theme_bw() +
  xlab("Coefficient of Variation") + ylab("Frequency (# of plots)") +
  theme(axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14, angle=90), 
        axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=12), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.25,0.85),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text=element_text(size=8, face="italic"),
        legend.background = element_rect(colour = NA),
        legend.key.height = unit(1, "mm")
  ) 
  


#Plot density of observed CV
hist(df.cv.id$cv, freq=FALSE, col="grey70", border="white", 
     xlab="Temporal C.V.", main="", breaks=4, xlim=c(0,1))
lines(density(df.cv.id$cv, adjust=1), lwd=4, col="black")
# text(x=1.5, y=c(1.7,1.5,1.3,1.1,0.9), c("ARTR", "HECO", "POSE", "PSSP", "rare"), cex=0.8)
# text(x=1.5, y=1.9, "Species", cex=1.2)
# lines(x=rep(median(df.cv.id$cv),3), y=seq(0,2), lty="dashed", lwd=2, col="white")
box()
