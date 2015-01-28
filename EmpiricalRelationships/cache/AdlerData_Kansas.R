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
polydata.ks <- read.csv("/Users/atredenn/Documents/Projects/Diversity_Stability/Kansas_Data/KansasData_Adler_Clean.csv")
polydata.ks <- as.data.frame(polydata.ks)
polydata.ks$species <-as.character(polydata.ks$species)

polys.ks <- subset(polydata.ks, select=c(quad, year, species, area))

data.ks <- polys.ks
data.ks <- data.ks[order(data.ks$quad, data.ks$year),]
names(data.ks)

data.ks$graze <- NA
index.graze1 <- grep("e2qa", data.ks$quad)
index.graze2 <- grep("e2qo", data.ks$quad)
index.graze <- c(index.graze1, index.graze2)
data.ks[index.graze,5] <- 2
data.ks[is.na(data.ks$graze),5] <- 1

data.ks <- data.ks[data.ks[,5] == 1,]
data.ks <- subset(data.ks, select=c(quad, year, species, area))


#Separate out dominant species, skip this if not interested in just dominant species
dom.spp <- c("Bouteloua gracilis",
             "Buchloe dactyloides",
             "Schizachyrium scoparium",
             "Bouteloua curtipendula",
             "Bouteloua hirsuta",
             "Andropogon gerardii")

data.bogr <- data.ks[data.ks$species == dom.spp[1],]
data.buda <- data.ks[data.ks$species == dom.spp[2],]
data.scsc <- data.ks[data.ks$species == dom.spp[3],]
data.bocu <- data.ks[data.ks$species == dom.spp[4],]
data.bohi <- data.ks[data.ks$species == dom.spp[5],]
data.ange <- data.ks[data.ks$species == dom.spp[6],]

data.dom <- rbind(data.bogr, data.buda, data.scsc, data.bocu, data.bohi, data.ange)

data.dom <- data.dom[order(data.dom$quad, data.dom$year),]


#First sum areas over species within each quad for ALL species
df.agg.ks <- ddply(data.ks, .(quad, year, species), summarise, 
                   sum = sum(area),
                   numspecies = length(species))
df.agg.tot <- ddply(data.ks, .(quad, year), summarise, 
                    sum = sum(area))
totavg.df <- ddply(df.agg.tot, .(year), summarise,
                   avg = mean(sum))

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
                       tot = mean(sum))

ddd <- data.dom[data.dom$year == 32,]
sum(ddd$area)




yrs <- seq(32,72,1)
dat.all <- data.frame(year = NA,
                      species = NA,
                      avg = NA)

for (i in 1:length(yrs)){
  dat <- q1.df[q1.df$year == yrs[i],]
  for (j in 1:6){
    ifelse (length(grep(dom.spp[j], dat$species)) == 0,
            dat <- rbind(dat, data.frame(year = yrs[i],
                                         species = dom.spp[j],
                                         avg = 0)),
            dat <- dat)
  }
  dat.all <- rbind(dat.all, dat)
 }  

q1.df <- dat.all[2:nrow(dat.all),]


##Plot dominant species average cover through time
ggplot(data=q1.df, aes(x=year, y=(avg/100))) + 
  geom_line(aes(linetype=species), color="grey45") +
  geom_point(aes(shape=species), color="white", size=3) +
  geom_point(aes(shape=species), size=2) +
  geom_line(aes(x=df.agg.domavg$year, y=(df.agg.domavg$tot/100)), color="steelblue", size=1.1) + 
  geom_line(aes(x=totavg.df$year, y=(totavg.df$avg/100)), color="darkorange", size=1.5) + 
  theme_bw() +
  xlab("Year (19xx)") + ylab("Mean Cover (%)") +
  scale_y_continuous(limits=c(0,100)) +
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
