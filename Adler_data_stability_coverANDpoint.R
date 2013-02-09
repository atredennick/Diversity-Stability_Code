rm(list=ls())

#set working directory
wd = "/Users/atredenn/Documents/Projects/Diversity_Stability/"
setwd(wd)
getwd()

library(plyr)

#Arizona data
pointdata.az <- read.csv("/Users/atredenn/Documents/Projects/Diversity_Stability/Arizona_Data/AZ_allrecords_point_features.csv")
pointdata.az <- as.data.frame(pointdata.az)
polydata.az <- read.csv("/Users/atredenn/Documents/Projects/Diversity_Stability/Arizona_Data/allrecords_polygon_features.csv")
polydata.az <- as.data.frame(polydata.az)

points <- subset(pointdata.az, select=c(quad, year, Species, Canopy_cov))
polys <- subset(polydata.az, select=c(quad, year, Species, Area))

colnames(points) <- colnames(polys)

data.az <- rbind(points, polys)
data.az <- data.az[order(data.az$quad, data.az$year),]
names(data.az)

data.az <- data.az[data.az$Area != 0,]

df.agg.az <- ddply(data.az, .(quad, year, Species), summarise, 
                   sum = sum(Area)
)

df.spp.az <- ddply(df.agg.az, .(quad, year), summarise,
                   sum = sum(sum),
                   spp = length(Species)
)

df.cv.az <- ddply(df.spp.az, .(quad), summarise,
                  cv = mean(sum)/sd(sum),
                  spp = mean(spp)
)

plot(df.cv.az$spp, df.cv.az$cv)
model <- lm(df.cv.az$cv ~ df.cv.az$spp)
summary(model)


spp <- as.character(unique(df.agg.az$Species))

par(mfrow=c(2,1))
df.A1P <- df.spp.az[df.spp.az$quad=="A1P",]
df.A1P.spp <- df.agg.az[df.agg.az$quad=="A1P",]
spp <- as.character(unique(df.A1P.spp$Species))
plot(df.A1P$year, df.A1P$sum, type="lines", col="red", lwd=1.6,
     ylab="cover (%)", xlab="year (19xx)")
points(df.A1P$year, df.A1P$sum, pch=21, bg="white", col="white", cex=1.5)
points(df.A1P$year, df.A1P$sum, pch=21, bg="white", col="red")
for(i in 1:length(spp)){
  lines(df.A1P.spp[df.A1P.spp$Species==spp[i],2], df.A1P.spp[df.A1P.spp$Species==spp[i],4],
        col="grey50", lwd=1.3)
  points(df.A1P.spp[df.A1P.spp$Species==spp[i],2], df.A1P.spp[df.A1P.spp$Species==spp[i],4],
         pch=21, col="white", bg="white", cex=1)
  points(df.A1P.spp[df.A1P.spp$Species==spp[i],2], df.A1P.spp[df.A1P.spp$Species==spp[i],4],
         col="grey50", cex=0.5)
}

df.WD5 <- df.spp.az[df.spp.az$quad=="WD5",]
df.WD5.spp <- df.agg.az[df.agg.az$quad=="WD5",]
spp <- as.character(unique(df.WD5.spp$Species))
plot(df.WD5$year, df.WD5$sum, type="lines", col="red", lwd=1.6, ylim=c(0,0.1),
     ylab="cover (%)", xlab="year (19xx)")
points(df.WD5$year, df.WD5$sum, pch=21, bg="white", col="white", cex=1.5)
points(df.WD5$year, df.WD5$sum, pch=21, bg="white", col="red")
for(i in 1:length(spp)){
  lines(df.WD5.spp[df.WD5.spp$Species==spp[i],2], df.WD5.spp[df.WD5.spp$Species==spp[i],4],
        col="grey50", lwd=1.3)
  points(df.WD5.spp[df.WD5.spp$Species==spp[i],2], df.WD5.spp[df.WD5.spp$Species==spp[i],4],
         pch=21, col="white", bg="white", cex=1)
  points(df.WD5.spp[df.WD5.spp$Species==spp[i],2], df.WD5.spp[df.WD5.spp$Species==spp[i],4],
         col="grey50", cex=0.5)
}

mod = lm (df.WD5$sum ~ df.WD5.spp[df.WD5.spp$Species==spp[1],4])
mod2 = lm (df.A1P$sum ~ df.A1P.spp[df.A1P.spp$Species==spp[1],4])

##Plot all quads with species and total cover
quads <- as.character(unique(df.spp.az$quad))

pdf("/users/atredenn/desktop/arizona_plots.pdf", height=5, width=5)
for(i in 1:length(quads)){
  df.quad <- df.spp.az[df.spp.az$quad==quads[i],]
  df.quad.spp <- df.agg.az[df.agg.az$quad==quads[i],]
  spp <- as.character(unique(df.quad.spp$Species))
  
  plot(df.quad$year, df.quad$sum, type="lines", col="red", lwd=1.6, ylim=c(0,0.1),
       ylab="cover (%)", xlab="year (19xx)", main=quads[i])
  points(df.quad$year, df.quad$sum, pch=21, bg="white", col="white", cex=1.5)
  points(df.quad$year, df.quad$sum, pch=21, bg="white", col="red")
  for(j in 1:length(spp)){
    lines(df.quad.spp[df.quad.spp$Species==spp[j],2], df.quad.spp[df.quad.spp$Species==spp[j],4],
          col="grey50", lwd=1.3)
    points(df.quad.spp[df.quad.spp$Species==spp[j],2], df.quad.spp[df.quad.spp$Species==spp[j],4],
           pch=21, col="white", bg="white", cex=1)
    points(df.quad.spp[df.quad.spp$Species==spp[j],2], df.quad.spp[df.quad.spp$Species==spp[j],4],
           col="grey50", cex=0.5)
  }
}
dev.off()
