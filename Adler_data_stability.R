rm(list=ls())

#set working directory
wd = "/Users/atredenn/Documents/Projects/Diversity_Stability/"
setwd(wd)
getwd()

library(plyr)

#Arizona data
data.az <- read.csv("/Users/atredenn/Documents/Projects/Diversity_Stability/AZ_allrecords_point_features.csv")
data.az <- as.data.frame(data.az)

df.agg.az <- ddply(data.az, .(quad, year, Species), summarise, 
                sum = length(Species)
               )

df.spp.az <- ddply(df.agg.az, .(quad, year), summarise,
                sum = sum(sum),
                spp = length(Species)
                )

df.cv.az <- ddply(df.spp.az, .(quad), summarise,
               cv = mean(sum)/sd(sum),
               spp = mean(spp),
               mu = mean(sum)
              )
df.cv.az <- df.cv.az[!is.na(df.cv.az$cv),]
df.cv.az <- df.cv.az[!is.infinite(df.cv.az$cv),]

# y.az <- df.cv.az$cv
# qqnorm(y.az)
y.az <- log10(df.cv.az$cv)
x.az <- df.cv.az$spp
model.az <- lm(y.az~x.az)
summary(model.az)

#create prediction vector
newx.az=seq(1,3.5, 0.01)
prd.az <- predict(model.az, newdata=data.frame(x.az=newx.az), interval=c("confidence"), 
              level=0.95, type="response")

#back-calculate prediction vector from log
prd.az <- 10^prd.az



#Montana
data.mt <- read.csv("/Users/atredenn/Documents/Projects/Diversity_Stability/Montana_allrecords_cover.csv")
data.mt <- as.data.frame(data.mt)

df.agg.mt <- ddply(data.mt, .(quad, year, Species), summarise, 
                sum = length(Species)
)

df.spp.mt <- ddply(df.agg.mt, .(quad, year), summarise,
                sum = sum(sum),
                spp = length(Species)
)

df.cv.mt <- ddply(df.spp.mt, .(quad), summarise,
               cv = mean(sum)/sd(sum),
               spp = mean(spp),
               mu = mean(sum)
)

df.cv.mt <- df.cv.mt[!is.na(df.cv.mt$cv),]
df.cv.mt <- df.cv.mt[!is.infinite(df.cv.mt$cv),]


y.mt <- df.cv.mt$cv
x.mt <- df.cv.mt$spp
model.mt = lm(y.mt~x.mt)
summary(model.mt)

#create prediction vector
newx.mt=seq(1,8, 0.01)
prd.mt <- predict(model.mt, newdata=data.frame(x.mt=newx.mt), interval=c("confidence"), 
               level=0.95, type="response")



#Kansas
data.ks <- read.csv("/Users/atredenn/Documents/Projects/Diversity_Stability/KansasData_Adler_Clean.csv")
data.ks <- as.data.frame(data.ks)


#Chunk of code for splitting plot and year in original dataset
# data.ks <- read.csv("/users/atredenn/desktop/KS_allrecords.csv")
# data.ks <- data.ks[data.ks$species!="Bare ground",]
# 
# years = numeric(length(data.ks$plotyear))
# plot = character(length(data.ks$plotyear))
# 
# for(i in 1:length(plot)){
#   years[i] <- as.numeric(substr(data.ks$plotyear[i], 
#                                 nchar(as.character(data.ks$plotyear[i]))-1,
#                                 nchar(as.character(data.ks$plotyear[i]))))
#   plot[i] <- substr(data.ks$plotyear[i], 
#                     1,
#                     nchar(as.character(data.ks$plotyear[i]))-2)
# }
# 
# years <- as.numeric(substr(data.ks$plotyear, length(data.ks$plotyear)-1,length(data.ks$plotyear)))
# plot <- substr(data.ks$plotyear, 1,(length(data.ks$plotyear)-2))
# 
# data.ks["years"] <- years
# data.ks["plot"] <- plot
# 
# write.csv(data.ks, file="/users/atredenn/desktop/KansasData_Adler_Clean.csv")

# data2 <- data2[data2$Species != "Unknown Forb",]

df.agg.ks <- ddply(data.ks, .(plot, years, species), summarise, 
                 sum = length(species)
)


df.spp.ks <- ddply(df.agg.ks, .(plot,years), summarise,
                sum = sum(sum),
                spp = length(species)
)

df.cv.ks <- ddply(df.spp.ks, .(plot), summarise,
                 cv = mean(sum)/sd(sum),
                 spp = mean(spp)
)

df.cv.ks <- df.cv.ks[!is.na(df.cv.ks$cv),]
df.cv.ks <- df.cv.ks[!is.infinite(df.cv.ks$cv),]



y.ks <- df.cv.ks$cv
x.ks <- df.cv.ks$spp
model.ks <- lm(y.ks~x.ks)
summary(model.ks)

#create prediction vector
newx.ks=seq(4,24)
prd.ks <- predict(model.ks, newdata=data.frame(x.ks=newx.ks), interval=c("confidence"), 
               level=0.95, type="response")


pdf("AdlerData_Stability.pdf", height=6, width=6)
par(mfrow=c(2,2))
plot(df.cv.az$spp, df.cv.az$cv, col="#00000050", pch=19, las=1, xlab="average number of species",
     ylab=expression(paste("abundance stability (", mu/sigma, ")")), cex=0.8, main="Arizona")
lines(newx.az, prd.az[,1], lwd=2)
lines(newx.az, prd.az[,2], lty="dashed", lwd=1.5)
lines(newx.az, prd.az[,3], lty="dashed", lwd=1.5)
text(3, 7, "p = 0.00018", cex=0.8)

plot(df.cv.mt$spp, df.cv.mt$cv, col="#00000050", pch=19, las=1, xlab="average number of species",
     ylab=expression(paste("abundance stability (", mu/sigma, ")")), cex=0.8, main="Montana")
lines(newx.mt, prd.mt[,1], lwd=2)
lines(newx.mt, prd.mt[,2], lty="dashed", lwd=1.5)
lines(newx.mt, prd.mt[,3], lty="dashed", lwd=1.5)
text(7, 1, "p = 0.00014", cex=0.8)

plot(df.cv.ks$spp, df.cv.ks$cv, col="#00000050", pch=19, las=1, xlab="average number of species",
     ylab=expression(paste("abundance stability (", mu/sigma, ")")), cex=0.8, main="Kansas")
lines(newx.ks, prd.ks[,1], lwd=2)
lines(newx.ks, prd.ks[,2], lty="dashed", lwd=1.5)
lines(newx.ks, prd.ks[,3], lty="dashed", lwd=1.5)
text(20, 0.8, "p < 0.0001", cex=0.8)


plot(df.cv.az$spp, df.cv.az$cv, xlim=c(1, 25), ylim=c(0.4,8), col="darkorange", cex=0.8,
     main="All Sites", xlab="log[average number of species]", 
     ylab=expression(paste("log[abundance stability (", mu/sigma, ")]")), log="xy")
points(df.cv.mt$spp, df.cv.mt$cv, col="darkblue",cex=0.8)
points(df.cv.ks$spp, df.cv.ks$cv, cex=0.8)
lines(newx.az, prd.az[,1], lwd=1, col="darkorange")
lines(newx.mt, prd.mt[,1], lwd=1, col="darkblue")
lines(newx.ks, prd.ks[,1], lwd=1)
legend(8, 8.5, legend=c("Arizona", "Montana", "Kansas"), 
       col=c("darkorange", "darkblue", "black"), pch=1, cex=0.8)

y <- c(log10(df.cv.az$cv), log10(df.cv.mt$cv), log10(df.cv.ks$cv))
x <- c(log10(df.cv.az$spp), log10(df.cv.mt$spp), log10(df.cv.ks$spp))
# x <- c(df.cv$spp, df.cv2$spp, df.cv3$spp)
lm.all <- lm(y~x)
summary(lm.all)
newx.all <- log10(seq(1,60))
prd.all <- predict(lm.all, newdata=data.frame(x=newx.all))

lines(10^newx.all, 10^prd.all, lwd=2)

dev.off()

