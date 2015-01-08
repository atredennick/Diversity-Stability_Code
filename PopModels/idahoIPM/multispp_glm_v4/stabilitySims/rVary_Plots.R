####
#### PLOT RESULTS FROM VARYING-r SIMULATIONS
####
rm(list=ls()) 
library(ggplot2)
library(reshape2)
library(plyr)
library(gridExtra)

####
#### Read in processed RDS data (see script below)
####
monD <- readRDS("monocultureCutTS_rVary.rds")
mixD <- readRDS("mixtureCutTS_rVary.rds")

####
#### Calculate delta Y: net biodiversity effect
####
# total cover from mixture for species relative cover
totMix <- apply(X = mixD[,1:4], MARGIN = 1, FUN = sum)
totMixD <- data.frame(totCov = totMix,
                      simulation = mixD$simulation)
totMixAvg <- ddply(totMixD, .(simulation), summarise,
                   meanCov = mean(totCov))

# average species cover in mixture
mixM <- melt(mixD[,c(1,2,3,4,6)], id.vars = "simulation")
sppMixAvg <- ddply(mixM, .(variable, simulation), summarise,
                   avgCov = mean(value))
sppMixAvg$totCov <- rep(totMixAvg$meanCov, times = 4)
sppMixAvg$relCov <- with(sppMixAvg, avgCov/totCov)

# average cover in monoculture by species and simulation
monM <- melt(monD[,c(1,2,3,4,6)], id.vars = "simulation")
sppMonAvg <- ddply(mixM, .(variable, simulation), summarise,
                   monCov = mean(value))

# combine data frames and sum across species
deltaYD <- merge(sppMixAvg, sppMonAvg, by = c("variable", "simulation"))
deltaYD$expSpp <- with(deltaYD, relCov*monCov)
eCov <- ddply(deltaYD, .(simulation), summarise,
              expY = sum(expSpp))
deltaY <- totMixAvg$meanCov - eCov$expY

####
#### Calculate time series metrics
####
# temporal stability (community)
mixD$totC <- apply(mixD[,1:4], MARGIN = 1, sum)
tsD <- ddply(mixD[,5:7], .(simulation), summarise,
             ts = mean(totC)/sd(totC))

# community synchrony
sigC <- ddply(mixD[,5:7], .(simulation), summarise,
              value = var(totC))
sigS <- ddply(mixM, .(simulation, variable), summarise,
              value = sd(value))
sigSsum <- ddply(sigS, .(simulation), summarise,
                 value = (sum(value))^2)
psiC <- sigC$value / sigSsum$value

# environmental synchrony (from monocultures)
monD$totC <- apply(monD[,1:4], MARGIN = 1, sum)
sigE <- ddply(monD[,5:7], .(simulation), summarise,
              value = var(totC))
sigSe <- ddply(monM, .(simulation, variable), summarise,
              value = sd(value))
sigSesum <- ddply(sigSe, .(simulation), summarise,
                 value = (sum(value))^2)
psiE <- sigE$value / sigSesum$value

psiD <- resid(lm(psiC~psiE))

outD <- data.frame(psiC=psiC, psiE=psiE, psiD=psiD, deltaY=deltaY, TS=tsD$ts)

mod <- lm(log(TS)~deltaY+psiE+psiD, data=outD)
summary(mod)
# plot(mod)
modPE <- lm(log(TS)~psiE, data=outD)
predE <- as.data.frame(exp(predict(modPE, newdata = data.frame(psiE=seq(0.2, 1,0.05)))))
predE$psiE <- seq(0.2, 1,0.05)
colnames(predE)[1] <- "TS"

g1 <- ggplot(outD, aes(x=psiE, y=TS))+
  geom_point(shape=1, size=4)+
  geom_line(data=predE, aes(x=psiE, y=TS), color="purple", size=1)+
#   stat_smooth(method="lm", color="purple", fill="purple", size=1)+
#   scale_y_continuous(trans = "log")+
  theme_bw()
g2 <- ggplot(outD, aes(x=abs(psiD), y=TS))+
  geom_point(shape=1, size=4)+
  stat_smooth(method="lm", color="purple", fill="purple", size=1)+
#   scale_y_continuous(trans = "log")+
  theme_bw()
g3 <- ggplot(outD, aes(x=deltaY, y=TS))+
  geom_point(shape=1, size=4)+
  stat_smooth(method="lm", color="purple", fill="purple", size=1)+
  theme_bw()
g <- grid.arrange(g1,g2,g3,ncol=1)

#==========================================================#
#==========================================================#
#==========================================================#
#==========================================================#
#==========================================================#
# 
# ####
# #### Read in RDS data
# ####
# mixDfull <- as.data.frame(readRDS("output/outTimeSeries_rVary_mixture.rds"))
# monDfull <- as.data.frame(readRDS("output/outTimeSeries_rVary_monoculture.rds"))
# 
# ####
# #### Get relevant time steps (every 50) and output
# ####
# tVec <- seq(1, 1000, 50)
# nSims <- length(unique(mixDfull$simulation))
# monD <- data.frame(ARTR=NA,HECO=NA,POSE=NA,PSSP=NA,timeStep=NA,simulation=NA)
# mixD <- data.frame(ARTR=NA,HECO=NA,POSE=NA,PSSP=NA,timeStep=NA,simulation=NA)
# for(i in 1:nSims){
#   tmpD <- subset(monDfull, simulation==i)
#   tmpD2 <- tmpD[tVec, ]
#   monD <- rbind(monD, tmpD2)
# }
# for(i in 1:nSims){
#   tmpD <- subset(mixDfull, simulation==i)
#   tmpD2 <- tmpD[tVec, ]
#   mixD <- rbind(mixD, tmpD2)
# }
# 
# monD <- monD[2:nrow(monD),]
# mixD <- mixD[2:nrow(mixD),]
# saveRDS(monD, "monocultureCutTS_rVary.rds")
# saveRDS(mixD, "mixtureCutTS_rVary.rds")


#### PREVIOUS PLOTS ####
# cv <- read.csv("CV_rVary.csv")
# sdr <- read.csv("rSD_rVary.csv")
# pgrs <- read.csv("pgrs_rVary.csv")
# monoSynch <- read.csv("monoSynch.csv")
# commSynch <- read.csv("commSynchrony_rVary.csv")
# 
# plot(monoSynch[,2], commSynch[,2])
# mod <- lm(commSynch[,2]~monoSynch[,2])
# abline(mod)
# resmod <- resid(mod)
# 
# plot(resmod, log(1/cv[,2]), xlab=expression(phi[D]), ylab="ln(1/CV)")
# lnCV <- log(1/cv[,2])
# abline(lm(lnCV~resmod), lwd=3)
# plot(commSynch[,2], log(1/cv[,2]))
# plot(monoSynch[,2], log(1/cv[,2]))
# 
# dfCV <- data.frame(TS = lnCV,
#                    phiD = resmod)
# ggplot(data=dfCV, aes(x=phiD, y=TS))+
#   geom_point(shape=1, size=4)+
#   stat_smooth(method = "lm", color="purple", fill="purple", size=1)+
#   ylab("ln(TS)")+
#   xlab(expression(phi[D]))+
#   theme_bw()
# 
# plot(density(pgrs[,2], adjust=5), col=1, lwd=3)
# for(i in 3:5){
#   lines(density(pgrs[,i], adjust=5), col=i, lwd=3)
# }
# 
# tmpI <- which(pgrs<0, arr.ind = TRUE)[,1]
# 
# plotD <- as.data.frame(cbind(cv[,2], sdr[,2], apply(pgrs[,2:5], MARGIN = 1, FUN = mean)))
# colnames(plotD) <- c("cv", "sdr", "meanR")
# 
# ggplot(plotD[-tmpI,], aes(x=sdr, y=cv))+
#   geom_point(shape=19, size=3)+
#   geom_smooth(color="purple", size=1, fill="purple", method="loess")+
#   theme_bw()+
#   xlab("sd(r)")+
#   ylab(expression(CV[T]))
# 
