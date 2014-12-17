####
#### PLOT RESULTS FROM VARYING-r SIMULATIONS
####
rm(list=ls()) 
library(ggplot2)
library(reshape2)
library(plyr)
library(gridExtra)

####
#### Read in RDS data
####
mixDfull <- as.data.frame(readRDS("output/outTimeSeries_rVary_mixture.rds"))
monDfull <- as.data.frame(readRDS("output/outTimeSeries_rVary_monoculture.rds"))

####
#### Get relevant time steps (every 50)
####
tVec <- seq(1, 1000, 50)
nSims <- length(unique(mixDfull$simulation))
for(i in 1:nSims){
  tmpD <- subset(monDfull, simulation==i)
  tmpD2 <- tmpD[tVec, ]
}





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
