rm(list=ls())

#set working directory
wd = "/Users/atredenn/Documents/Projects/Diversity_Stability/"
setwd(wd)
getwd()

hivar.file <- "/Users/atredenn/Documents/Projects/Diversity_Stability/Simulation_Output/stability_simulations_hivar.csv"
lovar.file <- "/Users/atredenn/Documents/Projects/Diversity_Stability/Simulation_Output/stability_simulations_lovar.csv"

cv.hi <- read.csv(hivar.file)
cv.lo <- read.csv(lovar.file)

cv.hi <- cv.hi[complete.cases(cv.hi),]
cv.lo <- cv.lo[complete.cases(cv.lo),]

hi.index <- sample(seq(1,length(cv.hi[,1])), 100)
lo.index <- sample(seq(1,length(cv.lo[,1])), 100)

cv.hi <- cv.hi[hi.index,]
cv.lo <- cv.lo[lo.index,]


N = matrix(ncol=5, nrow=length(cv.hi[,1]))
N[,1]=2
N[,2]=4
N[,3]=8
N[,4]=16
N[,5]=32
# N=log(N)

plot(N[,1], cv.lo[,2], xlim=c(min(N)-0.4,max(N)+0.4), pch=21, col="#00000050", cex=0.3,
     xlab="Number of Species", ylab=expression(paste("Stability (", mu/sigma, ")")),
     ylim = c(min(cv.hi[,2:6]), max(cv.lo[,2:6])), axes=FALSE)
axis(1, at=N[1,], labels=c("2", "4", "8", "16", "32"), cex.axis=0.8)
axis(2)
box()
points(N[,2], cv.lo[,3], pch=21, col="#00000050", cex=0.3)
points(N[,3], cv.lo[,4], pch=21, col="#00000050", cex=0.3)
points(N[,4], cv.lo[,5], pch=21, col="#00000050", cex=0.3)
points(N[,5], cv.lo[,6], pch=21, col="#00000050", cex=0.3)


N.sim=c(2,4,8,16,32)
cv.lo.long = c(cv.lo[,2], cv.lo[,3], cv.lo[,4], cv.lo[,5], cv.lo[,6])
spp = c(rep(N.sim[1], length(cv.lo[,1])),
        rep(N.sim[2], length(cv.lo[,1])),
        rep(N.sim[3], length(cv.lo[,1])),
        rep(N.sim[4], length(cv.lo[,1])),
        rep(N.sim[5], length(cv.lo[,1])))

mod <- lm(cv.lo.long ~ spp)
summary(mod)
pred = predict(mod)
pred = unique(pred)
lines(N[1,], pred)





plot(N[,1], cv.hi[,2], xlim=c(min(N)-0.4,max(N)+0.4), pch=21, col="#00000050", cex=0.3,
     xlab="Number of Species", ylab=expression(paste("Stability (", mu/sigma, ")")),
     ylim = c(min(cv.hi[,2:6]), max(cv.hi[,2:6])), axes=FALSE)
axis(1, at=N[1,], labels=c("2", "4", "8", "16", "32"), cex.axis=0.8)
axis(2)
box()
points(N[,2], cv.hi[,3], pch=21, col="#00000050", cex=0.3)
points(N[,3], cv.hi[,4], pch=21, col="#00000050", cex=0.3)
points(N[,4], cv.hi[,5], pch=21, col="#00000050", cex=0.3)
points(N[,5], cv.hi[,6], pch=21, col="#00000050", cex=0.3)

N.sim=c(2,4,8,16,32)
cv.hi.long = c(cv.hi[,2], cv.hi[,3], cv.hi[,4], cv.hi[,5], cv.hi[,6])
spp = c(rep(N.sim[1], length(cv.hi[,1])),
        rep(N.sim[2], length(cv.hi[,1])),
        rep(N.sim[3], length(cv.hi[,1])),
        rep(N.sim[4], length(cv.hi[,1])),
        rep(N.sim[5], length(cv.hi[,1])))

mod <- lm(cv.hi.long ~ spp)
summary(mod)
pred = predict(mod)
pred = unique(pred)
lines(N[1,], pred)






means = c(mean(cv.hi[,2]), mean(cv.hi[,3]), mean(cv.hi[,4]), mean(cv.hi[,5]), mean(cv.hi[,6]))
lines(N[1,], means, col="red")

points(N[1,1], mean(cv[,1]), pch=17, col="white", cex=2)
points(N[1,1], mean(cv[,1]), pch=17, col="red", cex=1)
points(N[1,2], mean(cv[,2]), pch=17, col="white", cex=2)
points(N[1,2], mean(cv[,2]), pch=17, col="red", cex=1)
points(N[1,3], mean(cv[,3]), pch=17, col="white", cex=2)
points(N[1,3], mean(cv[,3]), pch=17, col="red", cex=1)
points(N[1,4], mean(cv[,4]), pch=17, col="white", cex=2)
points(N[1,4], mean(cv[,4]), pch=17, col="red", cex=1)
points(N[1,5], mean(cv[,5]), pch=17, col="white", cex=2)
points(N[1,5], mean(cv[,5]), pch=17, col="red", cex=1)



