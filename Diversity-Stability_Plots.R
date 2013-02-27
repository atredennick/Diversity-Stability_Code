rm(list=ls())

#set working directory
wd = "/Users/atredenn/Documents/Projects/Diversity_Stability/"
setwd(wd)
getwd()

hivar.file <- "/Users/atredenn/Documents/Projects/Diversity_Stability/Simulation_Output/stability_simulations_hivar.csv"
lovar.file <- "/Users/atredenn/Documents/Projects/Diversity_Stability/Simulation_Output/stability_simulations_lovar.csv"

# hivar.file <- "/Users/atredenn/desktop/stability_simulations_hivar.csv"
# lovar.file <- "/Users/atredenn/desktop/stability_simulations_lovar.csv"

cv.hi <- read.csv(hivar.file)
cv.lo <- read.csv(lovar.file)

cv.hi <- cv.hi[complete.cases(cv.hi),]
cv.lo <- cv.lo[complete.cases(cv.lo),]

# hi.index <- sample(seq(1,length(cv.hi[,1])), 100)
# lo.index <- sample(seq(1,length(cv.lo[,1])), 100)
# 
# cv.hi <- cv.hi[hi.index,]
# cv.lo <- cv.lo[lo.index,]


N = matrix(ncol=5, nrow=length(cv.lo[,1]))
N[,1]=2
N[,2]=4
N[,3]=8
N[,4]=16
N[,5]=32
# N=log(N)

# pdf(file="/users/atredenn/desktop/sims.pdf", height=8, width = 4)
par(mfrow=c(1,1))
par(mar=c(3.1,4.1,4.1,2.1), oma=c(0,0,0,0),mgp = c(1.8, 0.5, 0),tck=-0.02)
plot(N[,1], cv.lo[,2], xlim=c(min(N)-0.4,max(N)+0.4), pch=21, col="#00000050", cex=0.3,
     xlab="Number of Species", ylab=expression(paste("Stability (", mu/sigma, ")")),
     ylim = c(min(cv.lo[,2:6]), max(cv.lo[,2:6])), axes=FALSE, cex.lab=0.8)
axis(1, at=N[1,], labels=c("2", "4", "8", "16", "32"), cex.axis=0.8)
axis(2, cex.axis=0.8)
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

mod12 <- lm(cv.lo.long[1:length(cv.lo[,2])*2] ~ spp[1:length(cv.lo[,2])*2])
summary(mod12)

mod <- lm(cv.lo.long ~ spp)
summary(mod)
pred = predict(mod)
pred = unique(pred)
lines(N[1,], pred, lwd=2)

x = log10(spp)
# x=spp
mod2 <- lm(cv.lo.long ~ x)
summary(mod2)
newx=((seq(2,32, 0.1)))
pred = predict(mod2, newdata=data.frame(x=log10(newx)))
pred = unique(pred)
lines(newx, pred, col="darkorange", lwd=2)
eq1 = c(expression(paste(S[T]), "=", N[spp], "*1.1 + 65.2"))
legend(2,152,legend=c("",""), lwd=2, col=c("black", "darkorange"), bty="n", cex=0.7)




N = matrix(ncol=5, nrow=length(cv.hi[,1]))
N[,1]=2
N[,2]=4
N[,3]=8
N[,4]=16
N[,5]=32
plot(N[,1], cv.hi[,2], xlim=c(min(N)-0.4,max(N)+0.4), pch=21, col="#00000050", cex=0.3,
     xlab="Number of Species", ylab=expression(paste("Stability (", mu/sigma, ")")),
     ylim = c(min(cv.hi[,2:6]), max(cv.hi[,2:6]+5)), axes=FALSE, cex.lab=0.8)
axis(1, at=N[1,], labels=c("2", "4", "8", "16", "32"), cex.axis=0.8)
axis(2, cex.axis=0.8)
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


mod12 <- lm(cv.hi.long[1:length(cv.hi[,2])*2] ~ spp[1:length(cv.hi[,2])*2])
summary(mod12)

mod <- lm(cv.hi.long ~ spp)
summary(mod)
pred = predict(mod)
pred = unique(pred)
lines(N[1,], pred, lwd=2)



x=log10(spp)
mod1 <- lm(cv.hi.long ~ x)
summary(mod1)
newx=((seq(2,32, 0.1)))
pred = predict(mod1, newdata=data.frame(x=log10(newx)))
pred = unique(pred)
lines(newx, pred, col="darkorange", lwd=2)
legend(2,40,legend=c("",""), lwd=2, col=c("black", "darkorange"), bty="n", cex=0.7)


# dev.off()
AIC(mod, mod1)
# 
# 
# means = c(mean(cv.hi[,2]), mean(cv.hi[,3]), mean(cv.hi[,4]), mean(cv.hi[,5]), mean(cv.hi[,6]))
# lines(N[1,], means, col="red")
# 
# points(N[1,1], mean(cv[,1]), pch=17, col="white", cex=2)
# points(N[1,1], mean(cv[,1]), pch=17, col="red", cex=1)
# points(N[1,2], mean(cv[,2]), pch=17, col="white", cex=2)
# points(N[1,2], mean(cv[,2]), pch=17, col="red", cex=1)
# points(N[1,3], mean(cv[,3]), pch=17, col="white", cex=2)
# points(N[1,3], mean(cv[,3]), pch=17, col="red", cex=1)
# points(N[1,4], mean(cv[,4]), pch=17, col="white", cex=2)
# points(N[1,4], mean(cv[,4]), pch=17, col="red", cex=1)
# points(N[1,5], mean(cv[,5]), pch=17, col="white", cex=2)
# points(N[1,5], mean(cv[,5]), pch=17, col="red", cex=1)
# 
# 

