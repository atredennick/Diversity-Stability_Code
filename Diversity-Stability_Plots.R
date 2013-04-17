rm(list=ls())

#set working directory
wd = "/Users/atredenn/Documents/Projects/Diversity_Stability/"
setwd(wd)
getwd()

completes <- function(x){
  y <- x[complete.cases(x),]
  return(y)
}

hivar.file <- "/Users/atredenn/Documents/Projects/Diversity_Stability/Simulation_Output/stability_simulations_hivar.csv"
lovar.file <- "/Users/atredenn/Documents/Projects/Diversity_Stability/Simulation_Output/stability_simulations_lovar.csv"
nodom.file <- "/Users/atredenn/Documents/Projects/Diversity_Stability/Simulation_Output/stability_simulations_nodom.csv"
nodomhivar.file <- "/Users/atredenn/Documents/Projects/Diversity_Stability/Simulation_Output/stability_simulations_nodomhivar.csv"
domhivarall.file <- "/Users/atredenn/Documents/Projects/Diversity_Stability/Simulation_Output/stability_simulations_dom_hivar_ALL.csv"
domhivarrare.file <- "/Users/atredenn/Documents/Projects/Diversity_Stability/Simulation_Output/stability_simulations_dom_hivar_rare.csv"

cv.hi <- read.csv(hivar.file)
cv.lo <- read.csv(lovar.file)
cv.no <- read.csv(nodom.file)
cv.nohivar <- read.csv(nodomhivar.file)
cv.domhiall <- read.csv(domhivarall.file)
cv.domhirare <- read.csv(domhivarrare.file)

cv.hi <- completes(cv.hi)
cv.lo <- completes(cv.lo)
cv.no <- completes(cv.no)
cv.nohivar <- completes(cv.nohivar)
cv.domhiall <- completes(cv.domhiall)
cv.domhirare <- completes(cv.domhirare)

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
legend(2,164,legend=c("",""), lwd=2, col=c("black", "darkorange"), bty="n", cex=0.7)




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
legend(2,37,legend=c("",""), lwd=2, col=c("black", "darkorange"), bty="n", cex=0.7)



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

vectors <- function(x){
  y <- c(mean(x[,2]),
         mean(x[,3]),
         mean(x[,4]),
         mean(x[,5]),
         mean(x[,6])
         )
  return(y)
}


no.vec <- vectors(cv.no)
no.vec.s <- no.vec/no.vec
hi.vec.mean <- vectors(cv.hi)
hi.vec <- hi.vec.mean/no.vec
lo.vec.mean <- vectors(cv.lo)
lo.vec <- lo.vec.mean/no.vec
domhiall.vec.mean <- vectors(cv.domhiall)
domhiall.vec <- domhiall.vec.mean/no.vec
domhirare.vec.mean <- vectors(cv.domhirare)
domhirare.vec <- domhirare.vec.mean/no.vec
nohivar.vec.mean <- vectors(cv.nohivar)
nohivar.vec <- nohivar.vec.mean/no.vec

df <- data.frame(nodom.novar = no.vec.s,
                 nodom.hivar = nohivar.vec,
                 dom.hivar = hi.vec,
                 dom.lovar = lo.vec,
                 dom.hiall = domhiall.vec,
                 dom.hirare = domhirare.vec
                 )
library(reshape2)
df <- melt(df)
df$N <- rep(c(2,4,8,16,32),6)


library(ggplot2)
library(RColorBrewer)
colors <- brewer.pal(10, "Blues")
labels.sim <- c("NO DOM/LOW VAR",
                "NO DOM/HIGH VAR",
                "DOM/HIGH VAR DOM",
                "DOM/LOW VAR DOM",
                "DOM/HIGH VAR ALL",
                "DOM/HIGH VAR RARE")

ggplot(data = df) +
  geom_line(aes(x=N, y=value, color=variable)) +
  geom_point(aes(x=N, y=value), color="white", size=6) +
  geom_point(aes(x=N, y=value, color=variable, shape=variable), size=3) +
  theme_bw() +
  xlab("Number of Speices (N)") + 
  ylab(expression(paste("Relative Stability ( ", S[N]/S[1], ")"))) +
  scale_y_continuous(limits=c(0,1.2)) +
  scale_color_discrete(name = "Simulation", labels=labels.sim)+
  scale_shape_discrete(name = "Simulation", labels=labels.sim) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16, angle=90), 
        axis.text.x = element_text(size=14), 
        axis.text.y = element_text(size=14), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.key = element_blank(),
        legend.title = element_text(size=14)
  )


df.means <- data.frame(nodom.novar = no.vec,
                       nodom.hivar = nohivar.vec.mean,
                       dom.hivar = hi.vec.mean,
                       dom.lovar = lo.vec.mean,
                       dom.hiall = domhiall.vec.mean,
                       dom.hirare = domhirare.vec.mean
                       )

df.means <- melt(df.means)
df.means$N <- rep(c(1,2,3,4,5),6)

colors <- brewer.pal(8, "Greys")
ggplot(data=df.means) +
  geom_bar(aes(y=value, x=N, fill=variable), stat="identity", position="dodge") +
  theme_bw() +
  xlab("Number of Speices (N)") + 
  ylab(expression(paste("Average Temporal Stability (", mu/sigma, ")"))) +
  scale_fill_manual(name = "Simulation", 
                    labels=labels.sim,
                    values=rev(colors[2:8])) +
  scale_x_discrete(limits=c(1,2,3,4,5), labels=c(2,4,8,16,32)) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16, angle=90), 
        axis.text.x = element_text(size=14), 
        axis.text.y = element_text(size=14), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.key = element_blank(),
        legend.title = element_text(size=14)
  )



##CALCULATE DSR REGRESSION COEFFICIENTS

#set up functions
long.vector <- function(x){
  y <- c(x[,2], x[,3], x[,4], x[,5], x[,6])
  return(y)
}

N.regs <- function(x){
  N.sim <- c(2,4,8,16,32)
  spp <- c(rep(N.sim[1], length(x[,1])),
           rep(N.sim[2], length(x[,1])),
           rep(N.sim[3], length(x[,1])),
           rep(N.sim[4], length(x[,1])),
           rep(N.sim[5], length(x[,1]))
           )
  return(spp)
}

coef.get <- function(vec, spp){
  mod1 <- lm(vec ~ spp)
  aa <- mod1$coeff
  mod2 <- lm(vec ~ log(spp))
  bb <- mod2$coeff
  cc = matrix(ncol=3, nrow=2)
  cc[1,1:2] <- aa
  cc[2,1:2] <- bb
  aa.p <- summary.lm(mod1)$coefficients["spp","Pr(>|t|)"] 
  bb.p <- summary.lm(mod2)$coefficients["log(spp)","Pr(>|t|)"] 
  cc[1,3] <- aa.p
  cc[2,3] <- bb.p
  return(cc)
}

#No Dom / Low Var
vec1 <- long.vector(cv.no)
spp1 <- N.regs(cv.no)
coefs1 <- coef.get(vec1, spp1)

#No dom / high var
vec2 <- long.vector(cv.nohivar)
spp2 <- N.regs(cv.nohivar)
coefs2 <- coef.get(vec2, spp2)

#Dom / high dom var
vec3 <- long.vector(cv.hi)
spp3 <- N.regs(cv.hi)
coefs3 <- coef.get(vec3, spp3)

#Dom / lo dom var
vec4 <- long.vector(cv.lo)
spp4 <- N.regs(cv.lo)
coefs4 <- coef.get(vec4, spp4)

#Dom / high all var
vec5 <- long.vector(cv.domhiall)
spp5 <- N.regs(cv.domhiall)
coefs5 <- coef.get(vec5, spp5)

#Dom / high dom var
vec6 <- long.vector(cv.domhirare)
spp6 <- N.regs(cv.domhirare)
coefs6 <- coef.get(vec6, spp6)

coefs1
coefs2
coefs3
coefs4
coefs5
coefs6

st <- c(vec1, vec2, vec3, vec4, vec5, vec6)
N <- c(spp1, spp2, spp3, spp4, spp5, spp6)
sim <- c(rep(1,length(vec1)),
         rep(2,length(vec2)),
         rep(3,length(vec3)),
         rep(4,length(vec4)),
         rep(5,length(vec5)),
         rep(6,length(vec6)))

df.all <- data.frame(st <- st,
                     N <- N,
                     sim <- sim)

qplot(data=df.all, x=N....N, y=st....st, color=sim, geom=c("smooth", "point"), method='lm')


p <- ggplot(df.all, aes(N....N,st....st)) + geom_point(size=1.2)
p + facet_grid(. ~ sim....sim) + stat_smooth(method="lm", color="steelblue") + 
  stat_smooth(formula=y~log(x), method="lm", color="darkorange") +
  theme_bw() +
  xlab("Number of Speices (N)") + 
  ylab(expression(paste("Temporal Stability (", mu/sigma, ")"))) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16, angle=90), 
        axis.text.x = element_text(size=14), 
        axis.text.y = element_text(size=14), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.key = element_blank(),
        legend.title = element_text(size=14)
  )


