rm(list=ls())

#set working directory
wd = "/Users/atredenn/Documents/Projects/Diversity_Stability/"
setwd(wd)
getwd()

model <- function(rm1, rm2, 
                  N1.start, N2.start, 
                  K1, K2, 
                  dvar1, dvar2, 
                  evar1, evar2, 
                  beat12, beta21, 
                  time){
  r1 = numeric(time)
  r2 = numeric(time)
  N1 = numeric(time)
  N2 = numeric(time)
  Ntot = numeric(time)
  
  N1[1] <- N1.start
  N2[1] <- N2.start
  Ntot[1] <- N1[1] + N2[1]
  
  for(t in 2:time){
    #   whitevar1 <- rnorm(1, mean=0)
    whitevar2 <- rnorm(1, mean=0)
    whitevar1 <- rnorm(1,mean=whitevar2)
    dnoise1 <- rnorm(1, mean=0)
    dnoise2 <- rnorm(1, mean=0)
    
    sum.non1 = (beta12 * N2[t-1])/K2
    sum.non2 = (beta21 * N1[t-1])/K1
    
    r1[t] = rm1 * (1-(N1[t-1]/K1) - sum.non1) + evar1 * whitevar1 + (dvar1 * dnoise1/sqrt(N1[t-1]))
    
    r2[t] = rm2 * (1-(N2[t-1]/K2) - sum.non2) + evar2 * whitevar2 + (dvar2 * dnoise2/sqrt(N2[t-1]))
    
    
    N1[t] = N1[t-1] + N1[t-1]*r1[t]
    N2[t] = N2[t-1] + N2[t-1]*r2[t]
    Ntot[t] = N1[t] + N2[t]
  }
  
  pops <- matrix(ncol=3, nrow=time)
  pops[,1] <- N1
  pops[,2] <- N2
  pops[,3] <- Ntot
  return(pops)
}


##Run null model with low environmental response and no competition

rm1 = 0.5
rm2 = 0.8
K1 = 1000
K2 = 1500
evar1 = 0.02
evar2 = 0.02
dvar1 = 1
dvar2 = 1
beta12 = 0
beta21 = 0
time = 600


N1.start <- 1000
N2.start <- 1500


model.null <- model(rm1, rm2, 
                    N1.start, N2.start, 
                    K1, K2, 
                    dvar1, dvar2, 
                    evar1, evar2, 
                    beat12, beta21, 
                    time)

N1 <- model.null[,1]
N2 <- model.null[,2]
Ntot <- model.null[,3]

cv.tot = sd(Ntot)/mean(Ntot)
cv.tot

# pdf("/users/atredenn/desktop/2species_comp.pdf", width=6, height=6)
# par(mfrow=c(2,2))

time.plot = seq(1,time,1)
plot(time.plot, Ntot, type="l", col="grey45", ylim=c(600,3000), ylab="Biomass", xlab="Time")
lines(time.plot, N1, col="lightblue")
lines(time.plot, N2, col="orange")
abline(h=mean(Ntot), col="black", lwd=2)
abline(h=(sd(Ntot)+mean(Ntot)), col="black", lty="dashed")
abline(h=(mean(Ntot)-sd(Ntot)), col="black", lty="dashed")

abline(h=mean(N1), col="cadetblue", lwd=2)
abline(h=(sd(N1)+mean(N1)), col="cadetblue", lty="dashed")
abline(h=(mean(N1)-sd(N1)), col="cadetblue", lty="dashed")

abline(h=mean(N2), col="darkorange", lwd=2)
abline(h=(sd(N2)+mean(N2)), col="darkorange", lty="dashed")
abline(h=(mean(N2)-sd(N2)), col="darkorange", lty="dashed")

text(100, 3000, paste("c.v. = ", round(cv.tot, 3)))



#Run model with high environmental response for species 2
rm1 = 0.5
rm2 = 0.8
K1 = 1000
K2 = 1500
evar1 = 0.02
evar2 = 0.1
dvar1 = 1
dvar2 = 1
beta12 = 0
beta21 = 0
time = 600


N1.start <- 1000
N2.start <- 1500


model.env <- model(rm1, rm2, 
                   N1.start, N2.start, 
                   K1, K2, 
                   dvar1, dvar2, 
                   evar1, evar2, 
                   beat12, beta21, 
                   time)

N1 <- model.env[,1]
N2 <- model.env[,2]
Ntot <- model.env[,3]

cv.tot = sd(Ntot)/mean(Ntot)
cv.tot

time.plot = seq(1,time,1)
plot(time.plot, Ntot, type="l", col="grey45", ylim=c(600,3100), ylab="Biomass", xlab="Time")
lines(time.plot, N1, col="lightblue")
lines(time.plot, N2, col="orange")
abline(h=mean(Ntot), col="black", lwd=2)
abline(h=(sd(Ntot)+mean(Ntot)), col="black", lty="dashed")
abline(h=(mean(Ntot)-sd(Ntot)), col="black", lty="dashed")

abline(h=mean(N1), col="cadetblue", lwd=2)
abline(h=(sd(N1)+mean(N1)), col="cadetblue", lty="dashed")
abline(h=(mean(N1)-sd(N1)), col="cadetblue", lty="dashed")

abline(h=mean(N2), col="darkorange", lwd=2)
abline(h=(sd(N2)+mean(N2)), col="darkorange", lty="dashed")
abline(h=(mean(N2)-sd(N2)), col="darkorange", lty="dashed")

text(100, 3100, paste("c.v. = ", round(cv.tot, 3)))


#Run model with high competition
rm1 = 0.5
rm2 = 0.8
K1 = 1000
K2 = 1500
evar1 = 0.02
evar2 = 0.02
dvar1 = 1
dvar2 = 1
beta12 = 0.7
beta21 = 0.9
time = 600


N1.start <- 800
N2.start <- 300


model.env <- model(rm1, rm2, 
                   N1.start, N2.start, 
                   K1, K2, 
                   dvar1, dvar2, 
                   evar1, evar2, 
                   beat12, beta21, 
                   time)

N1 <- model.env[,1]
N2 <- model.env[,2]
Ntot <- model.env[,3]

cv.tot = sd(Ntot)/mean(Ntot)
cv.tot

time.plot = seq(1,time,1)
plot(time.plot, Ntot, type="l", col="grey45", ylim=c(200,1500), ylab="Biomass", xlab="Time")
lines(time.plot, N1, col="lightblue")
lines(time.plot, N2, col="orange")
abline(h=mean(Ntot), col="black", lwd=2)
abline(h=(sd(Ntot)+mean(Ntot)), col="black", lty="dashed")
abline(h=(mean(Ntot)-sd(Ntot)), col="black", lty="dashed")

abline(h=mean(N1), col="cadetblue", lwd=2)
abline(h=(sd(N1)+mean(N1)), col="cadetblue", lty="dashed")
abline(h=(mean(N1)-sd(N1)), col="cadetblue", lty="dashed")

abline(h=mean(N2), col="darkorange", lwd=2)
abline(h=(sd(N2)+mean(N2)), col="darkorange", lty="dashed")
abline(h=(mean(N2)-sd(N2)), col="darkorange", lty="dashed")

text(100, 1500, paste("c.v. = ", round(cv.tot, 3)))



#Run model with high competition and environmental response
rm1 = 0.5
rm2 = 0.8
K1 = 1000
K2 = 1500
evar1 = 0.02
evar2 = 0.1
dvar1 = 1
dvar2 = 1
beta12 = 0.7
beta21 = 0.9
time = 600


N1.start <- 800
N2.start <- 300


model.env <- model(rm1, rm2, 
                   N1.start, N2.start, 
                   K1, K2, 
                   dvar1, dvar2, 
                   evar1, evar2, 
                   beat12, beta21, 
                   time)

N1 <- model.env[,1]
N2 <- model.env[,2]
Ntot <- model.env[,3]

cv.tot = sd(Ntot)/mean(Ntot)
cv.tot

time.plot = seq(1,time,1)
plot(time.plot, Ntot, type="l", col="grey45", ylim=c(200,1500), ylab="Biomass", xlab="Time")
lines(time.plot, N1, col="lightblue")
lines(time.plot, N2, col="orange")
abline(h=mean(Ntot), col="black", lwd=2)
abline(h=(sd(Ntot)+mean(Ntot)), col="black", lty="dashed")
abline(h=(mean(Ntot)-sd(Ntot)), col="black", lty="dashed")

abline(h=mean(N1), col="cadetblue", lwd=2)
abline(h=(sd(N1)+mean(N1)), col="cadetblue", lty="dashed")
abline(h=(mean(N1)-sd(N1)), col="cadetblue", lty="dashed")

abline(h=mean(N2), col="darkorange", lwd=2)
abline(h=(sd(N2)+mean(N2)), col="darkorange", lty="dashed")
abline(h=(mean(N2)-sd(N2)), col="darkorange", lty="dashed")

text(100, 1500, paste("c.v. = ", round(cv.tot, 3)))

# dev.off()



#Run model with dominant species
rm1 = 0.5
rm2 = 0.8
K1 = 1000
K2 = 1500
evar1 = 0.02
evar2 = 0.1
dvar1 = 1
dvar2 = 1
beta12 = 0.8
beta21 = 0.2
time = 600


N1.start <- 200
N2.start <- 1300


model.env <- model(rm1, rm2, 
                   N1.start, N2.start, 
                   K1, K2, 
                   dvar1, dvar2, 
                   evar1, evar2, 
                   beat12, beta21, 
                   time)

N1 <- model.env[,1]
N2 <- model.env[,2]
Ntot <- model.env[,3]

cv.tot = sd(Ntot)/mean(Ntot)
cv.tot

time.plot = seq(1,time,1)
plot(time.plot, Ntot, type="l", col="grey45", ylim=c(200,2500), ylab="Biomass", xlab="Time")
lines(time.plot, N1, col="steelblue")
lines(time.plot, N2, col="orange")
abline(h=mean(Ntot), col="black", lwd=2)
abline(h=(sd(Ntot)+mean(Ntot)), col="black", lty="dashed")
abline(h=(mean(Ntot)-sd(Ntot)), col="black", lty="dashed")

abline(h=mean(N1), col="steelblue", lwd=2)
abline(h=(sd(N1)+mean(N1)), col="steelblue", lty="dashed")
abline(h=(mean(N1)-sd(N1)), col="steelblue", lty="dashed")

abline(h=mean(N2), col="darkorange", lwd=2)
abline(h=(sd(N2)+mean(N2)), col="darkorange", lty="dashed")
abline(h=(mean(N2)-sd(N2)), col="darkorange", lty="dashed")

text(100, 2500, paste("c.v. = ", round(cv.tot, 3)))
