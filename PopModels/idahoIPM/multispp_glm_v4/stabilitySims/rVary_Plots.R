####
#### PLOT RESULTS FROM VARYING-r SIMULATIONS
####

library(ggplot2)

cv <- read.csv("CV_rVary.csv")
sdr <- read.csv("rSD_rVary.csv")
pgrs <- read.csv("pgrs_rVary.csv")

plot(density(pgrs[,2], adjust=5), col=1, lwd=3)
for(i in 3:5){
  lines(density(pgrs[,i], adjust=5), col=i, lwd=3)
}

tmpI <- which(pgrs<0, arr.ind = TRUE)[,1]

plotD <- as.data.frame(cbind(cv[,2], sdr[,2], apply(pgrs[,2:5], MARGIN = 1, FUN = mean)))
colnames(plotD) <- c("cv", "sdr", "meanR")

ggplot(plotD[-tmpI,], aes(x=meanR, y=cv))+
  geom_point(shape=19, size=3)+
  geom_smooth(color="purple", size=1, fill="purple", method="loess")+
  theme_bw()+
  xlab("sd(r)")+
  ylab(expression(CV[T]))

