####
#### PLOT RESULTS FROM VARYING-r SIMULATIONS
####

library(ggplot2)

cv <- read.csv("CV_rVary.csv")
sdr <- read.csv("rSD_rVary.csv")

plotD <- as.data.frame(cbind(cv[,2], sdr[,2]))
colnames(plotD) <- c("cv", "sdr")

ggplot(plotD, aes(x=sdr, y=cv))+
  geom_point(shape=19, size=3)+
  geom_smooth(color="purple", size=1, fill="purple", method="loess")+
  theme_bw()+
  xlab("sd(r)")+
  ylab(expression(CV[T]))
