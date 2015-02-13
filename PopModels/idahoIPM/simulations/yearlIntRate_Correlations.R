#Look at species correlation in intrinsic yearly growth rates

library(corrplot)
library(reshape2)

allD <- read.csv("RmaxYearly_ipm.csv")
M <- cor(allD[,1:4])
corrplot.mixed(M, lower="ellipse", upper="number", col=c("tomato3","skyblue4"),tl.col="grey35")
M

cov(allD[,1:4])


matplot(allD$yearParams, as.matrix(allD[,1:4]), type="l")
matplot(allD$yearParams, as.matrix(allD[,1:4]), type="p", pch=1, add=TRUE)
synch <- sd(rowSums(allD[,1:4]))/(sd(allD[,1])+sd(allD[,2])+sd(allD[,3])+sd(allD[,4]))^2

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
library(ggplot2)
library(reshape2)
df <- melt(allD, id.vars = "yearParams")
ggplot(df, aes(x=yearParams, y=value, color=variable))+
  geom_line()+
  geom_point(size=4)+
  geom_hline(aes(yintercept=0), color="grey59")+
  xlab("Year")+
  ylab("Yearly intrinsic growth rate")+
  scale_color_manual(name="Species", values=cbPalette[1:4])+
  ggtitle(round(synch,2))
