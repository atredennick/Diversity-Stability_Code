#Look at species correlation in intrinsic yearly growth rates

library(corrplot)
library(reshape2)

allD <- read.csv("RmaxYearly_ipm.csv")
M <- cor(allD[,1:4])
corrplot.mixed(M, lower="ellipse", upper="number", col=c("tomato3","skyblue4"),tl.col="grey35")
M

