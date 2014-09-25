#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))
library(ggplot2)

#bring in growth rates
pgrMax <- read.csv("Rmax_ipm.csv")
pgrMax

ggplot(data=pgrMax, aes(x=species, y=maxR))+
  geom_bar(stat="identity")+
  xlab("Species")+
  ylab("Intrinsic Growth Rate (r)")+
  theme_bw()
