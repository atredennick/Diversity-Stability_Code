#Estimate sensitivity to competition for each species

#Sensitivity to competition: S = 1 - (r[invasion]-r[alone])

invadeD <- read.csv("invasionGrowthRates_YrFluct.csv")
intrinsicD <- read.csv("intrinsicGrowthRates_YrFluct.csv")

S <- 1 - (invadeD$maxR/intrinsicD$maxR)
S <- as.data.frame(S)
S$species <- invadeD$species

library(ggplot2)
library(reshape2)
ggplot(S, aes(x=species, y=S))+
  geom_bar(stat="identity")+
  ylab("Sensitivity to Competition (S)")+
  theme_bw()


allD <- merge(invadeD, intrinsicD, by = "species")
colnames(allD) <- c("Species", "Invasion", "Intrinsic")
mD <- melt(allD)
ggplot(mD, aes(x=Species, y=value, fill=variable))+
  geom_bar(stat="identity", position="dodge", width=0.8, color="white")+
  scale_fill_manual(labels=c("Invasion growth rate", "Intrinsic growth rate"),
                    values=c("coral", "steelblue"), name="")+
  ylab("Growth rate")+
  theme_bw()



