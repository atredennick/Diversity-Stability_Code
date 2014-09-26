#Estimate sensitivity to competition for each species

#Sensitivity to competition: S = 1 - (r[invasion]-r[alone])

invadeD <- read.csv("invasionGrowthRates_YrFluct.csv")
intrinsicD <- read.csv("intrinsicGrowthRates_YrFluct.csv")

S <- 1 - (invadeD$maxR/intrinsicD$maxR)
S <- as.data.frame(S)
S$species <- invadeD$species

library(ggplot2)
ggplot(S, aes(x=species, y=S))+
  geom_bar(stat="identity")+
  ylab("Sensitivity to Competition (S)")+
  theme_bw()
