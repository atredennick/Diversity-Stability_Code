#Estimate sensitivity to competition for each species

#Sensitivity to competition: S = 1 - (r[invasion]-r[alone])

####
#### Community wide ----------------------------------------
####

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


####
#### Pair wise ----------------------------------------
####

invadePairs <- read.csv("invasionGrowthRates_YrFluct_PairWise.csv")
pairS <- matrix(NA, 4, 4)
for(i in 1:nrow(pairS)){
  pairS[i,] <- as.numeric(1 - (invadePairs[i,] / intrinsicD[i,1]))
}
pairS[which(pairS==1)] <- 0
pD <- as.data.frame(pairS)
# pD$species <- colnames(invadePairs)
colnames(pD) <- colnames(invadePairs)
rownames(pD) <- colnames(invadePairs)
library(xtable)
print(xtable(pD), type="html", file="pairwiseCompetitionEffects.html", digits=3)

