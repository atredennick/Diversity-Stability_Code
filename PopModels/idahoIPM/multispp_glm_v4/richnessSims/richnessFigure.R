####
#### Script to make prelim richness-stability figure
####

setwd("./output")
nFiles <- length(list.files())
out <- matrix(nrow=nFiles, ncol=2)
for(i in 1:nFiles){
  sing <- readRDS(list.files()[i])
  tmpComm <- apply(X = sing, MARGIN = 1, FUN = sum)
  tmpTS <- mean(tmpComm)/sd(tmpComm)
  tmpSpp <- apply(X = sing, MARGIN = 2, FUN = sum)
  sppR <- length(which(tmpSpp!=0))
  out[i,] <- c(tmpTS, sppR)
}

plot(out[,2], out[,1])

