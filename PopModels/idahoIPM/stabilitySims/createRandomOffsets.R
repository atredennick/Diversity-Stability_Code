sims <- 200
spp <- 4
offs <- matrix(NA, ncol=spp, nrow=sims)
for(i in 1:sims){
  ttt <- runif(4,0.8,1.2)
  offs[i,] <- ttt
}

write.csv(offs, "randomOffsets_rNEW.csv", row.names=FALSE)
