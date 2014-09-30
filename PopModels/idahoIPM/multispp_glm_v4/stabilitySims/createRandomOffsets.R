sims <- 100
spp <- 4
offs <- matrix(NA, ncol=spp, nrow=sims)
for(i in 1:100){
  ttt <- rnorm(4,1,0.5)
  offs[i,] <- ttt
}

write.csv(offs, "randomOffsets_r.csv", row.names=FALSE)
