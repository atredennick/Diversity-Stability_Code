####
#### Get random years for 2000 year simulation
####

allYrs=c(1:22)
yrVector=sample(allYrs,2000, replace = TRUE)
saveRDS(yrVector, "yrVector.rds")