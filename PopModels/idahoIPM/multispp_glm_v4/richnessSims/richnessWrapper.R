####
#### Wrapper for species richness simulations
####
setwd("~/Repos/Diversity_Stability/PopModels/idahoIPM/multispp_glm_v4/richnessSims")
for(i in 1:(length(list.files())-2)){
  setwd("~/Repos/Diversity_Stability/PopModels/idahoIPM/multispp_glm_v4/richnessSims")
  source(list.files()[i])
}#next file