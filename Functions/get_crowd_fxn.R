get_crowding <- function(vital_data, dist_data, n_spp, do_spp, 
                         spp_list, alpha_effect){
  crowd <- matrix(NA,dim(vital_data)[1],n_spp)
  for(i in 1:dim(grow_data)[1]){
    tmp_data <- subset(dist_data,year==vital_data$year[i] & quad==vital_data$quad[i])
    focal <- which(tmp_data$trackID==vital_data$trackID[i] & tmp_data$nbSpp==do_spp)
    xx <- tmp_data$x[focal]
    yy <- tmp_data$y[focal]
    tmp_data$distance <- sqrt((xx-tmp_data$x)^2+(yy-tmp_data$y)^2)
    tmp_data=subset(tmp_data,distance>0)
    for(j in 1:n_spp){
      nbTmp <- which(tmp_data$nbSpp==spp_list[j])
      if(length(nbTmp)>0){
        crowd[i,j] <- sum(exp(-1*alpha_effect[j]*tmp_data$distance[nbTmp]^2)*tmp_data$area[nbTmp])
      }else{
        crowd[i,j]=0
      }
    }
  }
  return(crowd)
}