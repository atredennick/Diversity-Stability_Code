env.stoch <- function(Nspecies, sigE, synch.targ){
  S <- Nspecies
  sigE <- sigE
  PhiE <- synch.targ
  
  if(PhiE!=1/S)
    betaz=(PhiE-1+sqrt(PhiE*(1-PhiE)*(S-1)))/(S*PhiE-1) else
      betaz=0.5
  
  EsigE=sigE
  AsigE=matrix(-1/S,nrow=S, ncol=S)
  diag(AsigE)=betaz-1/S
  AsigE=AsigE * EsigE/sqrt(betaz^2-2*betaz/S+1/S)
  
  #Finally, you get the vector of environmental stochasticity epsilon=sigE uei(t) with the following line of code:
  epsilon=AsigE %*% rnorm(S)
  return(epsilon)
}