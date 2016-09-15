dBECS_dt <- function(t, y, parms){
  with(as.list(c(y, parms)), {
  uptake <- vmax_up * S * B / (km_up + S)
  eProd <- enz_prod*uptake
  enzKin <- vmax_enz * C * E / (km_enz + C) #MM
  
  dB <- uptake*cue - eProd - B/death
  dE <- eProd*cue_enz - E/E_turnover
  dS <- enzKin - uptake
  dC <- inputC - enzKin + E/E_turnover + B/death - loss*C
  dCO2 <- uptake*(1-cue) + eProd*(1-cue_enz)
  #ans is the change in y
  list(c(dB, dE, dS, dC))#, dCO2))
  })
}

SS_BECS <- function(parms){
  #enz_prod, km_up, vmax_up, death, cue, inputC
  with(as.list(parms),{
  S <- km_up/death/(cue*(1-enz_prod)*vmax_up-1/death)
  C <- km_enz/E_turnover/(vmax_enz*cue*cue_enz*enz_prod-1/E_turnover)
  B <- (inputC-loss*C)*(km_up+S)/(vmax_up*S)
  uptake <-  vmax_up * S * B / (km_up + S)
  E <- E_turnover*enz_prod*uptake*cue_enz
  list(S=S, B=B, E=E, C=C)
  })
}

max_Bs_BECS <- function(parms){
  with(as.list(parms),{
    S <- km_up / (vmax_up * death * (cue-enz_prod) - 1)
    a <- death*vmax_up
    enz_prodOpt <- (-a*km_up^2-a*km_up*S+a*km_up+S*km_up^2+S^2)/(S*(km_up+S))
    B_opt <- -inputC *( 1/death+(enz_prodOpt-1)*vmax_up*S/(km_up+S))^-1
    return(list(enz_prodOpt=enz_prodOpt, B_opt=B_opt))
  })
#   enz_prodResponse <- adply(seq(0, 1, by=0.001), c(1), function(xx){
#     parms$enz_prod <- xx
#     ans <- SS_BECS(parms)
#     return(data.frame(enz_prod=xx, B=ans$B, S=ans$S))
#   })
#   return(enz_prodResponse[which.max(enz_prodResponse$B), 'enz_prod'])
}