dBECS_dt <- function(t, y, parms){
  with(as.list(c(y, parms)), {
  uptake <- vmax_up * S * B / (km_up + S)
  eProd <- enz_prod*uptake
  enzKin <- vmax_enz * C * E / (km_enz + C) #MM
  
  dB <- uptake*cue - eProd - B/death
  dE <- eProd - E/E_turnover
  dS <- enzKin - uptake
  dC <- inputC - enzKin + E/E_turnover + B/death
  dCO2 <- uptake*(1-cue)
  #ans is the change in y
  list(c(dB, dE, dS, dC))#, dCO2))
  })
}

SS_BECS <- function(parms){
  with(as.list(parms),{
  S <- km_up / (vmax_up * death * (cue-enz_prod) - 1)
  B <- -inputC *( 1/death+(enz_prod-1)*vmax_up*S/(km_up+S))^-1
  list(S=S, B=B)
  })
}