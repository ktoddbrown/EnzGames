dBE_dt <- function(t, y, parms){
  with(as.list(c(y, parms)), {
    uptake <- v * E / (k + E)
    prod <- p * cue * uptake
    
    dB <- (1-p)*cue*uptake - B/tau_B
    dE <- prod-E/tau_E
    list(c(dB, dE))
  })
}

BE_SS <- function(parms){
  with(as.list(c(parms)), {
  list(E = v*tau_E*p*cue-k, 
       B = tau_B/tau_E*(v*tau_E*cue*p-k)*(1/p-1))
  })
}