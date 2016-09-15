dBEC_dt <- function(t, y, parms){
  with(as.list(c(y, parms)), {
    uptake <- (v * E / (k + E)) * C
    prod <- p * cue * uptake
    
    dB <- (1-p)*cue*uptake - B/tau_B
    dE <- prod-E/tau_E
    dC <- input + B/tau_B + E/tau_E - uptake
    list(c(dB, dE, dC))
  })
}