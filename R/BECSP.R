dBECSP_dt <- function(t, y, parms, P_str='uptake'){
  with(as.list(c(y, parms)), {
    uptake <- vmax_up * S * B / (km_up + S)
    if(P_str %in% 'uptake'){
      eProd <- enz_prod*uptake
    }else{
      eProd <- enz_prod*B
    }
    enzKin <- vmax_enz * C * E / (km_enz + C) #MM
    absorb <- (1-(P_S+P_C)/P_sat)/tau_a #Rs+P/tau_d
    #{\tau_a}(1-\frac{P_S+P_C}{P_{sat})})
    
    dB <- uptake*cue - eProd - eProd*enz_cost - B/tau_B
    dE <- eProd - E/tau_E
    dS <- inputC*(1-rho_c) + enzKin - uptake - S*absorb + P_S/tau_d
    dC <- inputC*rho_c - enzKin + E/tau_E + B/tau_B - C*absorb + P_C/tau_d
    dP_S <- S*absorb - P_S/tau_d
    dP_C <- C*absorb - P_C/tau_d

    list(c(dB, dE, dS, dC, dP_S, dP_C))#, dCO2))
  })
}

# $$B = I \tau_B (\epsilon_{cue}-p) $$
# $$E = I p \tau_E $$
# $$S = \frac{k}{\tau_B (\epsilon_{cue}-p) v - 1}$$
# $$C = \frac{k_m \rho_c}{p \tau_E v_m - \rho_c}$$
# $$P = (S+C) \frac{P_{sat}\tau_d}{(S+C)\tau_d + \tau_a P_{sat}}$$
# $$P_s = S\frac{\tau_d}{\tau_a}(1-\frac{P}{P_{sat}})$$
# $$P_c = C\frac{\tau_d}{\tau_a}(1-\frac{P}{P_{sat}})$$
SS_BECSP <- function(parms, P_str='uptake'){
  if(!grepl('^uptake$', P_str)){
    stop('Non-uptake production not coded for the steady state.')
  }
  #enz_prod, km_up, vmax_up, death, cue, inputC
  with(as.list(parms),{
    B <- inputC * tau_B * (cue-enz_prod) #Input and biomass param
    E <- inputC * tau_E * enz_prod # input and enzyme param
    S <- km_up / (vmax_up * tau_B * (cue-enz_prod) - 1) #biomass, enzyme prod
    C <- (km_enz * rho_c) / (p * tau_E * vmax_enz - rho_c) #enzyme, M-M conversion, C allocation
    #C+S with sorption kinetics
    P <- (S + C) * (P_sat * tau_d) / ( (S + C) * tau_d + tau_a * P_sat )
    P_s <- S * tau_d / tau_a * (1 - P / P_sat)
    P_c <- C * tau_d / tau_a * (1 - P / P_sat)
    
    list(B=B, E=E, S=S, C=C, P_s=P_s, P_c=P_c)
  })
}
