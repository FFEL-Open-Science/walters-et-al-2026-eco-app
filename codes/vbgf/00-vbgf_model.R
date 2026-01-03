vbgf_model_code <- nimbleCode({

  # Von Bertalanffy growth function ----
  
  ## Likelihood ----
  
  for (i in 1:I) {
    
    mu_L[i, 1] <- L_inf * (1 - exp(-K[i] * (1 - a0)))
    
    L[i, 1] ~ dnorm(mu_L[i, 1], sd = sig_L)
    
    res_L[i, 1] <- L[i, 1] - mu_L[i, 1]
    
    for (a in 2:A[i]) {
      
      mu_L[i, a] <- L_inf * (1 - exp(-K[i] * (a - a0)))
      
      L[i, a] ~ dnorm(mu_L[i, a] + rho * res_L[i, a - 1], sd = sig_L)
      
      res_L[i, a] <- L[i, a] - mu_L[i, a]
      
    }
  }
  
  for (i in 1:I) {
    K[i] ~ T(dnorm(mu_K, sd = sig_K), 0, )
  }
  
  ## Priors ----
  sig_L ~ T(dnorm(0, sd = 5), 0, )
  
  L_inf ~ T(dnorm(52, sd = 5), 0, )
  
  a0 ~ T(dnorm(0, sd = 1), 0, )
  
  mu_K ~ T(dnorm(0, sd = 1), 0, )
  
  sig_K ~ T(dnorm(0, sd = 1), 0, )
  
  rho ~ dbeta(1, 1)
  
  ## Derived parameters ----
  
  ### Derived age at 10, 30 and 50 cm ----
  
  a10 <- a0 - log((L_inf - 10) / L_inf) / mu_K
  
  a30 <- a0 - log((L_inf - 30) / L_inf) / mu_K
  
  a50 <- a0 - log((L_inf - 50) / L_inf) / mu_K
  
  ### Quantities for model check ----
  
  for (i in 1:I) {
    
    L_rep[i, 1] ~ dnorm(mu_L[i, 1], sd = sig_L)
    
    res_L_rep[i, 1] <- L_rep[i, 1] - mu_L[i, 1]
    
    for (a in 2:A[i]) {
      
      L_rep[i, a] ~ dnorm(mu_L[i, a] + rho * res_L[i, a - 1], sd = sig_L)
      
      res_L_rep[i, a] <- L_rep[i, a] - mu_L[i, a]
      
    }
  }
})