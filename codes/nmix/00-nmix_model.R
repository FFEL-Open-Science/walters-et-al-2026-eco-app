nmix_model_code <- nimbleCode({

  # Stage-structured and dynamic N-mixture model ----

  ## Likelihood ----

  # Initial population size
  for (c in 1:sC) {
    M[c, 1] ~ dpois(nu * mu_M[c])
  }

  # Model projects population abundance to the start of the time step in year
  # 2023 (t = 36, so loops below end in T + 1), so that we can get an interpolated
  # estimate of the population when the surveys occur in the summer (including
  # summer of 2022).
  
  ### Productivity model (gamma) ----
  for (t in 2:(T + 1)) {
    
    if (dd == "gamma" | dd == "all") {
      gamma[t - 1] <- exp(z_gamma[1] * beta_gamma[1] * ztair[t + 8 - 1] +
                            z_gamma[2] * beta_gamma[2] * zdsch[t - 1] +
                            z_gamma[3] * beta_gamma[3] * zskt2[t + 8 - 1] +
                            eps_gamma[t - 1]) * alpha_gamma /
        (1 + beta_gamma_DD * (sum(M[2:3, t - 1]) / 2))
    } else {
      gamma[t - 1] <- exp(z_gamma[1] * beta_gamma[1] * ztair[t + 8 - 1] +
                            z_gamma[2] * beta_gamma[2] * zdsch[t - 1] +
                            z_gamma[3] * beta_gamma[3] * zskt2[t + 8 - 1] +
                            eps_gamma[t - 1]) * alpha_gamma
    }
    
    eps_gamma[t - 1] ~ dnorm(0, sig_gamma)
    
  } 
  
  ### Survival model (phi) ----
  for (t in 2:(T + 1)) {
    for (c in 1:sC) {
      
      if (dd == "phi" | dd == "all") {
        phi[c, t - 1] <- ilogit(beta0_phi[c] +
                                z_phi[c, 1] * beta_phi[c, 1] * ztair[t + 8 - 1] +
                                z_phi[c, 2] * beta_phi[c, 2] * zskt2[t + 8 - 1] +
                                z_phi[c, 3] * beta_phi[c, 3] * zskt1[t + 8 - 1] +
                                eps_phi[c, t - 1]) / (1 + beta_phi_DD[c] * M[c, t - 1])
      } else {
        phi[c, t - 1] <- ilogit(beta0_phi[c] +
                                  z_phi[c, 1] * beta_phi[c, 1] * ztair[t + 8 - 1] +
                                  z_phi[c, 2] * beta_phi[c, 2] * zskt2[t + 8 - 1] +
                                  z_phi[c, 3] * beta_phi[c, 3] * zskt1[t + 8 - 1] +
                                  eps_phi[c, t - 1])
      }

      eps_phi[c, t - 1] ~ dnorm(0, sig_phi[c])

    }
  }

  
  for (t in 2:(T + 1)) {
    ### Numbers recruiting to a size class from t-1 to t
    R[1, t] ~ dpois(gamma[t - 1] * (sum(M[2:3, t - 1]) / 2))
    
    for (c in 2:sC) {
      R[c, t] ~ dbin(phi[c - 1, t - 1] * psi[c - 1, t + 8 - 1], M[c - 1, t - 1])
    }
    
    ### Numbers surviving and remaining in a size class from t-1 to t
    for (c in 1:sC) {
      S[c, t] ~ dbin(phi[c, t - 1] * (1 - psi[c, t + 8 - 1]), M[c, t - 1])
    }
  
  }
  
  ### Total population size at t ----
  for (t in 2:(T + 1)) {
    for (c in 1:sC) {
      M[c, t] <- R[c, t] + S[c, t]
    }
  }

  ### Total population size at time of survey ----
  for (t in 2:(T + 1)) {
    for (c in 1:sC){
      M_int[c, t - 1] <- round((1 - w[t - 1]) * M[c, t - 1] + w[t - 1] * M[c, t])
    }
  }


  ### Number of fish subject to being surveyed ----
  for (t in 1:(T - 1)) {
    for (j in 1:J[t]) {

      # Probability of being included in the snorkel survey
      omega[j, t] <- exp(-beta_omega / d[j, t])

      for (c in 1:sC) {
        # sy is the survey year index
        N[j, t, c] ~ dbin(omega[j, t], M_int[c, sy[t]])
      }
    }
  }

  ### Observation model ----
  for (t in 1:(T - 1)) {
    for (c in 1:sC) {
      for (j in 1:J[t]) {
        if (obs == "binomial") {
          Y[j, t, c] ~ dbin(p[c, t], N[j, t, c])
        } else if (obs == "beta-binomial") {
          Y[j, t, c] ~ dBetaBinom_One(N = N[j, t, c],
                                      shape1 = alpha_B[c, t],
                                      shape2 = beta_B[c, t])
        } else if (obs == "poisson") {
          Y[j, t, c] ~ dpois(N[j, t, c] * p[c, t])
        }
      }
    }
  }

  ### Detection probability model ----
  
  for (t in 1:(T - 1)) {
    eps_p[t] ~ dnorm(0, sig_p)
    
    for (c in 1:sC) {
      if (obs == "binomial") {
        logit(p[c, t]) <- beta0_p[c] +
          z_p[1] * beta_p[1] * zvis[t] +
          z_p[2] * beta_p[2] * zdis[t] +
          eps_p[t]
      } else if (obs == "beta-binomial") {
        logit(p[c, t]) <- beta0_p[c] +
          z_p[1] * beta_p[1] * zvis[t] +
          z_p[2] * beta_p[2] * zdis[t] +
          eps_p[t]
        
        alpha_B[c, t] <- theta_p[c] * p[c, t]
        
        beta_B[c, t] <- theta_p[c] * (1 - p[c, t])
      } else if (obs == "poisson") {
        log(p[c, t]) <- beta0_p[c] +
          z_p[1] * beta_p[1] * zvis[t] +
          z_p[2] * beta_p[2] * zdis[t] +
          eps_p[t]
      } 
    }
  }

  ## Priors ----

  nu ~ dgamma(1, rate = 1)

  ### Survival model parameters ----

  for (c in 1:sC) {
    beta0_phi[c] ~ dnorm(0, sd = 2.5)
    sig_phi[c] ~ T(dnorm(0, sd = 2.5), 0, )
    for (b in 1:3) {
      zeta_phi[c, b] ~ dbeta(1, 1)
      z_phi[c, b] ~ dbern(zeta_phi[c, b])
      beta_phi[c, b] ~ dnorm(0, sd = 2.5)
    }
    if (dd == "phi" | dd == "all") {
      beta_phi_DD[c] ~ T(dnorm(0, sd = 0.1), 0, )  
    }
  }

  ### Productivity parameters ----

  alpha_gamma ~ T(dnorm(1, sd = 1), 1, )

  for (b in 1:3) {
    zeta_gamma[b] ~ dbeta(1, 1)
    z_gamma[b] ~ dbern(zeta_gamma[b])
    beta_gamma[b] ~ dnorm(0, sd = 2.5)
  }

  if (dd == "gamma" | dd == "all") {
    beta_gamma_DD ~ T(dnorm(0, sd = 0.1), 0, )
  }

  sig_gamma ~ T(dnorm(0, sd = 2.5), 0, )

  ### Observation model ----

  for (c in 1:sC) {
    beta0_p[c] ~ dnorm(0, sd = 2.5)
  }
  
  for (b in 1:2) {
    zeta_p[b] ~ dbeta(1, 1)
    z_p[b] ~ dbern(zeta_p[b])
    beta_p[b] ~ dnorm(0, sd = 2.5)
  }
  
  sig_p ~ T(dnorm(0, sd = 2.5), 0, )
  
  if (obs == "beta-binomial") {
    for (c in 1:sC) {
      theta_p[c] ~ dgamma(0.01, rate = 0.01)
    }
  } 

  ### Temporary emigration ----

  beta_omega ~ T(dnorm(0, sd = 2.5), 0, )
  
  ## Derived quantities ----

  ### Quantities for model check ----

      for (t in 1:(T - 1)) {
        for (c in 1:sC) {
          for (j in 1:J[t]) {
          if (obs == "binomial") {
              Y_exp[j, t, c] <- N[j, t, c] * p[c, t]

              Y_var[j, t, c] <- N[j, t, c] * p[c, t] * (1 - p[c, t])

              Y_rep[j, t, c] ~ dbin(p[c, t], N[j, t, c])
          } else if (obs == "beta-binomial") {
            Y_exp[j, t, c] <- (N[j, t, c] * alpha_B[c, t]) / (alpha_B[c, t] + beta_B[c, t])

            Y_var[j, t, c] <- (N[j, t, c] * alpha_B[c, t] * beta_B[c, t] *
                                 (alpha_B[c, t] + beta_B[c, t] + N[j, t, c])) /
              (((alpha_B[c, t] + beta_B[c, t]) ^ 2) * (alpha_B[c, t] + beta_B[c, t] + 1))

            Y_rep[j, t, c] ~ dBetaBinom_One(N = N[j, t, c],
                                            shape1 = alpha_B[c, t],
                                            shape2 = beta_B[c, t])
          } else if (obs == "poisson") {
            Y_exp[j, t, c] <- N[j, t, c] * p[c, t]

            Y_var[j, t, c] <- N[j, t, c] * p[c, t]

            Y_rep[j, t, c] ~ dpois(N[j, t, c] * p[c, t])
          }
       }
    }
  }
  
  ### Population growth rate and population structure ----
  
  for (t in 2:T) {
    lam_M[t - 1] <- sum(M[1:3, t]) / sum(M[1:3, t - 1])
  }
  
  for (t in 1:T) {
    for (c in 1:sC) {
      eta_M[c, t] <- M[c, t] / sum(M[1:3, t])
    }
  }
  
  # Size class transition model ----

  ## Likelihood ----
  
  ### Transition probability model (psi) ----
  for (t in 1:sT) {
    eps_psi[t] ~ dnorm(0, sd = sig_psi)
    psi[3, t] <- 0
    for (c in 1:2) {
      logit(psi[c, t]) <- beta0_psi[c] +
                          z_psi[1] * beta_psi[1] * ztair[t] +
                          z_psi[2] * beta_psi[2] * zskt2[t] +
                          z_psi[3] * beta_psi[3] * zskt1[t] +
                          eps_psi[t]
      }
  }

  ### Matrix for size class transitions ----
  for (i in 1:I) {
    for (a in (fA[i] + 1):lA[i]) {
      Psi[1, 1, i, a - 1] <- 1 - psi[1, t_mat[i, a - 1]]
      Psi[1, 2, i, a - 1] <- psi[1, t_mat[i, a - 1]]
      Psi[1, 3, i, a - 1] <- 0

      Psi[2, 1, i, a - 1] <- 0
      Psi[2, 2, i, a - 1] <- 1 - psi[2, t_mat[i, a - 1]]
      Psi[2, 3, i, a - 1] <- psi[2, t_mat[i, a - 1]]

      Psi[3, 1, i, a - 1] <- 0
      Psi[3, 2, i, a - 1] <- 0
      Psi[3, 3, i, a - 1] <- 1
    }
  }

  for (i in 1:I) {
    for (a in (fA[i] + 1):lA[i]) {
      C[i, a] ~ dcat(Psi[C[i, a - 1], 1:3, i, a - 1])
    }
  }

  ## Priors ----
  for(c in 1:2) {
    beta0_psi[c] ~ dnorm(0, sd = 2.5)
  }
  
  for (b in 1:3) {
    zeta_psi[b] ~ dbeta(1, 1)
    z_psi[b] ~ dbern(zeta_psi[b])
    beta_psi[b] ~ dnorm(0, sd = 2.5)
  }
  
  sig_psi ~ T(dnorm(0, sd = 2.5), 0, )
  
  ## Derived quantities ----

  ### Quantities for model check ----

  for (i in 1:I) {
    C_rep[i, fA[i]] <- C[i, fA[i]] # set initial class to the initial observed class
    for (a in (fA[i] + 1):lA[i]) {
      C_rep[i, a] ~ dcat(Psi[C_rep[i, a - 1], 1:3, i, a - 1])
      for (c in 1:3) {
        Pnum_tmp[i, a, c] <- (equals(C_rep[i, a], c) - Psi[C_rep[i, a - 1], c, i, a - 1]) ^ 2
        Pden_tmp[i, a, c] <- Psi[C_rep[i, a - 1], c, i, a - 1]
      }
      Pres[i, a] <- sum(Pnum_tmp[i, a, 1:3] / Pden_tmp[i, a, 1:3]) 
    }
  }
  
})
