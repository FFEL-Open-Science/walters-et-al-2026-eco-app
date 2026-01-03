compute_rte <- function(par, ref_par, mu_ref_par, dA, ref_Mn, ref_A, ref_lam,
                        lmd_mr_par, lsd_sr_par, level, C, T) {
  
  if (!(level %in% c("vr", "ep"))) {
    stop("level must be either 'vr' or 'ep'")
  }
  
  if (level == "vr") {
    allowed_pars <- c("gamma", "phi1", "phi2", "phi3", "psi1", "psi2")  
  } else if(level == "ep") {
    allowed_pars <- c("gamma_alpha", "gamma_tair", "gamma_dsch", "gamma_skt2",
                      "gamma_est", "gamma_dd", "phi1_tair", "phi1_skt2", 
                      "phi1_skt1", "phi1_est", "phi1_dd", "phi2_tair", 
                      "phi2_skt2", "phi2_skt1", "phi2_est", "phi2_dd",
                      "phi3_tair", "phi3_skt2", "phi3_skt1", "phi3_est", 
                      "phi3_dd", "psi1_tair", "psi1_skt2", "psi1_skt1", 
                      "psi1_est", "psi2_tair", "psi2_skt2", "psi2_skt1", 
                      "psi2_est")  
  }
  

  if (!(par %in% allowed_pars)) {
    if(level == "vr") {
      stop("Check the value of par. It must be one of 'gamma', 'phi1', 'phi2',
         'phi3', 'psi1' or 'psi2'.")
    } else if (level == "ep") {
      stop("Check the value of par. It must be one of 'gamma_alpha', 
                      'gamma_tair', 'gamma_dsch', 'gamma_skt2',
                      'gamma_est', 'gamma_dd', 'phi1_tair', 'phi1_skt2', 
                      'phi1_skt1', 'phi1_est', 'phi1_dd', 'phi2_tair', 
                      'phi2_skt2', 'phi2_skt1', 'phi2_est', 'phi2_dd',
                      'phi3_tair', 'phi3_skt2', 'phi3_skt1', 'phi3_est', 
                      'phi3_dd', 'psi1_tair', 'psi1_skt2', 'psi1_skt1', 
                      'psi1_est', 'psi2_tair', 'psi2_skt2', 'psi2_skt1' or 
                      'psi2_est'.")
      
    }
  }
  
  # Extract list derivatives of A with respect of parameter provided as input to
  # the function
  dA_dpar <- purrr::map(dA, ~.[stringr::str_starts(names(.), paste0("d", par))])
  
  #  Identity matrix for computation of real-time elasticity of indirect effects
  #  of past changes in vital rates
  I <- diag(C)
  
  #  Vector of 1's for computation of real-time elasticity of indirect effects of
  #  past changes in vital rates
  e <- matrix(1, C, 1)
  
  # Empty objects to store computations
  eA_mu <- eA_sig <- eM_mu <- eM_sig <- array(NA, dim = c(C, C, T))
  
  tot_eA_mu <- tot_eA_sig <- tot_eM_mu <- tot_eM_sig <- rep(NA, T)
  
  P2 <- P3 <- array(0, dim = c(C ^ 2, C ^ 2, T))
  
  w_mu <- w_sig <- array(0, dim = c(C, C, C, T + 1)) 
  
  for (t in 1:T) {
    for (i in 1:C) {
      for (j in 1:C) {
        # Direct effects of changes in vital rates ----
        ## Compute real time elasticity in relation to the mean (mu) ----
        eA_mu[i, j, t] <- mu_ref_par * dA_dpar[[t]][[1]][i, j] * ref_Mn[t, j] / ref_lam[t]
        
        ## Compute real time elasticity in relation to the SD (sig) ----
        eA_sig[i, j, t] <- (ref_par[t] - mu_ref_par) * dA_dpar[[t]][[1]][i, j] * 
          ref_Mn[t, j] / ref_lam[t]
        
        # Indirect effects of past changes in vital rates ----
        
        ## Compute perturbation matrices ----
        row <- (C + 1) * i - C
        col <- (C + 1) * j - C
        
        ### Perturbation matrix for difference in mean ----
        P2[row, col, t] <- mu_ref_par * dA_dpar[[t]][[1]][i, j]
        
        ### Perturbation matrix for differences in standard deviation ----
        P3[row, col, t] <- (ref_par[t] - mu_ref_par) * dA_dpar[[t]][[1]][i, j]

        ## Compute real time elasticity in relation to the mean (mu) ----
        row <- (C * i - (C - 1)):(C * i)
        col <- (C * j - (C - 1)):(C * j)
        
        ## w_mu represents the perturbation to population structure caused by
        ## historical changes in the mean of a vital rate
        w_mu[, i, j, t + 1] <- (I - ref_Mn[t + 1, ] %*% t(e)) %*%
                               (P2[row, col, t] %*% ref_Mn[t, ] +
                               ref_A[[t]] %*% w_mu[, i, j, t]) / ref_lam[t]
        
        eM_mu[i, j, t] <- (t(e) %*% ref_A[[t]] %*% w_mu[, i, j, t]) / ref_lam[t]
      
        ## Compute real time elasticity in relation to the SD (sig) ----
        ## w_sig represents the perturbation to population structure caused by
        ## historical changes in the SD of a vital rate
        w_sig[, i, j, t + 1] <- (I - ref_Mn[t + 1, ] %*% t(e)) %*%
                                (P3[row, col, t] %*% ref_Mn[t, ] +
                                ref_A[[t]] %*% w_sig[, i, j, t]) / ref_lam[t]
        
        eM_sig[i, j, t] <- (t(e) %*% ref_A[[t]] %*% w_sig[, i, j, t]) / ref_lam[t]
        
      }
    }
    tot_eA_mu[t] <- sum(eA_mu[, , t])
    tot_eA_sig[t] <- sum(eA_sig[, , t])
    tot_eM_mu[t] <- sum(eM_mu[, , t])
    tot_eM_sig[t] <- sum(eM_sig[, , t])
  }
  
  # Compute temporal average of each real-time elasticity ----
  avg_eA_mu <- mean(tot_eA_mu)
  avg_eA_sig <- mean(tot_eA_sig)
  avg_eM_mu <- mean(tot_eM_mu)
  avg_eM_sig <- mean(tot_eM_sig)
  
  # Compute contributions to the difference log or ratio of geometric mean of lambda ----
  if (level == "vr") {
    contA_mu <- lmd_mr_par * avg_eA_mu
    contA_sig <- lsd_sr_par * avg_eA_sig
    contM_mu <- lmd_mr_par * avg_eM_mu
    contM_sig <- lsd_sr_par * avg_eM_sig
    tcont <- contA_mu + contA_sig + contM_mu + contM_sig
  }
  
  if (level == "ep") {
    contA_mu <- lmd_mr_par ^ avg_eA_mu
    contA_sig <- lsd_sr_par ^ avg_eA_sig
    contM_mu <- lmd_mr_par ^ avg_eM_mu
    contM_sig <- lsd_sr_par ^ avg_eM_sig
    tcont <- contA_mu * contM_mu + contA_sig * contM_sig
  }
  
  
  # Compute total contributions to the difference in geometric mean of lambda ----
  out <- tibble(par = par, contA_mu = contA_mu, contA_sig = contA_sig,
                contM_mu = contM_mu, contM_sig = contM_sig, 
                tcont = tcont, avg_eA_mu = avg_eA_mu, avg_eA_sig = avg_eA_sig,
                avg_eM_mu = avg_eM_mu, avg_eM_sig = avg_eM_sig)
  
  return(out)
  
}