# This function computes the contribution of vital rates to the differences in
# geometric mean of population growth rates in two periods. The computations are
# done at the level of vital rates. Code was written based on codes included in
# supplementary material from Koons et al (2016, 2017).

ltre_lamg <- function(spl, dat, year.range, per1, per2, dir_diff, nS) {
  
  # Initial setup ----
  library(tidyverse)
  source("codes/ltre/00-function_get_samples.R")
  source("codes/ltre/00-function_extract_sample.R")
  source("codes/ltre/00-function_compute_mat_derivatives.R")
  source("codes/ltre/00-function_compute_log_diff.R")
  source("codes/ltre/00-function_compute_rte.R")
  
  per1 <- per1 - (min(year.range) - 1)
  per2 <- per2 - (min(year.range) - 1)

  if (length(per1) != length(per2)) {
    stop("The two periods must have the same number of years.")
  }
  
  T <- length(per1)
  C <- 3
  
  # Get samples ----
  spl <- get_samples(spl, year.range)

  # Compute LTRE ----
  lamg_per1 <- lamg_per2 <- rep(NA, nS)
  
  rte_gamma <- rte_phi1 <- rte_phi2 <- rte_phi3 <- tibble()
  rte_psi1 <- rte_psi2 <- diff_vr <- tibble()
  
  print("// Computing transient LTRE for differences in geometric mean of lambda//")
  bp <- txtProgressBar(min = 0, max = nS, style = 3, char = "+", width = 50)
  
  for (s in 1:nS) {
    ## Extract samples for iteration s of MCMC ----
    sl <- extract_sample(spl, dat, year.range, s)
    
    ## Compute geometric mean of lambda in the two periods ----
    lamg_per1[s] <- prod(filter(sl$lam, t %in% per1)$lam) ^ (1 / T)

    lamg_per2[s] <- prod(filter(sl$lam, t %in% per2)$lam) ^ (1 / T)
    
    ## Compute differences between the two periods ----
    ## Computations are based on a reference population that has: (1) an initial 
    ## structure that is defined by the mean of the first year structures for the 
    ## two periods; and (2) vital rates that are also defined by the mean vital
    ## rates of each time step for the two periods (Koons et al 2016).
    
    ### Define initial size for reference population ----
    ref_Mn <- matrix(NA, nrow = T + 1, ncol = 3)
    ref_Mn[1, ] <- c(mean(c(sl$Mn$`1`[1], sl$Mn$`1`[1 + T])),
                     mean(c(sl$Mn$`2`[1], sl$Mn$`2`[1 + T])),
                     mean(c(sl$Mn$`3`[1], sl$Mn$`3`[1 + T])))
    
    ### Compute dynamics of reference population ----
    ref_gamma <- ref_phi1 <- ref_phi2 <- ref_phi3 <- rep(NA, T)
    ref_psi1 <- ref_psi2 <- ref_lam <- rep(NA, T)
    ref_A <- ref_vr_list <- dA <- list()
    
    for (t in 1:T) {
      ### Per time-step mean of vital rates for the two periods
      ref_gamma[t] <- mean(c(sl$gamma[t + 0], sl$gamma[t + T]))
      ref_phi1[t] <- mean(c(sl$phi1[t + 0], sl$phi1[t + T]))
      ref_phi2[t] <- mean(c(sl$phi2[t + 0], sl$phi2[t + T]))
      ref_phi3[t] <- mean(c(sl$phi3[t + 0], sl$phi3[t + T]))
      ref_psi1[t] <- mean(c(sl$psi1[t + 0], sl$psi1[t + T]))
      ref_psi2[t] <- mean(c(sl$psi2[t + 0], sl$psi2[t + T]))
      
      ref_vr_list[[t]] <- list(gamma = ref_gamma[t], phi1 = ref_phi1[t], 
                               phi2 = ref_phi2[t], phi3 = ref_phi3[t], 
                               psi1 = ref_psi1[t], psi2 = ref_psi2[t])
      
      ref_A[[t]] <- matrix(c(ref_phi1[t] * (1 - ref_psi1[t]), 0.5 * ref_gamma[t], 0.5 * ref_gamma[t],
                             
                             ref_phi1[t] * ref_psi1[t], ref_phi2[t] * (1 - ref_psi2[t]), 0,
                             
                             0, ref_phi2[t] * ref_psi2[t], ref_phi3[t]), byrow = TRUE, nrow = 3)
      
      ### Compute lambda and project abundances for the reference population
      ref_lam[t] <- sum(ref_A[[t]] %*% ref_Mn[t, ])
      
      ref_Mn[t + 1, ] <- (ref_A[[t]] %*% ref_Mn[t, ]) / sum(ref_A[[t]] %*% ref_Mn[t, ])
      
      ### Compute derivatives of transition matrix A with respect to vital rates
      dA[[t]] <- compute_mat_derivatives(ref_vr_list[[t]], "vr")
    }
    
    ### Compute mean of reference vital rates ----
    mu_ref_gamma <- mean(ref_gamma)
    mu_ref_phi1 <- mean(ref_phi1)
    mu_ref_phi2 <- mean(ref_phi2)
    mu_ref_phi3 <- mean(ref_phi3)
    mu_ref_psi1 <- mean(ref_psi1)
    mu_ref_psi2 <- mean(ref_psi2)
    
    ### Compute log differences between vital rates for the two periods ----
    
    #### Difference between log mean ----
    lmd_gamma <- compute_log_diff(sl$gamma, "mu", dir_diff, T)
    lmd_phi1 <- compute_log_diff(sl$phi1, "mu", dir_diff, T)
    lmd_phi2 <- compute_log_diff(sl$phi2, "mu", dir_diff, T)
    lmd_phi3 <- compute_log_diff(sl$phi3, "mu", dir_diff, T)
    lmd_psi1 <- compute_log_diff(sl$psi1, "mu", dir_diff, T)
    lmd_psi2 <- compute_log_diff(sl$psi2, "mu", dir_diff, T)
    
    #### Difference between log sd ----
    lsd_gamma <- compute_log_diff(sl$gamma, "sig", dir_diff, T)
    lsd_phi1 <- compute_log_diff(sl$phi1, "sig", dir_diff, T)
    lsd_phi2 <- compute_log_diff(sl$phi2, "sig", dir_diff, T)
    lsd_phi3 <- compute_log_diff(sl$phi3, "sig", dir_diff, T)
    lsd_psi1 <- compute_log_diff(sl$psi1, "sig", dir_diff, T)
    lsd_psi2 <- compute_log_diff(sl$psi2, "sig", dir_diff, T)
    
    diff_vr <- rbind(diff_vr, 
                     tibble(par = c("gamma", "phi1", "phi2", "phi3", "psi1", "psi2"),
                            ln_mu_diff = c(lmd_gamma, lmd_phi1, lmd_phi2, lmd_phi3,
                                           lmd_psi1, lmd_psi2),
                            ln_sig_diff = c(lsd_gamma, lsd_phi1, lsd_phi2, lsd_phi3,
                                            lsd_psi1, lsd_psi2),
                            iter = s))
        
    ## Compute real-time elasticities ----
    
    rte_gamma <- rbind(rte_gamma, cbind(tibble(iter = s),
                       compute_rte("gamma", ref_gamma, mu_ref_gamma, dA, ref_Mn,
                                   ref_A, ref_lam, lmd_gamma, lsd_gamma, 
                                   level ="vr", C, T)))
    
    rte_phi1 <- rbind(rte_phi1, cbind(tibble(iter = s),
                       compute_rte("phi1", ref_phi1, mu_ref_phi1, dA, ref_Mn,
                                   ref_A, ref_lam, lmd_phi1, lsd_phi1, 
                                   level ="vr", C, T)))
    
    rte_phi2 <- rbind(rte_phi2, cbind(tibble(iter = s),
                      compute_rte("phi2", ref_phi2, mu_ref_phi2, dA, ref_Mn,
                                  ref_A, ref_lam, lmd_phi2, lsd_phi2, 
                                  level ="vr", C, T)))
    
    rte_phi3 <- rbind(rte_phi3, cbind(tibble(iter = s),
                      compute_rte("phi3", ref_phi3, mu_ref_phi3, dA, ref_Mn,
                                  ref_A, ref_lam, lmd_phi3, lsd_phi3, 
                                  level ="vr", C, T)))
    
    rte_psi1 <- rbind(rte_psi1, cbind(tibble(iter = s),
                      compute_rte("psi1", ref_psi1, mu_ref_psi1, dA, ref_Mn,
                                  ref_A, ref_lam, lmd_psi1, lsd_psi1, 
                                  level ="vr", C, T)))
    
    rte_psi2 <- rbind(rte_psi2, cbind(tibble(iter = s),
                      compute_rte("psi2", ref_psi2, mu_ref_psi2, dA, ref_Mn,
                                  ref_A, ref_lam, lmd_psi2, lsd_psi2, 
                                  level ="vr", C, T)))
    
    setTxtProgressBar(bp, s)
    flush.console()
  }
  
  close(bp)
  
  cont <- rbind(rte_gamma, rte_phi1, rte_phi2, rte_phi3, rte_psi1, rte_psi2)
  
  lamg <- tibble(iter = 1:nS, lamg_per1 = lamg_per1, lamg_per2 = lamg_per2)
  
  output_ltre_lamg <- list(lamg = lamg, cont = cont, diff_vr = diff_vr)
  
  return(output_ltre_lamg)

}
