extract_sample <- function(spl, dat, year.range, s) {
  
  # Get environmental data ----
  ztair <- dat$ztair[9:42]
  zskt1 <- dat$zskt1[9:42]
  zskt2 <- dat$zskt2[9:42]
  zdsch <- dat$zdsch[1:34]
  
  slist <- list()
  
  slist$M <- filter(spl$M_spl, iter == s, t <= length(year.range) - 1) |>
    pivot_wider(names_from = c, values_from = M) |>
    arrange(t) |>
    select(`1`:`3`)
  
  slist$Mn <- filter(spl$Mn_spl, iter == s, t <= length(year.range) - 1) |>
    arrange(t) |>
    select(`1`:`3`) 
  
  ### Get recruits to each size class ----
  slist$R <- filter(spl$R_spl, iter == s, t <= length(year.range) - 1) |>
    arrange(t) |>
    select(`1`:`3`) 
  
  ### Get survivors in each size class ----
  slist$S <- filter(spl$S_spl, iter == s, t <= length(year.range) - 1) |>
    arrange(t) |>
    select(`1`:`3`) 
  
  ### Get population growth rate (lambda) ----
  slist$lam <- filter(spl$lam_spl, iter == s, t <= length(year.range) - 1) |>
    arrange(t)
  
  ### Extract gamma and compute predictor components ----
  slist$gamma <- unlist(spl$gamma_spl[s, ])
  
  # We get a large approximation error when using the linear predictor gamma,  
  # as shown in the command commented out below. Gamma was identified as the 
  # parameter contributing to the most to the approximation error by comparing
  # it's variance to its approximate variance, as shown below.
  #
  # gamma_lpc <- cbind(spl$beta_gamma1_spl[s] * ztair,
  #                    spl$beta_gamma2_spl[s] * zdsch,
  #                    spl$beta_gamma3_spl[s] * zskt2,
  #                    log(spl$alpha_gamma[s]) + unlist(spl$eps_gamma_spl[s, ]))
  #
  # var_gamma <- var(exp(rowSums(gamma_lpc)))
  
  # app_var_gamma <- (exp(sum(colMeans(gamma_lpc))) ^ 2) * var(rowSums(gamma_lpc))
  #
  # So we use a different approach where we factor the components of the 
  # linear predictor and approximate them separately, as suggested by Jonas
  # Knape. This results in a much smaller approximation error.
  
  slist$gamma_alpha <- rep(spl$alpha_gamma[s], length(year.range) - 1)               
  
  slist$gamma_ztair <- exp(spl$beta_gamma1_spl[s] * ztair)
  
  slist$gamma_zdsch <- exp(spl$beta_gamma2_spl[s] * zdsch)
  
  slist$gamma_zskt2 <- exp(spl$beta_gamma3_spl[s] * zskt2)
  
  slist$gamma_eps <- exp(unlist(spl$eps_gamma_spl[s, ]))
  
  slist$gamma_dd <- unlist(1 + 0.5 * spl$beta_gamma_DD_spl[s] * (slist$M[, 2] + slist$M[, 3]))
  
  ### Extract phi and compute predictor components ----
  slist$phi1 <- unlist(spl$phi1_spl[s, ])
  slist$phi2 <- unlist(spl$phi2_spl[s, ])
  slist$phi3 <- unlist(spl$phi3_spl[s, ])
  
  slist$phi1_lpc <- cbind(spl$beta_phi11_spl[s] * ztair,
                          spl$beta_phi12_spl[s] * zskt2,
                          spl$beta_phi13_spl[s] * zskt1,
                          spl$beta0_phi1_spl[s] + unlist(spl$eps_phi1_spl[s, ]))
        
  slist$phi2_lpc <- cbind(spl$beta_phi21_spl[s] * ztair,
                          spl$beta_phi22_spl[s] * zskt2,
                          spl$beta_phi23_spl[s] * zskt1,
                          spl$beta0_phi2_spl[s] + unlist(spl$eps_phi2_spl[s, ]))
  
  slist$phi3_lpc <- cbind(spl$beta_phi31_spl[s] * ztair,
                          spl$beta_phi32_spl[s] * zskt2,
                          spl$beta_phi33_spl[s] * zskt1,
                          spl$beta0_phi3_spl[s] + unlist(spl$eps_phi3_spl[s, ]))
        
  slist$phi1_dd <- unlist(1 + spl$beta_phi_DD1_spl[s] * slist$M[, 1])
  
  slist$phi2_dd <- unlist(1 + spl$beta_phi_DD2_spl[s] * slist$M[, 2])
  
  slist$phi3_dd <- unlist(1 + spl$beta_phi_DD3_spl[s] * slist$M[, 3])
  
  ### Extract psi and compute predictor components ----
  slist$psi1 <- unlist(spl$psi1_spl[s, ])
  slist$psi2 <- unlist(spl$psi2_spl[s, ])
  
  slist$psi1_lpc <- cbind(spl$beta_psi1_spl[s] * ztair,
                          spl$beta_psi2_spl[s] * zskt2,
                          spl$beta_psi3_spl[s] * zskt1,
                          spl$beta0_psi1_spl[s] + unlist(spl$eps_psi_spl[s, ]))
  
  slist$psi2_lpc <- cbind(spl$beta_psi1_spl[s] * ztair,
                          spl$beta_psi2_spl[s] * zskt2,
                          spl$beta_psi3_spl[s] * zskt1,
                          spl$beta0_psi2_spl[s] + unlist(spl$eps_psi_spl[s, ]))
  
  return(slist)
}

