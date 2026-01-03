get_samples <- function(spl_df, year.range) {
  
  require(tidyverse)
  
  # Extract MCMC samples for each parameter ----
  # 
  # No need to consider indicator variables, as the sampled value of a parameter
  # will be zero if the associated random variable is zero (i.e. the parameter is
  # not sampled by the RJMCMC algorithm in Nimble).
  
  spl <- list()
  
  ## gamma parameters ----

  spl$gamma_spl <- select(spl_df, starts_with("gamma["))[, -length(year.range)] 
  
  spl$alpha_gamma_spl <- unlist(select(spl_df, starts_with("alpha_gamma")))
  
  spl$beta_gamma1_spl <- unlist(select(spl_df, starts_with("beta_gamma[1")))
  
  spl$beta_gamma2_spl <- unlist(select(spl_df, starts_with("beta_gamma[2")))
  
  spl$beta_gamma3_spl <- unlist(select(spl_df, starts_with("beta_gamma[3")))
  
  spl$eps_gamma_spl <- select(spl_df, starts_with("eps_gamma["))[-length(year.range)]
  
  spl$beta_gamma_DD_spl <- unlist(select(spl_df, starts_with("beta_gamma_DD")))
  
  ## phi parameters ----
  
  spl$phi1_spl <- select(spl_df, starts_with("phi[1,"))[, -length(year.range)]
  spl$phi2_spl <- select(spl_df, starts_with("phi[2,"))[, -length(year.range)]
  spl$phi3_spl <- select(spl_df, starts_with("phi[3,"))[, -length(year.range)]
  
  spl$beta0_phi1_spl <- unlist(select(spl_df, starts_with("beta0_phi[1")))
  spl$beta0_phi2_spl <- unlist(select(spl_df, starts_with("beta0_phi[2")))
  spl$beta0_phi3_spl <- unlist(select(spl_df, starts_with("beta0_phi[3")))
  
  spl$beta_phi11_spl <- unlist(select(spl_df, starts_with("beta_phi[1, 1")))
  spl$beta_phi21_spl <- unlist(select(spl_df, starts_with("beta_phi[2, 1")))
  spl$beta_phi31_spl <- unlist(select(spl_df, starts_with("beta_phi[3, 1")))
  
  spl$beta_phi12_spl <- unlist(select(spl_df, starts_with("beta_phi[1, 2")))
  spl$beta_phi22_spl <- unlist(select(spl_df, starts_with("beta_phi[2, 2")))
  spl$beta_phi32_spl <- unlist(select(spl_df, starts_with("beta_phi[3, 2")))
  
  spl$beta_phi13_spl <- unlist(select(spl_df, starts_with("beta_phi[1, 3")))
  spl$beta_phi23_spl <- unlist(select(spl_df, starts_with("beta_phi[2, 3")))
  spl$beta_phi33_spl <- unlist(select(spl_df, starts_with("beta_phi[3, 3")))
  
  spl$eps_phi1_spl <- select(spl_df, starts_with("eps_phi[1"))[-length(year.range)]
  spl$eps_phi2_spl <- select(spl_df, starts_with("eps_phi[2"))[-length(year.range)]
  spl$eps_phi3_spl <- select(spl_df, starts_with("eps_phi[3"))[-length(year.range)]
  
  spl$beta_phi_DD1_spl <- unlist(select(spl_df, starts_with("beta_phi_DD[1")))
  spl$beta_phi_DD2_spl <- unlist(select(spl_df, starts_with("beta_phi_DD[2")))
  spl$beta_phi_DD3_spl <- unlist(select(spl_df, starts_with("beta_phi_DD[3")))
  
  ## psi parameters ----
  # Get only years of surveys for psi (exclude earlier years with age data only)
  # and remove last year (psi predicted to use in interpolation for estimating M)
  
  spl$psi1_spl <- select(spl_df, starts_with("psi[1,", ignore.case = FALSE))[, 8 + 1:(length(year.range) - 1)]
  spl$psi2_spl <- select(spl_df, starts_with("psi[2,", ignore.case = FALSE))[, 8 + 1:(length(year.range) - 1)]
  
  spl$beta0_psi1_spl <- unlist(select(spl_df, starts_with("beta0_psi[1")))
  spl$beta0_psi2_spl <- unlist(select(spl_df, starts_with("beta0_psi[2")))
  
  spl$beta_psi1_spl <- unlist(select(spl_df, starts_with("beta_psi[1")))
  spl$beta_psi2_spl <- unlist(select(spl_df, starts_with("beta_psi[2")))
  spl$beta_psi3_spl <- unlist(select(spl_df, starts_with("beta_psi[3")))
  
  spl$eps_psi_spl <- select(spl_df, starts_with("eps_psi["))[8 + 1:(length(year.range) - 1)]
  
  spl$M_spl <- select(spl_df, starts_with("M[")) |>
    mutate(iter = 1:n()) |>
    pivot_longer(starts_with("M["), names_to = "idx", values_to = "M") |>
    mutate(c = as.numeric(str_extract_all(idx, "[[:digit:]]+", simplify = TRUE)[, 1]),
           t = as.numeric(str_extract_all(idx, "[[:digit:]]+", simplify = TRUE)[, 2]),
           year = 1987 + t) |>
    filter(year < 2023) |> # 2023 is projected in the model
    select(iter, t, c, M) |>
    arrange(iter, t)
  
  spl$Mn_spl <- mutate(spl$M_spl, M = M / sum(M), .by = c(iter, t)) |>
    pivot_wider(names_from = c, values_from = M)
  
  spl$R_spl <- select(spl_df, starts_with("R[")) |>
    mutate(iter = 1:n()) |>
    pivot_longer(starts_with("R["), names_to = "idx", values_to = "R") |>
    mutate(c = as.numeric(str_extract_all(idx, "[[:digit:]]+", simplify = TRUE)[, 1]),
           t = as.numeric(str_extract_all(idx, "[[:digit:]]+", simplify = TRUE)[, 2]),
           year = 1987 + t) |>
    filter(year < 2023) |> # 2023 is projected in the model
    select(iter, t, c, R) |>
    pivot_wider(names_from = c, values_from = R) |>
    arrange(iter, t) 
  
  spl$S_spl <- select(spl_df, starts_with("S[")) |>
    mutate(iter = 1:n()) |>
    pivot_longer(starts_with("S["), names_to = "idx", values_to = "S") |>
    mutate(c = as.numeric(str_extract_all(idx, "[[:digit:]]+", simplify = TRUE)[, 1]),
           t = as.numeric(str_extract_all(idx, "[[:digit:]]+", simplify = TRUE)[, 2]),
           year = 1987 + t) |>
    filter(year < 2023) |> # 2023 is projected in the model
    select(iter, t, c, S) |>
    pivot_wider(names_from = c, values_from = S) |>
    arrange(iter, t) 
  
  spl$lam_spl <- select(spl_df, starts_with("lam_M[")) |>
    mutate(iter = 1:n()) |>
    pivot_longer(starts_with("lam_M["), names_to = "idx", values_to = "lam") |>
    mutate(t = as.numeric(str_extract_all(idx, "[[:digit:]]+", simplify = TRUE)[, 1]),
           year = 1987 + t) |>
    filter(year < 2022) |> # lambda from 2022-2023 is projected in the model (i.e. not observed)
    select(iter, t, lam) |>
    arrange(iter, t)
  
  return(spl)
}