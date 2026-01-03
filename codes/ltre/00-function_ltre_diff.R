# This function computes the contribution of environmental predictors  
# and environmental stochasticity affecting vital rates, as well as that of 
# population structure to annual differences in population growth rate (lambda).

ltre_diff <- function(spl, dat, year.range, nS) {

  require(tidyverse)
  require(zoo)
  source("codes/ltre/00-function_get_samples.R")
  source("codes/ltre/00-function_extract_sample.R")
  source("codes/ltre/00-function_compute_sensitivities.R")
  
  # Get samples ----
  spl <- get_samples(spl, year.range)
  
  # Get environmental data ----
  ztair <- dat$ztair[9:42]
  zskt1 <- dat$zskt1[9:42]
  zskt2 <- dat$zskt2[9:42]
  zdsch <- dat$zdsch[1:34]
  
  # Compute LTRE ----
  ltre_diff <- list()
  diff_lam_es <- diff_tcon_es <- taylor_es <- data.frame()
  
  print("// Computing transient LTRE for annual difference in lambda//")
  bp <- txtProgressBar(min = 0, max = nS, style = 3, char = "+", width = 50)
  
  for (s in 1:nS) {
    ## Extract samples for iteration s of MCMC ----
    sl <- extract_sample(spl, dat, year.range, s)
    
    # Compute differences between consecutive years
    gamma_diff <- diff(sl$gamma)
    gamma_alpha_diff <- diff(sl$gamma_alpha)
    gamma_ztair_diff <- diff(sl$gamma_ztair)
    gamma_zdsch_diff <- diff(sl$gamma_zdsch)
    gamma_zskt2_diff <- diff(sl$gamma_zskt2)
    gamma_est_diff <- diff(sl$gamma_eps)
    gamma_dd_diff <- diff(sl$gamma_dd)
    
    phi1_diff <- diff(sl$phi1)
    phi2_diff <- diff(sl$phi2)
    phi3_diff <- diff(sl$phi3)
    
    phi1_ztair_diff <- diff(sl$phi1_lpc[, 1])
    phi1_zskt2_diff <- diff(sl$phi1_lpc[, 2])
    phi1_zskt1_diff <- diff(sl$phi1_lpc[, 3])
    phi1_est_diff <- diff(sl$phi1_lpc[, 4])
    phi1_dd_diff <- diff(sl$phi1_dd)
    
    phi2_ztair_diff <- diff(sl$phi2_lpc[, 1])
    phi2_zskt2_diff <- diff(sl$phi2_lpc[, 2])
    phi2_zskt1_diff <- diff(sl$phi2_lpc[, 3])
    phi2_est_diff <- diff(sl$phi2_lpc[, 4])
    phi2_dd_diff <- diff(sl$phi2_dd)
    
    phi3_ztair_diff <- diff(sl$phi3_lpc[, 1])
    phi3_zskt2_diff <- diff(sl$phi3_lpc[, 2])
    phi3_zskt1_diff <- diff(sl$phi3_lpc[, 3])
    phi3_est_diff <- diff(sl$phi3_lpc[, 4])
    phi3_dd_diff <- diff(sl$phi3_dd)
    
    psi1_diff <- diff(sl$psi1)
    psi2_diff <- diff(sl$psi2)
    
    psi1_ztair_diff <- diff(sl$psi1_lpc[, 1])
    psi1_zskt2_diff <- diff(sl$psi1_lpc[, 2])
    psi1_zskt1_diff <- diff(sl$psi1_lpc[, 3])
    psi1_est_diff <- diff(sl$psi1_lpc[, 4])
    
    psi2_ztair_diff <- diff(sl$psi2_lpc[, 1])
    psi2_zskt2_diff <- diff(sl$psi2_lpc[, 2])
    psi2_zskt1_diff <- diff(sl$psi2_lpc[, 3])
    psi2_est_diff <- diff(sl$psi2_lpc[, 4])
    
    Mn1_diff <- diff(sl$Mn$`1`)
    Mn2_diff <- diff(sl$Mn$`2`)
    Mn3_diff <- diff(sl$Mn$`3`)
    
    ## Compute means for sequential years in the time series ----
    
    mean_list <- list()
    
    # Compute means of consecutive years
    mean_list$gamma_mean <- rollmean(sl$gamma, k = 2, align = "left")
    mean_list$gamma_alpha_mean <- rollmean(sl$gamma_alpha, k = 2, align = "left")
    mean_list$gamma_ztair_mean <- rollmean(sl$gamma_ztair, k = 2, align = "left")
    mean_list$gamma_zdsch_mean <- rollmean(sl$gamma_zdsch, k = 2, align = "left")
    mean_list$gamma_zskt2_mean <- rollmean(sl$gamma_zskt2, k = 2, align = "left")
    mean_list$gamma_est_mean <- rollmean(sl$gamma_eps, k = 2, align = "left")
    mean_list$gamma_dd_mean <- rollmean(sl$gamma_dd, k = 2, align = "left")
    
    mean_list$phi1_mean <- rollmean(sl$phi1, k = 2, align = "left")
    mean_list$phi2_mean <- rollmean(sl$phi2, k = 2, align = "left")
    mean_list$phi3_mean <- rollmean(sl$phi3, k = 2, align = "left")
    
    mean_list$phi1_lpc_mean <- rowSums(apply(sl$phi1_lpc, 2, rollmean, k = 2, 
                                             align = "left"))
    mean_list$phi1_dd_mean <- rollmean(sl$phi1_dd, k = 2, align = "left")
    
    mean_list$phi2_lpc_mean <- rowSums(apply(sl$phi2_lpc, 2, rollmean, k = 2, 
                                             align = "left"))
    mean_list$phi2_dd_mean <- rollmean(sl$phi2_dd, k = 2, align = "left")
    
    mean_list$phi3_lpc_mean <- rowSums(apply(sl$phi3_lpc, 2, rollmean, k = 2, 
                                             align = "left"))
    mean_list$phi3_dd_mean <- rollmean(sl$phi3_dd, k = 2, align = "left")
    
    mean_list$psi1_lpc_mean <- rowSums(apply(sl$psi1_lpc, 2, rollmean, k = 2, 
                                             align = "left"))
    mean_list$psi2_lpc_mean <- rowSums(apply(sl$psi2_lpc, 2, rollmean, k = 2, 
                                             align = "left"))

    mean_list$Mn1_mean <- rollmean(sl$Mn$`1`, k = 2, align = "left")
    mean_list$Mn2_mean <- rollmean(sl$Mn$`2`, k = 2, align = "left")
    mean_list$Mn3_mean <- rollmean(sl$Mn$`3`, k = 2, align = "left")
    
    ## Compute sensitivities evaluated at the mean values of the parameters ----
    sensitivities <- compute_sensitivities(x = mean_list, ltre = "diff", 
                                           level = "ep")
    
    # Compute contributions of each parameter to differences in lambda ----
    ltre_diff[[s]] <- tibble(gamma_tair = gamma_ztair_diff * sensitivities$gamma_tair,
                             gamma_dsch = gamma_zdsch_diff * sensitivities$gamma_dsch,
                             gamma_skt2 = gamma_zskt2_diff * sensitivities$gamma_skt2,
                             gamma_est = gamma_est_diff * sensitivities$gamma_est,
                             gamma_dd = gamma_dd_diff * sensitivities$gamma_dd,
                             
                             phi1_tair = phi1_ztair_diff * sensitivities$phi1_tair,
                             phi1_skt2 = phi1_zskt2_diff * sensitivities$phi1_skt2,
                             phi1_skt1 = phi1_zskt1_diff * sensitivities$phi1_skt1,
                             phi1_est = phi1_est_diff * sensitivities$phi1_est,
                             phi1_dd = phi1_dd_diff * sensitivities$phi1_dd,
                            
                             phi2_tair = phi2_ztair_diff * sensitivities$phi2_tair,
                             phi2_skt2 = phi2_zskt2_diff * sensitivities$phi2_skt2,
                             phi2_skt1 = phi2_zskt1_diff * sensitivities$phi2_skt1,
                             phi2_est = phi2_est_diff * sensitivities$phi2_est,
                             phi2_dd = phi2_dd_diff * sensitivities$phi2_dd,
                             
                             phi3_tair = phi3_ztair_diff * sensitivities$phi3_tair,
                             phi3_skt2 = phi3_zskt2_diff * sensitivities$phi3_skt2,
                             phi3_skt1 = phi3_zskt1_diff * sensitivities$phi3_skt1,
                             phi3_est = phi3_est_diff * sensitivities$phi3_est,
                             phi3_dd = phi3_dd_diff * sensitivities$phi3_dd,
                             
                             psi1_tair = psi1_ztair_diff * sensitivities$psi1_tair,
                             psi1_skt2 = psi1_zskt2_diff * sensitivities$psi1_skt2,
                             psi1_skt1 = psi1_zskt1_diff * sensitivities$psi1_skt1,
                             psi1_est = psi1_est_diff * sensitivities$psi1_est,
                             
                             psi2_tair = psi2_ztair_diff * sensitivities$psi2_tair,
                             psi2_skt2 = psi2_zskt2_diff * sensitivities$psi2_skt2,
                             psi2_skt1 = psi2_zskt1_diff * sensitivities$psi2_skt1,
                             psi2_est = psi2_est_diff * sensitivities$psi2_est,
                             
                             M1 = Mn1_diff * sensitivities$M1,
                             M2 = Mn2_diff * sensitivities$M2,
                             M3 = Mn3_diff * sensitivities$M3) |>
        mutate(year = 1987 + 1:n(),
               iter = s) |>
        pivot_longer(gamma_tair:M3, names_to = "par", values_to = "cont_es") |>
        mutate(rel_cont_es = abs(cont_es) / sum(abs(cont_es)), .by = year) |>
        mutate(type_1 = case_when(str_starts(par, "gamma") ~ "Productivity",
                                  str_starts(par, "phi1") ~ "Survival (size class 1)",
                                  str_starts(par, "phi2") ~ "Survival (size class 2)",
                                  str_starts(par, "phi3") ~ "Survival (size class 3)",
                                  str_starts(par, "psi1") ~ "Transition (size class 1 to 2)",
                                  str_starts(par, "psi2") ~ "Transition (size class 2 to 3)",
                                  str_starts(par, "M1") ~ "M1",
                                  str_starts(par, "M2") ~ "M2",
                                  str_starts(par, "M3") ~ "M3"),
               type_2 = ifelse(str_starts(par, "M"), "Population structure", "Vital rate"))
    
    A <- R_mat <- S_mat <- list()
    R_es <- S_es <- matrix(NA, nrow = length(year.range), ncol = 3)
    lam_es <- lam_ds <- rep(NA, length(year.range) - 1)
    lam_tt <- filter(spl$lam_spl, iter == s)$lam
    
    for (t in 1:(length(year.range) - 1)) {
      A[[t]] <- matrix(c(sl$phi1[t] * (1 - sl$psi1[t]), 0.5 * sl$gamma[t], 0.5 * sl$gamma[t],
                         
                         sl$phi1[t] * sl$psi1[t], sl$phi2[t] * (1 - sl$psi2[t]), 0, 
                         
                         0, sl$phi2[t] * sl$psi2[t], sl$phi3[t]), nrow = 3, byrow = TRUE)
      
      R_mat[[t]] <- matrix(c(0, 0.5 * sl$gamma[t], 0.5 * sl$gamma[t],
                             
                             sl$phi1[t] * sl$psi1[t], 0, 0, 
                             
                             0, sl$phi2[t] * sl$psi2[t], 0), nrow = 3, byrow = TRUE)
      
      S_mat[[t]] <- matrix(c(sl$phi1[t] * (1 - sl$psi1[t]), 0, 0,
                             
                             0, sl$phi2[t] * (1 - sl$psi2[t]), 0, 
                             
                             0, 0, sl$phi3[t]), nrow = 3, byrow = TRUE)
      
      R_es[t + 1, ] <- R_mat[[t]] %*% unlist(sl$M[t, ])
      
      S_es[t + 1, ] <- S_mat[[t]] %*% unlist(sl$M[t, ])
      
      # Multiply transition matrix by normalized population vector, so no need
      # to divide by total population size to get lambda.
      lam_es[t] <- sum((A[[t]] %*% unlist(sl$Mn[t, ])))
      
      # Compute demographic stochasticity by subtracting environmental stochasticity
      # growth rate from overall population growth rate
      lam_ds[t] <- filter(spl$lam_spl, iter == s)$lam[t] - lam_es[t]
    }
    
    diff_lam_es <-  rbind(diff_lam_es, 
                          data.frame(year = 1987 + 1:(length(year.range) - 2),
                                     diff = diff(lam_es),
                                     iter = s))
    diff_tcon_es <- rbind(diff_tcon_es, 
                          cbind(summarize(ltre_diff[[s]], 
                                          cont_es = sum(cont_es), 
                                          .by = year), 
                                iter = s))
    
    taylor_es <- rbind(taylor_es, 
                       data.frame(year = 1987 + 1:(length(year.range) - 2),
                                  taylor = diff_tcon_es / diff_lam_es,
                                  iter = s))
    
    setTxtProgressBar(bp, s)
    flush.console()
  }  
  
  colnames(taylor_es) <- c("year", "taylor_es", "iter")
  
  output_ltre_diff <- list(ltre_diff = ltre_diff,
                           diff_lam_es = diff_lam_es,
                           diff_tcon_es = diff_tcon_es,
                           taylor_es = taylor_es)
  
  close(bp)
  
  return(output_ltre_diff)
  
}