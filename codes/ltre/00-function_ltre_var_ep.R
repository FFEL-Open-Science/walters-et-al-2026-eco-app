# This function computes the contribution of environmental predictors  
# and environmental stochasticity affecting vital rates, as well as that of 
# population structure to variation in population growth rate (lambda).

ltre_var_ep <- function(spl, dat, year.range, nS) {
  
  require(tidyverse)
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
  ltre_var <- list()
  tot_cont_es <- tot_cont_ds <- tot_cont_esds <- rep(NA, nS)
  rel_cont_ds <- data.frame()
  var_lam_tt <- var_lam_es <- var_lam_ds <- rep(NA, nS)
  taylor_tt <- taylor_es <- taylor_ds <- rep(NA, nS)
  rel_cont_lam_es <- rel_cont_lam_ds <- rep(NA, nS)
  
  print("// Computing transient LTRE for variance of lambda//")
  bp <- txtProgressBar(min = 0, max = nS, style = 3, char = "+", width = 50)
  
  for (s in 1:nS) {
    ## Extract samples for iteration s of MCMC ----
    sl <- extract_sample(spl, dat, year.range, s)
    
    # Combine all parameters in one matrix
    
    all_pars <- cbind(sl$gamma_alpha, sl$gamma_ztair, sl$gamma_zdsch, 
                      sl$gamma_zskt2, sl$gamma_eps, sl$gamma_dd, sl$phi1_lpc, 
                      sl$phi2_lpc, sl$phi3_lpc, sl$phi1_dd, sl$phi2_dd, 
                      sl$phi3_dd, sl$psi1_lpc, sl$psi2_lpc, sl$Mn)

    par_names <- c("gamma_alpha", "gamma_tair", "gamma_dsch", "gamma_skt2", 
                   "gamma_est", "gamma_dd", 
                   "phi1_tair", "phi1_skt2", "phi1_skt1", "phi1_est",
                   "phi2_tair", "phi2_skt2", "phi2_skt1", "phi2_est",
                   "phi3_tair", "phi3_skt2", "phi3_skt1", "phi3_est",
                   "phi1_dd", "phi2_dd", "phi3_dd",
                   "psi1_tair", "psi1_skt2", "psi1_skt1", "psi1_est",
                   "psi2_tair", "psi2_skt2", "psi2_skt1", "psi2_est",
                   "M1", "M2", "M3")
    
    names(all_pars) <- par_names
    
    # Compute temporal variance-covariance matrix among parameters
    all_pars_vcv <- cov(all_pars)
    
    ## Compute means for the time series ----
    
    mean_list <- list()
    
    ### --- gamma

    mean_list$gamma_mean <- mean(sl$gamma_alpha) * mean(sl$gamma_ztair) * 
      mean(sl$gamma_zdsch) * mean(sl$gamma_zskt2) * mean(sl$gamma_eps) / 
      mean(sl$gamma_dd)
    
    ### --- components of gamma
    
    mean_list$gamma_alpha_mean <- mean(sl$gamma_alpha)
    mean_list$gamma_ztair_mean <- mean(sl$gamma_ztair)
    mean_list$gamma_zdsch_mean <- mean(sl$gamma_zdsch)
    mean_list$gamma_zskt2_mean <- mean(sl$gamma_zskt2)
    mean_list$gamma_est_mean <- mean(sl$gamma_eps)
    
    ### --- gamma dd
    mean_list$gamma_dd_mean <- mean(sl$gamma_dd)
    
    ### --- phi
    mean_list$phi1_mean <- plogis(sum(colMeans(sl$phi1_lpc))) / mean(sl$phi1_dd)

    mean_list$phi2_mean <- plogis(sum(colMeans(sl$phi2_lpc))) / mean(sl$phi2_dd)

    mean_list$phi3_mean <- plogis(sum(colMeans(sl$phi3_lpc))) / mean(sl$phi3_dd)

    ### --- lpc of phi
    mean_list$phi1_lpc_mean <- sum(colMeans(sl$phi1_lpc))

    mean_list$phi2_lpc_mean <- sum(colMeans(sl$phi2_lpc))

    mean_list$phi3_lpc_mean <- sum(colMeans(sl$phi3_lpc))
    
    ### --- phi dd
    mean_list$phi1_dd_mean <- mean(sl$phi1_dd)
    
    mean_list$phi2_dd_mean <- mean(sl$phi2_dd)
    
    mean_list$phi3_dd_mean <- mean(sl$phi3_dd)
    
    ### --- psi
    mean_list$psi1_mean <- plogis(sum(colMeans(sl$psi1_lpc)))

    mean_list$psi2_mean <- plogis(sum(colMeans(sl$psi2_lpc)))
     
    ### --- lpc of psi
    mean_list$psi1_lpc_mean <- sum(colMeans(sl$psi1_lpc))

    mean_list$psi2_lpc_mean <- sum(colMeans(sl$psi2_lpc))
    
    ### --- normalized M
    mean_list$Mn1_mean <- mean(sl$Mn$`1`)
  
    mean_list$Mn2_mean <- mean(sl$Mn$`2`)
    
    mean_list$Mn3_mean <- mean(sl$Mn$`3`)
  
    ## Compute sensitivities evaluated at the mean values of the parameters ----
    sensitivities <- compute_sensitivities(x = mean_list, ltre = "var", 
                                           level = "ep")
    
    # Compute contributions of each parameter to variability in lambda ----
    
    D <- length(sensitivities)
    
    contmatrix <- matrix(NA, D, D)
    
    for (i in 1:D){
      for(j in 1:D){
        contmatrix[i, j] <- all_pars_vcv[i, j] * sensitivities[i] * sensitivities[j]
      }
    }
    
    contributions_es <- rowSums(contmatrix)
    names(contributions_es) <- par_names
    
    type_1 <- c(rep("Productivity", 6), rep("Survival (size class 1)", 4),
                rep("Survival (size class 2)", 4), rep("Survival (size class 3)", 4),
                "Survival (size class 1)", "Survival (size class 2)",
                "Survival (size class 3)", rep("Transition (size class 1 to 2)", 4),
                rep("Transition (size class 2 to 3)", 4), 
                "M1", "M2", "M3")
    
    type_2 <- c(rep("Vital rate", 29), rep("Population structure", 3))
    
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
    
    # Variance of overall growth rate (includes both environmental and demographic
    # stochasticities)
    var_lam_tt[s] <- var(lam_tt)
    
    # Variance of environmental stochasticity growth rate
    var_lam_es[s] <- var(lam_es)
    
    # Variance of demographic stochasticity growth rate
    var_lam_ds[s] <- var(lam_ds)
    
    # Relative contributions of lam_es and lam_ds to lam_tt
    rel_cont_lam_es[s] <- (var(lam_es) + cov(lam_es, lam_ds)) / var(lam_tt)
    rel_cont_lam_ds[s] <- 1 - rel_cont_lam_es[s]
    
    # Contribution of environmental predictors to environmental stochasticity growth rate
    tot_cont_es[s] <- sum(contributions_es)
    rel_cont_es <- abs(contributions_es) / sum(abs(contributions_es))
    
    contributions_esds <- contributions_es + sensitivities * apply(all_pars, 2, cov, y = lam_ds)
    tot_cont_esds[s] <- sum(contributions_esds)
    rel_cont_esds <- abs(contributions_esds) / sum(abs(contributions_esds))
    
    # Decomposition of demographic stochasticity growth rate 
    R_dev <- sl$R[2:(length(year.range) - 1), ] - R_es[2:(length(year.range) - 1), ]
    
    S_dev <- sl$S[2:(length(year.range) - 1), ] - S_es[2:(length(year.range) - 1), ]
    
    D <- cbind(R_dev / rowSums(sl$M[1:(length(year.range) - 2), ]),
               S_dev / rowSums(sl$M[1:(length(year.range) - 2), ]))
    
    colnames(D) <- c("R1", "R2", "R3", "S1", "S2", "S3")
    
    contributions_ds <- rowSums(cov(D))
    tot_cont_ds[s] <- sum(contributions_ds)
    rel_cont_ds <- rbind(rel_cont_ds, 
                         c(abs(contributions_ds) / sum(abs(contributions_ds)), s))
    
    # Approximation error based on Taylor statistics
    taylor_es[s] <- tot_cont_es[s] / var_lam_es[s]
    taylor_ds[s] <- tot_cont_ds[s] / var_lam_ds[s]
    taylor_tt[s] <- tot_cont_esds[s] / var_lam_tt[s]
    
    ltre_var[[s]] <- data.frame(par = names(contributions_es),
                                type_1 = type_1,
                                type_2 = type_2,
                                sens = sensitivities, 
                                cont_es = contributions_es, 
                                cont_esds = contributions_esds, 
                                rel_cont_es = rel_cont_es,
                                rel_cont_esds = rel_cont_esds,
                                iter = s)
    
    rownames(ltre_var[[s]]) <- NULL
    
    setTxtProgressBar(bp, s)
    flush.console()
  }  
  
  colnames(rel_cont_ds) <- c("R1", "R2", "R3", "S1", "S2", "S3", "iter")
  
  close(bp)

  output_ltre_var <- list(ltre_var = ltre_var,
                          rel_cont_ds = rel_cont_ds,
                          var_lam_tt = var_lam_tt, 
                          var_lam_es = var_lam_es, 
                          var_lam_ds = var_lam_ds,
                          tot_cont_esds = tot_cont_esds,
                          tot_cont_es = tot_cont_es,
                          tot_cont_ds = tot_cont_ds,
                          tot_cont_esds = tot_cont_esds,
                          rel_cont_lam_es = rel_cont_lam_es,
                          rel_cont_lam_ds = rel_cont_lam_ds, 
                          taylor_es = taylor_es, 
                          taylor_ds = taylor_ds,
                          taylor_tt = taylor_tt)
 
  return(output_ltre_var)

}
