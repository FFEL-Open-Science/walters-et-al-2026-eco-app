# This function computes the contribution of vital rates (implicitly including
# effects of environmental predictors and environmental stochasticity) and 
# population structure to variation in population growth rate (lambda).

ltre_var_vr <- function(spl, year.range, nS) {
  
  require(tidyverse)
  source("codes/ltre/00-function_get_samples.R")
  source("codes/ltre/00-function_extract_sample.R")
  source("codes/ltre/00-function_compute_sensitivities.R")
  
  # Get samples ----
  spl <- get_samples(spl, year.range)
  
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
    
    ## Combine all parameters in one matrix ----
    all_pars <- cbind(sl$gamma, sl$phi1, sl$phi2, sl$phi3, 
                      sl$psi1, sl$psi2, sl$Mn)
    
    ## Compute means for the time series ----
    
    mean_list <- list()

    mean_list$gamma_mean <- mean(sl$gamma)
    
    mean_list$phi1_mean <- mean(sl$phi1)
    
    mean_list$phi2_mean <- mean(sl$phi2)
    
    mean_list$phi3_mean <- mean(sl$phi3)
    
    mean_list$psi1_mean <- mean(sl$psi1)
    
    mean_list$psi2_mean <- mean(sl$psi2)
    
    mean_list$Mn1_mean <- mean(sl$Mn$`1`)
  
    mean_list$Mn2_mean <- mean(sl$Mn$`2`)
    
    mean_list$Mn3_mean <- mean(sl$Mn$`3`)
  
    sensitivities <- compute_sensitivities(x = mean_list, ltre = "var",
                                           level = "vr")
    
    # Compute temporal variance-covariance matrix among parameters
    all_pars_vcv <- cov(all_pars)
    
    D <- length(sensitivities)
    
    contmatrix <- matrix(NA, D, D)
    
    for (i in 1:D){
      for(j in 1:D){
        contmatrix[i, j] <- all_pars_vcv[i, j] * sensitivities[i] * sensitivities[j]
      }
    }
    
    contributions_es <- rowSums(contmatrix)
    
    par_names <- c("gamma", "phi1", "phi2", "phi3", "psi1", "psi2",
                   "M1", "M2", "M3")
    
    names(contributions_es) <- par_names
    
    type <- c(rep('Vital rate', 6), rep('Population structure', 3))
    
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
      lam_ds[t] <- lam_tt[t] - lam_es[t]
    }
    
    # Variance of overall growth rate (includes both environmental and demographic
    # stochasticity)
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
                                type = type, 
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
                          rel_cont_lam_es = rel_cont_lam_es,
                          rel_cont_lam_ds = rel_cont_lam_ds, 
                          taylor_es = taylor_es, 
                          taylor_ds = taylor_ds,
                          taylor_tt = taylor_tt)
 
  return(output_ltre_var)

}
