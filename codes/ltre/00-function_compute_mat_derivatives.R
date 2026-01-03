compute_mat_derivatives <- function(x, level) {
  # x: a named list of vital rates
  
  if(!(level %in% c("vr", "ep"))) {
    stop("Level must be either 'vr' or 'ep'.")
  }
  
  if (level == "vr") {
    # Symbolically define transition matrix ----
    A <- expression(phi1 * (1 - psi1), 0.5 * gamma, 0.5 * gamma,
                    
                    phi1 * psi1, phi2 * (1 - psi2), 0,
                    
                    0, phi2 * psi2, phi3)
      
    # Compute derivatives of A with respect to vital rates ----
    dgamma <- as.expression(sapply(A, D, "gamma"))
    dA_dgamma <- matrix(sapply(dgamma, eval, x), ncol = 3, nrow = 3, byrow = TRUE)
            
    dphi1 <- as.expression(sapply(A, D, "phi1"))
    dA_dphi1 <- matrix(sapply(dphi1, eval, x), ncol = 3, nrow = 3, byrow = TRUE)
    
    dphi2 <- as.expression(sapply(A, D, "phi2"))
    dA_dphi2 <- matrix(sapply(dphi2, eval, x), ncol = 3, nrow = 3, byrow = TRUE)
    
    dphi3 <- as.expression(sapply(A, D, "phi3"))
    dA_dphi3 <- matrix(sapply(dphi3, eval, x), ncol = 3, nrow = 3, byrow = TRUE)
    
    dpsi1 <- as.expression(sapply(A, D, "psi1"))
    dA_dpsi1 <- matrix(sapply(dpsi1, eval, x), ncol = 3, nrow = 3, byrow = TRUE)
    
    dpsi2 <- as.expression(sapply(A, D, "psi2"))
    dA_dpsi2 <- matrix(sapply(dpsi2, eval, x), ncol = 3, nrow = 3, byrow = TRUE)
    
    mat_deriv <- list(dgamma = dA_dgamma, dphi1 = dA_dphi1, dphi2 = dA_dphi2,
                     dphi3 = dA_dphi3, dpsi1 = dA_dpsi1, dpsi2 = dA_dpsi2)
  }
  
  if (level == "ep") {
    # Symbolically define transition matrix ----
    A <- expression(((1 / (1 + exp(-phi1_lpc))) / phi1_dd) *
                      (1 - ((1 / (1 + exp(-psi1_lpc))))),
                    
                    (0.5 * gamma_alpha * gamma_tair * gamma_dsch *
                      gamma_skt2 * gamma_est) / gamma_dd, 
                    
                    (0.5 * gamma_alpha * gamma_tair * gamma_dsch *
                      gamma_skt2 * gamma_est) / gamma_dd, 
                    
                    ((1 / (1 + exp(-phi1_lpc))) / phi1_dd) * 
                      ((1 / (1 + exp(-psi1_lpc)))),
                    
                    ((1 / (1 + exp(-phi2_lpc))) / phi2_dd) *
                      (1 - ((1 / (1 + exp(-psi2_lpc))))),
                    
                    0,
                    
                    0,
                    
                    ((1 / (1 + exp(-phi2_lpc))) / phi2_dd) *
                      ((1 / (1 + exp(-psi2_lpc)))),
                    
                    (1 / (1 + exp(-phi3_lpc))) / phi3_dd)
    
    # Compute derivatives of A with respect to predictors ----
    ## gamma ----
    dgamma_alpha <- as.expression(sapply(A, D, "gamma_alpha"))
    dA_dgamma_alpha <- matrix(sapply(dgamma_alpha, eval, x), 
                              ncol = 3, nrow = 3, byrow = TRUE)
    
    dgamma_tair <- as.expression(sapply(A, D, "gamma_tair"))
    dA_dgamma_tair <- matrix(sapply(dgamma_tair, eval, x), 
                              ncol = 3, nrow = 3, byrow = TRUE)
    
    dgamma_dsch <- as.expression(sapply(A, D, "gamma_dsch"))
    dA_dgamma_dsch <- matrix(sapply(dgamma_dsch, eval, x), 
                              ncol = 3, nrow = 3, byrow = TRUE)
    
    dgamma_skt2 <- as.expression(sapply(A, D, "gamma_skt2"))
    dA_dgamma_skt2 <- matrix(sapply(dgamma_skt2, eval, x), 
                              ncol = 3, nrow = 3, byrow = TRUE)
    
    dgamma_est <- as.expression(sapply(A, D, "gamma_est"))
    dA_dgamma_est <- matrix(sapply(dgamma_est, eval, x), 
                              ncol = 3, nrow = 3, byrow = TRUE)
    
    dgamma_dd <- as.expression(sapply(A, D, "gamma_dd"))
    dA_dgamma_dd <- matrix(sapply(dgamma_dd, eval, x), 
                              ncol = 3, nrow = 3, byrow = TRUE)
    
    ## phi1 ----
    dphi1_lpc <- as.expression(sapply(A, D, "phi1_lpc"))
    dA_dphi1_lpc <- matrix(sapply(dphi1_lpc, eval, x), 
                           ncol = 3, nrow = 3, byrow = TRUE)

    dphi1_dd <- as.expression(sapply(A, D, "phi1_dd"))
    dA_dphi1_dd <- matrix(sapply(dphi1_dd, eval, x), 
                           ncol = 3, nrow = 3, byrow = TRUE)
    
    ## phi2 ----
    dphi2_lpc <- as.expression(sapply(A, D, "phi2_lpc"))
    dA_dphi2_lpc <- matrix(sapply(dphi2_lpc, eval, x), 
                           ncol = 3, nrow = 3, byrow = TRUE)
    
    dphi2_dd <- as.expression(sapply(A, D, "phi2_dd"))
    dA_dphi2_dd <- matrix(sapply(dphi2_dd, eval, x), 
                          ncol = 3, nrow = 3, byrow = TRUE)
    
    ## phi3 ----
    dphi3_lpc <- as.expression(sapply(A, D, "phi3_lpc"))
    dA_dphi3_lpc <- matrix(sapply(dphi3_lpc, eval, x), 
                           ncol = 3, nrow = 3, byrow = TRUE)
    
    dphi3_dd <- as.expression(sapply(A, D, "phi3_dd"))
    dA_dphi3_dd <- matrix(sapply(dphi3_dd, eval, x), 
                          ncol = 3, nrow = 3, byrow = TRUE)
    
    ## psi1 ----
    dpsi1_lpc <- as.expression(sapply(A, D, "psi1_lpc"))
    dA_dpsi1_lpc <- matrix(sapply(dpsi1_lpc, eval, x), 
                           ncol = 3, nrow = 3, byrow = TRUE)
    
    ## psi2 ----
    dpsi2_lpc <- as.expression(sapply(A, D, "psi2_lpc"))
    dA_dpsi2_lpc <- matrix(sapply(dpsi2_lpc, eval, x), 
                           ncol = 3, nrow = 3, byrow = TRUE)
    
    mat_deriv <- list(dgamma_alpha = dA_dgamma_alpha, dgamma_tair = dA_dgamma_tair,
                     dgamma_dsch = dA_dgamma_dsch, dgamma_skt2 = dA_dgamma_skt2,
                     dgamma_est = dA_dgamma_est, dgamma_dd = dA_dgamma_dd,
                     
                     dphi1_tair = dA_dphi1_lpc, dphi1_skt2 = dA_dphi1_lpc,
                     dphi1_skt1 = dA_dphi1_lpc, dphi1_est = dA_dphi1_lpc,
                     dphi1_dd = dA_dphi1_dd,
                     
                     dphi2_tair = dA_dphi2_lpc, dphi2_skt2 = dA_dphi2_lpc,
                     dphi2_skt1 = dA_dphi2_lpc, dphi2_est = dA_dphi2_lpc,
                     dphi2_dd = dA_dphi2_dd,
                     
                     dphi3_tair = dA_dphi3_lpc, dphi3_skt2 = dA_dphi3_lpc,
                     dphi3_skt1 = dA_dphi3_lpc, dphi3_est = dA_dphi3_lpc,
                     dphi3_dd = dA_dphi3_dd,
                     
                     dpsi1_tair = dA_dpsi1_lpc, dpsi1_skt2 = dA_dpsi1_lpc,
                     dpsi1_skt1 = dA_dpsi1_lpc, dpsi1_est = dA_dpsi1_lpc,
                     
                     dpsi2_tair = dA_dpsi2_lpc, dpsi2_skt2 = dA_dpsi2_lpc,
                     dpsi2_skt1 = dA_dpsi2_lpc, dpsi2_est = dA_dpsi2_lpc)
  }
  
  
  return(mat_deriv)
}