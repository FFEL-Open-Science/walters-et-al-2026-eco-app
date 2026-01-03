compute_sensitivities <- function(x, ltre, level) {
  
  if(!(ltre %in% c("var", "diff"))) {
    stop("ltre must be either 'var' or 'diff'")
  }
  
  if(ltre == "var") {
    if(!(level %in% c("vr", "ep"))) {
      stop("Level must be either 'vr' or 'ep'.")
    }
  } else if (ltre == "diff") {
    if(level != "ep") {
      stop("Level must be 'ep'. Function hasn't been developed at vital rate 
           level when conducting LTRE for differences between consecutive years.")
    }
  }
  
  if (level == "vr") {
    # Sensitivities evaluated for vital rates (i.e. ignoring covariates) and 
    # population structure. The sensitivities are determined based on the following
    # expression for growth rate (lambda).
    
    # lam_ex <- expression((phi1 * (1 - psi1) * Mn1 + # survive and remain in 1
    #                      0.5 * gamma * (Mn2 + Mn3) + # recruit to 1 via reproduction
    #                      phi1 * psi1 * Mn1 + # survive and recruit to 2 via growth
    #                      phi2 * (1 - psi2) * Mn2 + # survive and remain in 2
    #                      phi2 * psi2 * Mn2 + # survive and recruit to 3 via growth
    #                      phi3 * Mn3) # survive and remain in 3)
    #                      / (Mn1 + Mn2 + Mn3)) # normalized total population size at previous time-step
    
    # Note: Although the derivatives will often have the total population
    # size at the previous time step in the denominator, we can omit them after 
    # differentiating the equation above because the normalized values sum to 1.
    
    # D(lam_ex, "gamma")
    # 0.5 * (Mn2_mean + Mn3_mean) / (Mn1_mean + Mn2_mean + Mn3_mean)
    gamma_sens <- 0.5 * (x$Mn2_mean + x$Mn3_mean)
    
    # D(lam_ex, "phi1")
    # Mn1_mean / (Mn1_mean + Mn2_mean + Mn3_mean)
    phi1_sens <- x$Mn1_mean
    
    # D(lam_ex, "phi2")
    # Mn2_mean / (Mn1_mean + Mn2_mean + Mn3_mean)
    phi2_sens <- x$Mn2_mean
    
    # D(lam_ex, "phi3")
    # Mn3_mean / (Mn1_mean + Mn2_mean + Mn3_mean)
    phi3_sens <- x$Mn3_mean
    
    # D(lam_ex, "psi1"), D(lam_ex, "psi2")
    psi1_sens <- psi2_sens <- 0 
    
    # The derivative with respect to M will have a term that is common to all 
    # three population sizes. We'll assign it to pd_M to avoid clutter below.
    pd_Mn <-  x$phi1_mean * x$Mn1_mean + 
              0.5 * x$gamma_mean * (x$Mn2_mean + x$Mn3_mean) + 
              x$phi2_mean * x$Mn2_mean + 
              x$phi3_mean * x$Mn3_mean
    
    # D(lam_exp, "Mn1")
    Mn1_sens <- x$phi1_mean - pd_Mn
    
    # D(lam_ex, "Mn2")
    Mn2_sens <- (0.5 * x$gamma_mean + x$phi2_mean) - pd_Mn
    
    # D(lam_ex, "Mn3")
    Mn3_sens <- (0.5 * x$gamma_mean + x$phi3_mean) - pd_Mn
    
    sensitivities <- c(gamma_sens, phi1_sens, phi2_sens, phi3_sens, psi1_sens, psi2_sens,
                       Mn1_sens, Mn2_sens, Mn3_sens)
  
    return(sensitivities)
  }
  
  if (level == "ep") {
    # Sensitivies evaluated for vital rates, predictors, density dependence and 
    # population structure. The sensitivities are determined using the following
    # expressions and the chain rule.
    
    lam_exp <- expression((phi1 * (1 - psi1) * Mn1 + # survive and remain in 1
                          0.5 * gamma * (Mn2 + Mn3) + # recruit to 1 via reproduction
                          phi1 * psi1 * Mn1 + # survive and recruit to 2 via growth
                          phi2 * (1 - psi2) * Mn2 + # survive and remain in 2
                          phi2 * psi2 * Mn2 + # survive and recruit to 3 via growth
                          phi3 * Mn3) # survive and remain in 3)
                          / (Mn1 + Mn2 + Mn3)) # normalized total population size at previous time-step
    
    # The expression below is simplified for computing derivatives. Folowing 
    # suggestion by Jonas Knape, we exponentiate each individual component of the
    # gamma equation and multiply them out(naming them Z1 for gamma, Z2 for ztair, 
    # Z3 for zdsch, Z4 for ztsk2, and Z5 for epsilon)
    gamma_exp <- expression(Z1 * Z2 * Z3 * Z4 * Z5 / dd)
    
    # The expression below is simplified for computing derivatives
    # expression((1 / (1 + exp(-lp))) / (1 + dd * M)) 
    phi_exp <- expression((1 / (1 + exp(-lp))) / dd) 
    
    # The expression below is simplified for computing derivatives
    # expression((1 / (1 + exp(-lp)))
    psi_exp <- expression(1 / (1 + exp(-lp)))
      
    # Note: Although the derivatives will often have the the total population
    # size at the previous time step in the denominator, we can omit them after 
    # differentiating the lam_exp equation above because the normalized values
    # sum to 1.

    # Derivative of lambda with respect to predictors of gamma using the chain 
    # rule. Symbolic computation below uses Z1, but the same applies to Z2, Z3,
    # and Z4.
    # D(lam_exp, "gamma") = 0.5 * (Mn2_mean + Mn3_mean)
    # 
    # D(gamma_exp, "Z1") =  Z2 * Z3 * Z4 / dd
    # 
    # D(lam_exp, "gamma") * D(gamma_exp, "Z1") = 0.5 *  (Mn2_mean + Mn3_mean) * Z2 * Z3 * Z4 / dd
    
    gamma_alpha_sens <- 0.5 * (x$Mn2_mean + x$Mn3_mean) * x$gamma_ztair_mean *
      x$gamma_zdsch_mean * x$gamma_zskt2_mean * x$gamma_est_mean / x$gamma_dd_mean
    
    gamma_ztair_sens <- 0.5 * (x$Mn2_mean + x$Mn3_mean) * x$gamma_alpha_mean * 
      x$gamma_zdsch_mean * x$gamma_zskt2_mean * x$gamma_est_mean / x$gamma_dd_mean
    
    gamma_zdsch_sens <- 0.5 * (x$Mn2_mean + x$Mn3_mean) * x$gamma_alpha_mean * 
      x$gamma_ztair_mean * x$gamma_zskt2_mean * x$gamma_est_mean / x$gamma_dd_mean
      
    gamma_zskt2_sens <- 0.5 * (x$Mn2_mean + x$Mn3_mean) * x$gamma_alpha_mean * 
      x$gamma_ztair_mean * x$gamma_zdsch_mean * x$gamma_est_mean / x$gamma_dd_mean
      
    gamma_est_sens <- 0.5 * (x$Mn2_mean + x$Mn3_mean) * x$gamma_alpha_mean * 
      x$gamma_ztair_mean * x$gamma_zdsch_mean * x$gamma_zskt2_mean / x$gamma_dd_mean
    
    # Derivative of lambda with respect to density dependence in gamma
    # using the chain rule
    # D(lam_exp, "gamma") = 0.5 * (Mn2_mean + Mn3_mean)
    # 
    # D(gamma_exp, "dd") =  -Z1 * Z2 * Z3 * Z4 * Z5 / dd ^ 2
    # 
    # D(lam_exp, "gamma") * D(gamma_exp, "dd") = -0.5 * (M2_mean + M3_mean) * Z1 * Z2 * Z3 * Z4 * Z5 / dd ^ 2
    
    gamma_dd_sens <- -0.5 * (x$Mn2_mean + x$Mn3_mean) * x$gamma_alpha_mean * 
      x$gamma_ztair_mean * x$gamma_zdsch_mean * x$gamma_zskt2_mean * 
      x$gamma_est_mean / x$gamma_dd_mean ^ 2
    
    # Derivative of lambda with respect to linear predictors of phi1, phi2, and
    # phi3 using the chain rule. Symbolic computation below uses phi1, but the
    # same applies to phi2 and phi3
    # D(lam_exp, "phi1") = Mn1
    # 
    # D(phi_exp, "lp") = (exp(-lp) / (1 + exp(-lp))^2) / dd
    # 
    # D(lam_exp, "phi1") * D(phi_exp, "lp") = Mn1 * (exp(-lp) / (1 + exp(-lp))^2) / dd
    
    phi1_lp_sens <- x$Mn1_mean * (exp(-x$phi1_lpc_mean) / (1 + exp(-x$phi1_lpc_mean)) ^ 2) / x$phi1_dd_mean
    phi2_lp_sens <- x$Mn2_mean * (exp(-x$phi2_lpc_mean) / (1 + exp(-x$phi2_lpc_mean)) ^ 2) / x$phi2_dd_mean
    phi3_lp_sens <- x$Mn3_mean * (exp(-x$phi3_lpc_mean) / (1 + exp(-x$phi3_lpc_mean)) ^ 2) / x$phi3_dd_mean
    
    # Derivative of lambda with respect to density dependence of phi1, phi2, and
    # phi3 using the chain rule. Symbolic computation below uses phi1, but the
    # same applies to phi2 and phi3
    # D(lam_exp, "phi1") = Mn1
    # 
    # D(phi_exp, "dd") = -(1 / (1 + exp(-lp))) / dd ^ 2
    # 
    # D(lam_exp, "phi1") * D(gamma_exp, "dd") = Mn1 * -(1 / (1 + exp(-lp))) / dd ^ 2
    
    phi1_dd_sens <- x$Mn1_mean * -(1 / (1 + exp(-x$phi1_lpc_mean))) / x$phi1_dd_mean ^ 2
    phi2_dd_sens <- x$Mn2_mean * -(1 / (1 + exp(-x$phi2_lpc_mean))) / x$phi2_dd_mean ^ 2
    phi3_dd_sens <- x$Mn3_mean * -(1 / (1 + exp(-x$phi3_lpc_mean))) / x$phi3_dd_mean ^ 2
    
    # Derivative of lambda with respect to linear predictors of psi1 and psi2 
    # using the chain rule. Symbolic computation below uses psi1, but the same 
    # applies to psi2. As the derivative of lam_exp with respect to psi1 or psi2
    # is zero, then the derivative of lam_exp with respect to the linear predictor
    # will be zero.
    # 
    # D(lam_exp, "psi1") = 0
    # 
    # D(psi_exp, "lp") = exp(-lp) / (1 + exp(-lp))^2
    # 
    # D(lam_exp, "psi1") * D(psi_exp, "lp") = 0 * exp(-lp) / (1 + exp(-lp))^2 = 0
    
    psi1_lp_sens <- psi2_lp_sens <- 0
    
    # The derivative with respect to Mn will have a term that is common to all 
    # three population sizes. We'll assign it to pd_M to avoid clutter below.

    pd_Mn <-  x$phi1_mean * x$Mn1_mean + 
              0.5 * x$gamma_mean * (x$Mn2_mean + x$Mn3_mean) + 
              x$phi2_mean * x$Mn2_mean + 
              x$phi3_mean * x$Mn3_mean
    
    # D(lam_exp, "Mn1")
    Mn1_sens <- x$phi1_mean - pd_Mn
    
    # D(lam_exp, "Mn2")
    Mn2_sens <- (0.5 * x$gamma_mean + x$phi2_mean) - pd_Mn
    
    # D(lam_exp, "Mn3")
    Mn3_sens <- (0.5 * x$gamma_mean + x$phi3_mean) - pd_Mn
    

    # Note that the sensitivities for lp have to be repeated as many times
    # as there are components in the linear predictor
    if (ltre == "var") {
      sensitivities <- c(gamma_alpha_sens, gamma_ztair_sens, gamma_zdsch_sens, 
                         gamma_zskt2_sens, gamma_est_sens, gamma_dd_sens, 
                         rep(phi1_lp_sens, 4), rep(phi2_lp_sens, 4), rep(phi3_lp_sens, 4),
                         phi1_dd_sens, phi2_dd_sens, phi3_dd_sens,
                         rep(psi1_lp_sens, 4), rep(psi2_lp_sens, 4),
                         Mn1_sens, Mn2_sens, Mn3_sens)                    
    } else if (ltre == "diff"){
      sensitivities <- as.data.frame(
        cbind(gamma_alpha_sens, gamma_ztair_sens, gamma_zdsch_sens, 
              gamma_zskt2_sens, gamma_est_sens, gamma_dd_sens, 
              matrix(phi1_lp_sens, length(phi1_lp_sens), 4), 
              matrix(phi2_lp_sens, length(phi2_lp_sens), 4), 
              matrix(phi3_lp_sens, length(phi3_lp_sens), 4), 
              phi1_dd_sens, phi2_dd_sens, phi3_dd_sens,
              # Use length of normalized vector, as we set sensitivities
              # of psi as a scalar (0)
              matrix(psi1_lp_sens, length(Mn1_sens), 4), 
              matrix(psi2_lp_sens, length(Mn2_sens), 4), 
              Mn1_sens, Mn2_sens, Mn3_sens))
      
      names(sensitivities) <- c("gamma_alpha", "gamma_tair", "gamma_dsch", "gamma_skt2", 
                                "gamma_est", "gamma_dd", 
                                "phi1_tair", "phi1_skt2", "phi1_skt1", "phi1_est",
                                "phi2_tair", "phi2_skt2", "phi2_skt1", "phi2_est",
                                "phi3_tair", "phi3_skt2", "phi3_skt1", "phi3_est",
                                "phi1_dd", "phi2_dd", "phi3_dd",
                                "psi1_tair", "psi1_skt2", "psi1_skt1", "psi1_est",
                                "psi2_tair", "psi2_skt2", "psi2_skt1", "psi2_est",
                                "M1", "M2", "M3")
    }
    
    return(sensitivities)
  }
} 