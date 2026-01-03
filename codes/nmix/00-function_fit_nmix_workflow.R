# This function is needed to run MCMC chain in parallel, but can also be used to
# run a single chain and see the progress.

fit_nmix_workflow <- function(seed, dat, cnt, ni, nb, nt, nc) {
  # Function arguments:
    # dat: a list with the data to send to Nimble 
    # cnt: a list with constants to send to Nimble
    # ni: the number of iterations of the MCMC chains
    # nb: the number of iterations to discard as burnin
    # nt: the thinning rate to apply to the MCMC chains
    # nc: the number of MCMC chains
  
  library(nimbleEcology)
  library(truncnorm)
  source("codes/nmix/00-nmix_model.R") 

  set.seed(seed)
  
  # Function to generate initial values (it's a long one) ----
  initial_values <- function() {
    
    Y <- dat$Y
    
    # Simulate parameters for observation model
    beta0_p <- runif(cnt$sC, -1, 0)
    sig_p <- runif(1, 0, 0.5)
    eps_p <- rep(NA, cnt$T - 1)

    alpha_B <- beta_B <- p <- matrix(NA, nrow = cnt$sC, ncol = cnt$T - 1)

    zeta_p <- runif(2, 0.4, 0.6)
    z_p <- rbinom(2, 1, zeta_p[1:2])
    beta_p <- rep(NA, 2)
    beta_p[1] <- ifelse(z_p[1] == 0, 0, runif(1, 0, 0.5))
    beta_p[2] <- ifelse(z_p[2] == 0, 0, runif(1, -0.5, 0))
    theta_p <- runif(3, 0, 3)
    
    for (t in 1:(cnt$T - 1)) {
      eps_p[t] <- rnorm(1, 0, sig_p)

      for (c in 1:cnt$sC) {
        if (cnt$obs == "binomial" | cnt$obs == "beta-binomial") {
          p[c, t] <- plogis(beta0_p[c] +
                              z_p[1] * beta_p[1] * dat$zvis[t] +
                              z_p[2] * beta_p[2] * dat$zdis[t] +
                              eps_p[t])
        } else if (cnt$obs == "poisson") {
          p[c, t] <- exp(beta0_p[c] +
                           z_p[1] * beta_p[1] * dat$zvis[t] +
                           z_p[2] * beta_p[2] * dat$zdis[t] +
                           eps_p[t])
        }
      }
    }
    
    # Simulate number of fish available to be counted. Base initial abundance
    # available on the observed counts and simulated detection probability.
    N <- array(NA, dim = dim(Y))

    for (c in 1:cnt$sC) {
      N[, 1, c] <- round(Y[, 1, c] / p[c, 1])
    }

    # Simulate total population size based on a simulated value for beta_omega
    # and total distance traversed by snorkelers at a occasion
    beta_omega <- runif(1, 0, 2)
    
    M_tmp <- N[, 1, ] / exp(-beta_omega / dat$d[, 1])

    M <- M_int <- matrix(NA, nrow = cnt$sC, ncol = cnt$T + 1)

    M[, 1] <- round(apply(M_tmp, 2, max))

    ## Initial values for phi ----
    beta0_phi <- runif(cnt$sC, 3, 4)

    zeta_phi <- z_phi <- beta_phi <- matrix(NA, nrow = cnt$sC, ncol = 3)
    
    beta_phi_DD <- rep(NA, 3)
    
    for (c in 1:cnt$sC) {
      for (b in 1:3) {
        zeta_phi[c, b] <- runif(1, 0.4, 0.6)
        z_phi[c, b] <- rbinom(1, 1, zeta_phi[c, b])
        beta_phi[c, b] <- ifelse(z_phi[c, b] == 0, 0, runif(1, -0.5, 0.5))
      }
      if (cnt$dd == "phi" | cnt$dd == "all") {
        beta_phi_DD[c] <- runif(1, 0.0001, 0.00025)
      }
    }

    sig_phi <-  runif(cnt$sC, 0, 0.5)

    # Initial values for gamma ----
    alpha_gamma <- runif(1, 2.5, 3.5)
    zeta_gamma <- z_gamma <- beta_gamma <- rep(NA, 3)
    
    for (b in 1:3) {
      zeta_gamma[b] <- runif(1, 0.4, 0.6)
      z_gamma[b] <- rbinom(1, 1, zeta_gamma[b])
      beta_gamma[b] <- ifelse(z_gamma[b] == 0, 0, runif(1, -0.5, 0.5))
    }
    
    if (cnt$dd == "gamma" | cnt$dd == "all") {
      beta_gamma_DD <- runif(1, 0.0001, 0.00025)
    }

    sig_gamma <- runif(1, 0, 0.1)
    
    ## Initial values for psi ----
    beta0_psi <- c(runif(1, -0.5, 0.5), runif(1, -2, -1))
    
    zeta_psi <- z_psi <- beta_psi <- rep(NA, 3)
    for (b in 1:3) {
      zeta_psi[b] <- runif(1, 0.4, 0.6)
      z_psi[b] <- rbinom(1, 1, zeta_psi[b])
      beta_psi[b] <- ifelse(z_psi[b] == 0, 0, runif(1, -0.1, 0.1))
    }

    sig_psi <-  runif(1, 0, 0.5)

    gamma <- eps_gamma <- rep(NA, cnt$T)
    phi <- eps_phi <- matrix(NA, nrow = cnt$sC, ncol = cnt$T)
    S <- R <- matrix(NA, nrow = cnt$sC, ncol = cnt$T + 1)
    eps_psi <- rep(NA, cnt$sT)
    psi <- matrix(NA, nrow = cnt$sC, ncol = cnt$sT)
    
    for (t in 1:cnt$sT) {
      psi[3, t] <- 0
      for (c in 1:(cnt$sC - 1)) {
        eps_psi[t] <- rnorm(1, 0, sig_psi)

        psi[c, t] <- plogis(beta0_psi[c] +
                            z_psi[1] * beta_psi[1] * dat$ztair[t] +
                            z_psi[2] * beta_psi[2] * dat$zskt2[t] +
                            z_psi[3] * beta_psi[3] * dat$zskt1[t] +
                            eps_psi[t])
      }
    }
        
    for (t in 2:(cnt$T + 1)) {

      eps_gamma[t - 1] <- rnorm(1, 0, sig_gamma)
      
      if (cnt$dd == "gamma" | cnt$dd == "all") {
        gamma[t - 1] <- exp(z_gamma[1] * beta_gamma[1] * dat$ztair[t + 8 - 1] +
                            z_gamma[2] * beta_gamma[2] * dat$zdsch[t - 1] +
                            z_gamma[3] * beta_gamma[3] * dat$zskt2[t + 8 - 1] +
                            eps_gamma[t - 1]) * alpha_gamma /
          (1 + beta_gamma_DD * (sum(M[2:3, t - 1]) / 2))
      } else {
        gamma[t - 1] <- exp(z_gamma[1] * beta_gamma[1] * dat$ztair[t + 8 - 1] +
                              z_gamma[2] * beta_gamma[2] * dat$zdsch[t - 1] +
                              z_gamma[3] * beta_gamma[3] * dat$zskt2[t + 8 - 1] +
                              eps_gamma[t - 1]) * alpha_gamma
      }

      for (c in 1:cnt$sC) {

        eps_phi[c, t - 1] <- rnorm(1, 0, sig_phi[c])
        
        if (cnt$dd == "phi" | cnt$dd == "all") {
          phi[c, t - 1] <- plogis(beta0_phi[c] +
                                    z_phi[c, 1] * beta_phi[c, 1] * dat$ztair[t + 8 - 1] +
                                    z_phi[c, 2] * beta_phi[c, 2] * dat$zskt2[t + 8 - 1] +
                                    z_phi[c, 3] * beta_phi[c, 3] * dat$zskt1[t + 8 - 1] +
                                    eps_phi[c, t - 1]) / (1 + beta_phi_DD[c] * M[c, t - 1])
        } else {
          phi[c, t - 1] <- plogis(beta0_phi[c] +
                                    z_phi[c, 1] * beta_phi[c, 1] * dat$ztair[t + 8 - 1] +
                                    z_phi[c, 2] * beta_phi[c, 2] * dat$zskt2[t + 8 - 1] +
                                    z_phi[c, 3] * beta_phi[c, 3] * dat$zskt1[t + 8 - 1] +
                                    eps_phi[c, t - 1])
        }

      }

      R[1, t] <- rpois(1, gamma[t - 1] * (sum(M[2:3, t - 1] / 2)))
      
      for (c in 2:cnt$sC) {
        R[c, t] <- rbinom(1, M[c - 1, t - 1], phi[c - 1, t - 1] * psi[c - 1, t + 8 - 1])
      }
      
      for (c in 1:cnt$sC) {
        S[c, t] <- rbinom(1, M[c, t - 1], phi[c, t - 1] * (1 - psi[c, t + 8 - 1]))
      }
      
      for (c in 1:cnt$sC) {

        M[c, t] <- R[c, t] + S[c, t]

      }
    }

    for (t in 2:(cnt$T + 1)) {
      for (c in 1:cnt$sC){
        M_int[c, t - 1] <- round((1 - dat$w[t - 1]) * M[c, t - 1] + dat$w[t - 1] * M[c, t])
      }
    }

    for (t in 2:(cnt$T - 1)) {
      for (c in 1:cnt$sC) {
        N[1:cnt$J[t], t, c] <- rbinom(cnt$J[t], M_int[c, cnt$sy[t]], exp(-beta_omega / dat$d[, t]))
      }

    }

    M[, 2:(cnt$T + 1)] <- NA

    ini_lst <- list(exp_factor = runif(1, 1, 2),
                    M = M,
                    S = S,
                    R = R,
                    beta_omega = beta_omega,
                    N = N,

                    beta0_p = beta0_p,
                    beta_p = beta_p,
                    z_p = z_p,
                    zeta_p = zeta_p,
                    sig_p = sig_p,
                    theta_p = theta_p,

                    alpha_gamma = alpha_gamma,
                    beta_gamma = beta_gamma,
                    sig_gamma = sig_gamma,

                    beta0_phi = beta0_phi,
                    beta_phi = beta_phi,
                    sig_phi = sig_phi,
                    
                    beta0_psi = beta0_psi,
                    beta_psi = beta_psi,
                    sig_psi = sig_psi,
                    psi = psi)
    
    if(cnt$dd == "phi") {
      ini_lst$beta_phi_DD <- beta_phi_DD
    } else if (cnt$dd == "gamma") {
      ini_lst$beta_gamma_DD <- beta_gamma_DD
    } else if (cnt$dd == "all") {
      ini_lst$beta_phi_DD <- beta_phi_DD
      ini_lst$beta_gamma_DD <- beta_gamma_DD
    }
    
    return(ini_lst)
    
  }
  
  # Parameters to monitor ----
  pars <- c("beta0_phi", "beta_phi", "z_phi", "zeta_phi",
            "sig_phi", "eps_phi", "phi",

            "alpha_gamma", "beta_gamma", "z_gamma", "zeta_gamma",
            "sig_gamma", "eps_gamma", "gamma",
            
            "beta0_p", "beta_p", "z_p", "zeta_p", "sig_p", "eps_p", "p",

            "beta_omega", "M", "N", "S", "R", "lam_M", "eta_M", "mu_M", "M_int",

            "Y_rep", "Y_exp", "Y_var",

            "beta0_psi", "beta_psi", "z_psi", "zeta_psi", "sig_psi", "psi", 
            "eps_psi", "C", "C_rep", "Psi", "Pnum_tmp", "Pden_tmp", "Pres")

  if(cnt$dd == "phi") {
    pars <- c(pars, "beta_phi_DD")
  } else if (cnt$dd == "gamma") {
    pars <- c(pars, "beta_gamma_DD")
  } else if (cnt$dd == "all") {
    pars <- c(pars, "beta_phi_DD", "beta_gamma_DD")
  }
  
  if (cnt$obs == "beta-binomial") {
    pars <- c(pars, c("theta_p", "alpha_B", "beta_B"))
  }
  
  # Send to Nimble ----
  model <- nimbleModel(nmix_model_code, constants = cnt, data = dat,
                       inits = initial_values(), calculate = FALSE)
  
  comp_model <- compileNimble(model)

  conf_mcmc <- configureMCMC(model, monitors = pars, useConjugacy = FALSE,
                             onlySlice = TRUE)
  
  target <- c("beta_phi", "beta_gamma", "beta_p", "beta_psi")

  indicator <- c("z_phi", "z_gamma", "z_p", "z_psi")

  configureRJ(conf = conf_mcmc,
              targetNodes = target,
              indicatorNodes = indicator,
              control = list(mean = NULL, scale = NULL))
  
  mcmc <- buildMCMC(conf_mcmc)
  
  comp_mcmc <- compileNimble(mcmc, project = comp_model, resetFunctions = TRUE)
  
  samples <- runMCMC(comp_mcmc, inits = initial_values, niter = ni, nburnin = nb, 
                       thin = nt, nchains = nc)
  
  return(samples)
  
}