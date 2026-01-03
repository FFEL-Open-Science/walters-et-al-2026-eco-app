library(tidyverse)
library(nimble)
library(parallel)

source("codes/vbgf/00-vbgf_model.R")

# # Process age-length data ----

L <- read_csv("data/raw/age_length.csv") |>
  select(`1`:`9`) |>
  as.matrix() / 10 # convert to cm

I <- nrow(L)

A <- apply(L, 1, function(x) max(which(!is.na(x))))

dat <- list(L = L, A = A)

cnt <- list(I = I)

# Function to generate initial values ----

ini <- function(){list(
  L_inf = runif(1, 50, 60),
  a0 = runif(1, 0, 0.3),
  rho = runif(1, 0.4, 0.6),
  sig_L = runif(1, 0, 1),
  mu_K = runif(1, 0.1, 0.5),
  sig_K = runif(1, 0, 0.1))}

# Parameters to monitor ----

par <- c("L_inf", "a0", "rho", "sig_L", "mu_K", "sig_K", "K",
         "a10", "a30", "a50", "res_L", "L_rep", "res_L_rep")

# MCMC setup ----

ni <- 6e4
nb <- 3e4
nt <- (ni - nb) / 1000
nc <- 4

# Send to Nimble ----

spl <- nimbleMCMC(code = vbgf_model_code, constants = cnt,
                  data = dat, inits = ini, nchains = nc,
                  niter = ni, nburnin = nb, thin = nt,
                  monitors = par)

saveRDS(list(spl = spl, dat = dat, cnt = cnt), "outputs/spl_dat_vbgf.rds")
