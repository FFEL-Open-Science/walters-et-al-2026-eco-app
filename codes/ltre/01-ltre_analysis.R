
# Initial setup ----
library(tidyverse)
source("codes/ltre/00-function_ltre_var_vr.R")
source("codes/ltre/00-function_ltre_var_ep.R")
source("codes/ltre/00-function_ltre_diff.R")
source("codes/ltre/00-function_ltre_lamg.R")

#out <- readRDS("outputs/spl_dat_beta-binomial_dd_all_10e6_2e6.rds")
out <- readRDS("outputs/spl_dat_nmix_beta-binomial_dd_all.rds")

spl <- as.data.frame(do.call(rbind, out$spl))

dat <- out$dat

ns <- nrow(spl)

# Run function to compute transient LTRE to determine the contribution of 
# vital rates and population structure to the variance of the growth rate lambda

ltre_var_vr_out <- ltre_var_vr(spl, 1988:2022, ns)
saveRDS(ltre_var_vr_out, file = "outputs/output_ltre_var_vr.rds")

ltre_var_ep_out <- ltre_var_ep(spl, dat, 1988:2022, ns)
saveRDS(ltre_var_ep_out, file = "outputs/output_ltre_var_ep.rds")

ltre_diff_out <- ltre_diff(spl, dat, 1988:2022, ns)
saveRDS(ltre_diff_out, file = "outputs/output_ltre_diff.rds")

ltre_lamg_out <- ltre_lamg(spl, dat, 1988:2022, 1988:2004, 2005:2021, 21, ns)
saveRDS(ltre_lamg_out, file = "outputs/output_ltre_lamg.rds")
