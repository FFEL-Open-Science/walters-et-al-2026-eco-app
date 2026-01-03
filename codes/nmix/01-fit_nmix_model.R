library(nimbleEcology)
library(nimble)
library(parallel)

source("codes/nmix/00-function_prepare_data.R")
source("codes/nmix/00-function_fit_nmix_workflow.R")
source("codes/nmix/00-nmix_model.R")

obs <- "beta-binomial" # alternatives are "binomial" or "poisson"
dd <- "all" # alternatives are "gamma" or "phi"
dat_cnt <- prepare_data()
dat_cnt$cnt$obs <- obs
dat_cnt$cnt$dd <- dd

# MCMC setup ----

parallel <- TRUE
ni <- 10e6
nb <- 2e6
nt <- (ni - nb) / 1000
nc <- 4

start <- Sys.time()

if (parallel == FALSE) {
  spl <- fit_nmix_workflow(seed = 1, dat = dat_cnt$dat, cnt = dat_cnt$cnt,
                           ni = ni, nb = nb, nt = nt, nc = nc)
} else {
  my_cluster <- makeCluster(nc)
  
  spl <- parLapply(cl = my_cluster, 
                   X = 1:nc, 
                   fun = fit_nmix_workflow,
                   dat = dat_cnt$dat, cnt = dat_cnt$cnt,
                   ni = ni, nb = nb, nt = nt, nc = 1)
  
  stopCluster(my_cluster)
}

end <- Sys.time()

saveRDS(list(spl = spl, dat = dat_cnt$dat, cnt = dat_cnt$cnt,
            rt = end - start),
       paste("outputs/spl_dat_nmix_", obs, "_dd_", dd, ".rds", sep = ""))

print(end - start)
