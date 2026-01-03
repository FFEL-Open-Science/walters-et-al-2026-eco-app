compute_log_diff <- function(par, type, dir_diff, T) {
  
  if(!(type %in% c("mu", "sig"))) {
    stop("type must be either 'mu' or 'sig'.")
  }

  if(!(dir_diff %in% c(12, 21))) {
    stop("dir must be either 12 or 21.")
  }
  
  if (type == "mu") {
    if (dir_diff == 12) {
      log_diff <- log(mean(par[1:T])) - log(mean(par[(T + 1):(T * 2)]))
    }

    if (dir_diff == 21) {
      log_diff <- log(mean(par[(T + 1):(T * 2)])) - log(mean(par[1:T]))
    }
  }

  if (type == "sig") {
    if (dir_diff == 12) {
      log_diff <- log(sd(par[1:T])) - log(sd(par[(T + 1):(T * 2)]))
    }
    
    if (dir_diff == 21) {
      log_diff <- log(sd(par[(T + 1):(T * 2)])) - log(sd(par[1:T]))
    }
  }

  return(log_diff)
  
}