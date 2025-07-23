## ---------------------------
##
## Script name: 06a_StripeSimFunc.R
##
## Purpose of script:
##
## Author: Trent VanHawkins
##
## Date Created: 2025-07-21
##
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## ---------------------------

## view outputs in non-scientific notation

options(scipen = 6, digits = 4) 

## ---------------------------


# Forward Model Fit -------------------------------------------------------

forward_fit <- function(blk, trt, vst, mod_dat, dist, kappa_try) {
  ## Extract Needed Model Data
  intensity <- mod_dat$intensity[, blk, trt, vst]
  intensity_prev <- mod_dat$intensity[, blk, trt, as.numeric(vst) - 1]
  wind <- mod_dat$wind[, , blk, trt, vst]
  pi <- mean(intensity == 0)
  non_zero <- which(intensity != 0)
  
  ## Set up strage for temproary results
  results_list <- list()
  
  for (kappa in kappa_try) {
    init_theta <- mod_dat$inits[, kappa, blk, trt, vst] #Current inits
    
    fit <- try(
      optim(
        par = init_theta,
        fn = neg_loglik,
        gr = neg_grad,
        method = "BFGS",
        hessian = TRUE,
        control = list(maxit = 5000, reltol = 1e-8),
        y_current = intensity[non_zero],
        y_prev = intensity_prev[non_zero],
        wind_matrix = wind[non_zero, non_zero],
        dist_matrix = dist[non_zero, non_zero]
      ),
      silent = TRUE
    )
    
    if (!inherits(fit, "try-error")) {
      se <- sqrt(diag(ginv(fit$hessian)))
      results_list[[length(results_list) + 1]] <- list(
        block = blk,
        treat = trt,
        visit = vst,
        iters = fit$counts[1],
        init_kappa = unname(init_theta["kappa"]),
        neg_loglik = fit$value,
        converged = fit$convergence == 0,
        theta = list(fit$par),
        theta_se = list(se),
        pi = pi
      )
    }
  }
  
  if (length(results_list) > 0) {
    visit_dt <- data.table::rbindlist(results_list)
    visit_dt <- visit_dt[neg_loglik > -1e5]
    if (nrow(visit_dt) > 0) {
      return(visit_dt[which.min(neg_loglik)])
    }
  }
  
  return(NULL)
}

# Backward Model Fit ------------------------------------------------------


