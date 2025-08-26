## ---------------------------
##
## Script name: 02c_ModelFun.R
##
## Purpose of script: Contains all necessary modeling functions 
##
## Author: Trent VanHawkins
##
## Date Created: 2025-08-20
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
library(here)
source(here("Code/02a_ForwardGradFun.R"))

# Forward Fits -------------------------------------------------------------
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
    #Generate initial values for current kappa in grid search
    init_theta <- initialize_theta(y_cur = intensity[non_zero], 
                                   y_prev = intensity_prev[non_zero], 
                                   wind_mat = wind[non_zero, non_zero], 
                                   dist_mat = dist[non_zero, non_zero], 
                                   kappa_try = as.numeric(kappa),
                                   d_0 = 0.01)
    
    fit <- tryCatch(
      optim(
        par = init_theta,
        fn = neg_loglik,
        gr = neg_grad,
        method = "BFGS",
        control = list(maxit = 5000, reltol = 1e-6),
        y_current = intensity[non_zero],
        y_prev = intensity_prev[non_zero],
        wind_matrix = wind[non_zero, non_zero],
        dist_matrix = dist[non_zero, non_zero]
      ),
      error = function(e) {
        message(sprintf("forward_fit failed: %s [blk=%s, trt=%s, vst=%s, kappa=%s]",
                        conditionMessage(e), blk, trt, vst, kappa))
        return(NULL)
      }
    )
    
    if (!is.null(fit)) {
      results_list[[length(results_list) + 1]] <- list(
        block = blk,
        treat = trt,
        visit = vst,
        iters = fit$counts[2],
        init_kappa = as.numeric(kappa),
        neg_loglik = fit$value,
        converged = fit$convergence == 0,
        theta = list(fit$par),
        pi = pi
      )
    }
    
  }
  
  if (length(results_list) > 0) {
    visit_dt <- rbindlist(results_list)
    visit_dt <- visit_dt[neg_loglik > -1e6 & converged == T]
    if (nrow(visit_dt) > 0) {
      return(visit_dt[which.min(neg_loglik)])
    }
  }
  
  return(NULL)
}


# Forward Model Fitted Values ---------------------------------------------
get_fitted <- function(blk, trt, vst, par, mod_dat, d0 = 0.01) {
  y_prev <- mod_dat$intensity[,blk,trt,as.numeric(vst) - 1]
  wind_matrix <- mod_dat$wind[,,blk,trt,vst]
  dist_matrix <- mod_dat$dist
  beta  <- par["beta"]
  delta <- par["delta"]
  gamma <- par["gamma"]
  kappa <- par["kappa"]
  phi   <- par["phi"]
  
  
  #Compute Covariates
  ## Dispersal 
  dispersal <- kappa_inner_sum(y_prev, wind_matrix, dist_matrix, d0, kappa)
  
  n <- length(y_prev)
  S <- ncol(dispersal)
  
  ## Auto-infection
  y_vec <- y_prev * (1 - y_prev)
  
  
  eta_mat <- beta + delta * y_vec + gamma * dispersal
  mu_mat  <- inv_logit(eta_mat)
  
  return(pmin(pmax(mu_mat, 1e-6), 1 - 1e-6))  # Clip for stability
}


# Forward Model Deviance Residuals ----------------------------------------
compute_deviance_resid <- function(blk, trt, vst, par, pi, fitted, data) {

  y <- data$intensity[,blk, trt, vst]
  mu_hat <- fitted
  phi_hat <- par[["phi"]] 
  
  is_zero <- y == 0
  is_pos <- !is_zero
  
  ll_fit <- numeric(length(y))
  ll_sat <- numeric(length(y))
  
  # Zero values
  ll_fit[is_zero] <- log(pi)
  ll_sat[is_zero] <- log(pi)  # same since model fits perfectly
  
  # Positive values
  ll_fit[is_pos] <- log(1 - pi) + loglik_beta(y[is_pos], mu_hat[is_pos], phi_hat, sum = FALSE)
  ll_sat[is_pos] <- log(1 - pi) + loglik_beta(y[is_pos], y[is_pos], phi_hat, sum = FALSE)
  
  sqrt_term <- pmax(0, 2 * (ll_sat - ll_fit))
  sign(y - (1 - pi) * mu_hat) * sqrt(sqrt_term)
}



