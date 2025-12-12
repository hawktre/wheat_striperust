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
forward_fit <- function(blk, trt, vst, mod_dat, kappa_try) {
  
  ## Extract needed data
  intensity <- mod_dat$intensity[,blk,trt,vst]
  intensity_prev <- mod_dat$intensity[,blk,trt,as.numeric(vst)-1]
  wind <- mod_dat$wind[,,blk,trt,vst]
  dist <- mod_dat$dist

  ## Extract Needed Model Data
  non_zero <- which(intensity > 0)
  
  ## Set up strage for temproary results
  results_list <- list()
  
  for (kappa in kappa_try) {
    #Generate initial values for current kappa in grid search
    init_theta <- initialize_theta(y = intensity, 
                                   y_prev = intensity_prev, 
                                   wind_mat = wind, 
                                   dist_mat = dist, 
                                   kappa = as.numeric(kappa),
                                   d_0 = 0.01)
    
    fit <- tryCatch(
      optim(
        par = init_theta,
        fn = neg_loglik,
        gr = neg_grad,
        method = "BFGS",
        control = list(maxit = 5000, 
                       reltol = 1e-6),
        y_current = intensity,
        y_prev = intensity_prev,
        wind_matrix = wind,
        dist_matrix = dist
      ),
      error = function(e) {
        message(sprintf("forward_fit failed: %s [blk=%s, trt=%s, vst=%s, kappa=%s]",
                        conditionMessage(e), blk, trt, vst, kappa))
        return(NULL)
      }
    )
    
    if (!is.null(fit)) {
      g <- neg_grad(
        fit$par,
        y_current = intensity,
        y_prev = intensity_prev,
        wind_matrix = wind,
        dist_matrix = dist)
      
      results_list[[length(results_list) + 1]] <- list(
        block = blk,
        treat = trt,
        visit = vst,
        iters = fit$counts[2],
        init_kappa = as.numeric(kappa),
        neg_loglik = fit$value,
        grad_norm = sqrt(sum(g^2)),
        converged = fit$convergence == 0,
        theta = list(fit$par),
        alpha = mean(intensity == 0)
      )
    }
    
  }
  
  if (length(results_list) > 0) {
    visit_dt <- rbindlist(results_list)
    visit_dt <- visit_dt[neg_loglik > -1e6 & converged == T]
    return(visit_dt[which.min(neg_loglik)])

  }else{
    print("No good fits. Quitting forward fit.")
  }
  
}



# Forward Model Fitted Values ---------------------------------------------
get_fitted <- function(blk, trt, vst, par, alpha, mod_dat, d0 = 0.01) {
  y_prev <- mod_dat$intensity[,blk,trt,as.numeric(vst) - 1]
  wind_matrix <- mod_dat$wind[,,blk,trt,vst]
  dist_matrix <- mod_dat$dist
  beta  <- par["beta"]
  delta <- par["delta"]
  gamma <- par["gamma"]
  kappa <- par["kappa"]

  
  #Compute Covariates
  ## Dispersal 
  dispersal <- kappa_inner_sum(y_prev, wind_matrix, dist_matrix, d0, kappa)
  
  n <- length(y_prev)
  S <- ncol(dispersal)
  
  ## Auto-infection
  auto_infection <- y_prev * (1 - y_prev)
  
  
  eta <- beta + delta * auto_infection + gamma * dispersal
  mu  <- inv_logit(eta)
  fitted <- (1-alpha) * mu
  
  return(fitted)  
}


# Forward Model Deviance Residuals ----------------------------------------
compute_deviance_resid <- function(blk, trt, vst, par, alpha, fitted, data) {
  y <- data$intensity[, blk, trt, vst]
  mu_hat <- fitted
  phi_hat <- exp(par[["phi"]])
  
  # Separate zero and nonzero cases
  nonzero <- which(y > 0)
  
  # Fitted log-likelihood for all data (ZIBeta handles zeros fine)
  ll_fit <- loglik_zibeta(y, mu_hat, phi_hat, sum = FALSE)
  
  # Saturated log-likelihood:
  ll_sat <- ll_fit
  if (length(nonzero) > 0) {
    # For nonzero observations, mu_tilde = y_t (but clip for numerical safety)
    y_nz <- pmin(pmax(y[nonzero], 1e-6), 1 - 1e-6)
    ll_sat[nonzero] <- loglik_zibeta(y_nz, y_nz, phi_hat, sum = FALSE)
  }
  
  # Deviance residuals
  sqrt_term <- pmax(0, 2 * (ll_sat - ll_fit))
  r_d <- sign(y - mu_hat) * sqrt(sqrt_term)
  
  return(r_d)
}



