## ---------------------------
##
## Script name: 02b_ForwardModelFun.R
##
## Purpose of script: Contains all necessary modeling functions
##                    Natural scale parameters, no Hessians, returning all fits per loop.
##
## Author: Trent VanHawkins
##
## Date Created: 2025-08-20
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

## view outputs in non-scientific notation
options(scipen = 6, digits = 4)

## ---------------------------
library(here)
library(data.table)
source(here("Code/R/forward_model/02a_ForwardGradFun.R"))

# Forward Fits -------------------------------------------------------------
forward_fit <- function(blk, trt, vst, mod_dat, kappa_try) {
  ## Extract needed data
  intensity <- mod_dat$intensity[, blk, trt, vst]
  intensity_prev <- mod_dat$intensity[, blk, trt, as.numeric(vst) - 1]
  wind <- mod_dat$wind[,, blk, trt, vst]
  dist <- mod_dat$dist

  ## Extract Needed Model Data
  non_zero <- which(intensity > 0)

  ## Set up storage for temporary results
  results_list <- list()

  for (kappa in kappa_try) {
    # Generate initial values (delta, gamma, kappa are natural scale; phi is log scale)
    init_theta <- initialize_theta(
      y = intensity,
      y_prev = intensity_prev,
      wind_mat = wind,
      dist_mat = dist,
      kappa = as.numeric(kappa),
      d_0 = 0.01
    )

    fit <- tryCatch(
      optim(
        par = init_theta,
        fn = neg_loglik,
        gr = neg_grad,
        method = "BFGS",
        control = list(maxit = 5000, reltol = 1e-6),
        hessian = FALSE, # Completely stop tracking Hessians
        y_current = intensity,
        y_prev = intensity_prev,
        wind_matrix = wind,
        dist_matrix = dist
      ),
      error = function(e) {
        message(sprintf(
          "forward_fit failed: %s [blk=%s, trt=%s, vst=%s, kappa=%s]",
          conditionMessage(e),
          blk,
          trt,
          vst,
          kappa
        ))
        return(NULL)
      }
    )

    if (!is.null(fit)) {
      g <- neg_grad(
        fit$par,
        y_current = intensity,
        y_prev = intensity_prev,
        wind_matrix = wind,
        dist_matrix = dist
      )

      # Extract optimized parameters and convert only phi back to natural scale
      natural_theta <- fit$par
      natural_theta["phi"] <- exp(natural_theta["phi"])

      results_list[[length(results_list) + 1]] <- list(
        block = blk,
        treat = trt,
        visit = vst,
        kappa_init = kappa,
        neg_loglik = fit$value,
        grad_norm = sqrt(sum(g^2)),
        converged = fit$convergence == 0,
        theta = list(natural_theta),
        alpha = mean(intensity == 0)
      )
    }
  }

  if (length(results_list) > 0) {
    visit_dt <- rbindlist(results_list)
    visit_dt <- visit_dt[converged == TRUE]
    return(visit_dt) # Returns all valid fits for inspection across kappa values
  } else {
    print("No good fits. Quitting forward fit.")
    return(NULL)
  }
}


# Forward Model Fitted Values ---------------------------------------------
get_fitted <- function(blk, trt, vst, par, alpha, mod_dat, d0 = 0.01) {
  y_prev <- mod_dat$intensity[, blk, trt, as.numeric(vst) - 1]
  wind_matrix <- mod_dat$wind[,, blk, trt, vst]
  dist_matrix <- mod_dat$dist
  
  beta <- par["beta"]
  delta <- par["delta"]
  gamma <- par["gamma"]
  kappa <- par["kappa"]

  # Compute Covariates
  ## Dispersal
  dispersal <- kappa_inner_sum(y_prev, wind_matrix, dist_matrix, d0, kappa)

  ## Auto-infection
  auto_infection <- y_prev * (1 - y_prev)

  eta <- beta + delta * auto_infection + gamma * dispersal
  mu <- inv_logit(eta)
  mu <- pmin(pmax(mu, 1e-6), 1 - 1e-6) # Guardrail clipping included
  
  fitted <- (1 - alpha) * mu

  return(fitted)
}


# Forward Model Deviance Residuals ----------------------------------------
compute_deviance_resid <- function(blk, trt, vst, par, alpha, fitted, data) {
  y <- data$intensity[, blk, trt, vst]
  mu_hat <- fitted
  
  # par comes back on the natural scale out of forward_fit
  phi_hat <- par[["phi"]]

  # Separate zero and nonzero cases
  nonzero <- which(y > 0)

  # Fitted log-likelihood for all data (ZIBeta handles zeros fine)
  ll_fit <- loglik_zibeta(y, mu_hat, phi_hat, sum = FALSE)

  # Saturated log-likelihood:
  ll_sat <- ll_fit
  if (length(nonzero) > 0) {
    # For nonzero observations, mu_tilde = y_t (clipped for numerical safety)
    y_nz <- pmin(pmax(y[nonzero], 1e-6), 1 - 1e-6)
    ll_sat[nonzero] <- loglik_zibeta(y_nz, y_nz, phi_hat, sum = FALSE)
  }

  # Deviance residuals
  sqrt_term <- pmax(0, 2 * (ll_sat - ll_fit))
  r_d <- sign(y - mu_hat) * sqrt(sqrt_term)

  return(r_d)
}