## ---------------------------
##
## Script name: 04_ModelFit_Hurdle_Logistic.R
##
## Purpose of script: Fit model with logistic autoinfection term
##
## Author: Trent VanHawkins
##
## Date Created: 2025-05-16
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

## load up the packages we will need:  (uncomment as required)
library(here)
library(purrr)
library(data.table)
library(MASS)
# Read in the data --------------------------------------------------------
mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat_arrays.rds"))
source(here("Code/02a_GradDescentFun.R"))


# Set up indices ----------------------------------------------------------
blocks <- dimnames(mod_dat$intensity)[["block"]]
treats <- dimnames(mod_dat$intensity)[["treat"]]
visits <- dimnames(mod_dat$intensity)[["visit"]]
kappa_try <- dimnames(mod_dat$inits)[["kappa_try"]]

combos <- expand.grid(
  block = blocks,
  treat = treats,
  visit = visits[2:length(visits)],
  stringsAsFactors = FALSE
)
# Set up data -------------------------------------------------------------
dist <- mod_dat$dist #doesn't change

# Fit the model -----------------------------------------------------------

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
    visit_dt <- rbindlist(results_list)
    visit_dt <- visit_dt[neg_loglik > -1e5]
    if (nrow(visit_dt) > 0) {
      return(visit_dt[which.min(neg_loglik)])
    }
  }
  
  return(NULL)
}

library(furrr)
plan(multisession, workers = 4)

start <- Sys.time()
free_fits <- future_pmap(
  combos,
  ~forward_fit(..1, ..2, ..3, mod_dat, dist, kappa_try)
) %>% rbindlist()
end <- Sys.time()
runtime <- difftime(end, start, units = "mins")  # could be "mins", "hours", etc.
message("Runtime = ", round(runtime, 2), " minutes")
plan(sequential)


# Get fitted values -------------------------------------------------------
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

free_fits$fitted <- pmap(free_fits, ~get_fitted(blk = ..1, trt = ..2, vst = ..3, par = ..8, mod_dat))
saveRDS(free_fits, here("DataProcessed/results/forward_model/forward_fits.rds"))

