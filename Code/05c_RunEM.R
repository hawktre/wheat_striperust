## ---------------------------
##
## Script name: 05c_RunEM.R
##
## Purpose of script: Run the EM algorithm using functions from 05b
##
## Author: Trent VanHawkins
##
## Date Created: 2025-06-11
##
##
## ---------------------------

## view outputs in non-scientific notation

options(scipen = 6, digits = 4) 

## ---------------------------

## load up the packages we will need:  (uncomment as required)

library(tidyverse)
library(here)

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
library(tidyverse)
library(here)

# Read in the data --------------------------------------------------------
em_dat <- readRDS(here("DataProcessed/experimental/em_dat.rds"))

## Read in the needed functions
source(here("Code/05b_EMgradfun.R"))


# Specifications for the model -------------------------------------------------
dat <- em_dat$stripe_4$A1$visit4
n_groups <- length(unique(dat$group_id))

# Fit the model -----------------------------------------------------------

  
## Figure out which ones are zero
pi <- mean(dat$y_cur == 0)
non_zero <- which(dat$y_cur != 0)

# Parameters
max_em_iter <- 20
tol <- 1e-6
theta <- dat$inits
q_track <- numeric(max_em_iter)
theta_track <- vector("list", max_em_iter)

# Setup for inputs
em_args <- list(
  y         = dat$y_cur[non_zero],
  y_prev    = dat$y_prev[non_zero],
  wind_matrix = dat$wind_mat[non_zero, non_zero],
  dist_matrix = dat$dist_mat[non_zero, non_zero],
  group_id  = dat$group_id[non_zero]
)

# Initial E-step
p_mat <- do.call(e_step, c(list(par = theta, prior = 1/n_groups), em_args))

# Loop
for (em_iter in seq_len(max_em_iter)) {
  
  # M-step
  fit <- try(
    optim(
      par = theta,
      fn = wrapped_obj,
      gr = mstep_grad_em,
      method = "BFGS",
      control = list(maxit = 1000, reltol = 1e-8),
      p_mat = p_mat,
      !!!em_args
    ),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error")) {
    warning(paste("EM step", em_iter, "failed"))
    break
  }
  
  theta_new <- fit$par
  theta_track[[em_iter]] <- theta_new
  
  # Update E-step
  p_mat <- do.call(e_step, c(list(par = theta_new, prior = 1/n_groups), em_args))
  
  # Compute Q
  mu_mat <- do.call(get_mu, c(list(par = theta_new), em_args))
  q_val <- Q_fun(y = em_args$y, mu_mat = mu_mat, phi = theta_new["phi"], p_mat = p_mat)
  q_track[em_iter] <- -q_val
  
  # Convergence check
  if (em_iter > 1 && abs(q_track[em_iter] - q_track[em_iter - 1]) < tol) {
    message("EM converged at iteration ", em_iter)
    break
  }
  
  theta <- theta_new
}

# Trim results
q_track <- q_track[1:em_iter]
theta_track <- theta_track[1:em_iter]

# Wrap up
em_summary <- tibble(
  plot_id = "A1",
  visit = "visit4",
  em_iters = em_iter,
  converged = em_iter < max_em_iter,
  final_neg_loglik = q_track[em_iter],
  final_theta = list(theta_new)
)








