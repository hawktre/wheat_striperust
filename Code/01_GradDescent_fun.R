## ---------------------------
##
## Script name: 01_GradDescent_R
##
## Purpose of script: Functions to run gradient descent go estimate parameter vector for this project.
##
## Author: Trent VanHawkins
##
## Date Created: 2025-04-07
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


# --- Initialization ---
initialize_theta <- function() {
  list(
    beta = 0.1,
    delta = 0.1,
    gamma = 0.1,
    kappa = 1.0,
    phi = 10.0
  )
}

# --- Logit and Inverse Logit ---
logit <- function(p) log(p / (1 - p))
inv_logit <- function(x) 1 / (1 + exp(-x))

# --- special function for kappa gradient ---
kappa_inner_sum <- compute_kappa_inner_sums <- function(y_prev, wind_matrix, dist_matrix, d0, kappa) {
  n <- length(y_prev)
  
  # Compute shifted distance matrix and log
  dist_shifted <- dist_matrix + d0
  log_dist <- log(dist_shifted)
  kernel <- dist_shifted^(-kappa)
  
  # Broadcast y_prev across rows
  y_mat <- matrix(y_prev, nrow = n, ncol = n, byrow = TRUE)
  
  # Elementwise product
  spread_matrix <- y_mat * wind_matrix * kernel * log_dist
  
  # Zero out diagonal to exclude j = i
  diag(spread_matrix) <- 0
  
  # Row sums = vector of kappa inner sums for each i
  rowSums(spread_matrix)
}


# --- Main Gradient Descent Function ---
fit_beta_model <- function(y_current, y_prev, wind_matrix, dist_matrix,
                           alpha = 1e-3, max_iter = 5000, tol = 1e-6, d0 = 0.01) {
  
  n <- length(y_current)
  theta <- initialize_theta()
  param_history <- list()
  
  for (iter in 1:max_iter) {
    # --- Compute linear predictor and mu_i ---
    eta <- numeric(n)
    mu <- numeric(n)
    
    ## Compute inner sum for dispersal (vector)
    for (i in 1:n) {
      spread_sum <- numeric(n)
      for (j in 1:n) {
        if (j != i) {
          spread_sum[i] <- spread_sum[i] + y_prev[j] * wind_matrix[i, j] / (dist_matrix[i, j] + d0)^theta$kappa
        }
      }
      eta[i] <- theta$beta + theta$delta * y_prev[i] + theta$gamma * spread_sum[i]
      mu[i] <- inv_logit(eta[i])
    }
    # --- Compute weights ---
    y_star <- logit(y_current)
    mu_star <- digamma(mu*theta$phi) - digamma((1-mu)*theta$phi)
    theta_weight <- theta$phi * (y_star - mu_star) * mu * (1-mu)
    # --- Compute gradients (see 00a_ModelingConsiderations.pdf)---
    
    grad_beta <- sum(theta_weight)
    grad_delta <- sum(theta_weight * y_prev)
    grad_gamma <- sum(theta_weight * spread_sum)
    grad_kappa <- sum(theta_weight * -theta$gamma*kappa_inner_sum(y_prev, wind_matrix, dist_matrix, d0, theta$kappa))
    grad_phi <- sum(
      digamma(theta$phi)
      -mu*digamma(mu*theta$phi)
      -(1-mu)*digamma(theta$phi*(1-mu))
      + mu*log(y_current)
      +(1-mu)*log(1- y_current)
      )
    
    # --- Update parameters ---
    theta$beta  <- theta$beta  - alpha * grad_beta
    theta$delta <- theta$delta - alpha * grad_delta
    theta$gamma <- theta$gamma - alpha * grad_gamma
    theta$kappa <- theta$kappa - alpha * grad_kappa
    theta$phi   <- theta$phi   - alpha * grad_phi
    
    # --- Check convergence ---
    grad_norm <- sqrt(grad_beta^2 + grad_delta^2 + grad_gamma^2 + grad_kappa^2 + grad_phi^2)
    if (grad_norm < tol) {
      cat("Converged at iteration:", iter, "\n")
      break
    }
    
    # Optional: Track history
    param_history[[iter]] <- unlist(theta)
  }
  
  return(list(theta = theta, history = param_history))
}




