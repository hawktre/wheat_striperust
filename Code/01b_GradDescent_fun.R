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

# --- Logit and Inverse Logit ---
logit <- function(p) log(p / (1 - p))
inv_logit <- function(x) 1 / (1 + exp(-x))

# --- special function for kappa gradient ---
kappa_inner_sum <- function(y_prev, wind_matrix, dist_matrix, d0, kappa, derivative = F) {
  n <- length(y_prev)
  
  # Compute shifted distance matrix and log
  dist_shifted <- dist_matrix + d0
  log_dist <- log(dist_shifted)
  kernel <- dist_shifted^(-kappa)
  
  # Broadcast y_prev across rows
  y_mat <- matrix(y_prev, nrow = n, ncol = n, byrow = TRUE)
  
  if(derivative){
    # Elementwise product
    spread_matrix <- y_mat * wind_matrix * kernel * log_dist
  }
  else{
    spread_matrix <- y_mat * wind_matrix * kernel
  }
  
  # Zero out diagonal to exclude j = i
  diag(spread_matrix) <- 0
  
  # Row sums = vector of kappa inner sums for each i
  rowSums(spread_matrix)
}

# --- Initialization ---

initialize_theta <- function(y_cur, y_prev, wind_mat, dist_mat, d_0, kappa_try) {

 #Construct the design matrix
  X1 <- y_prev
  # Get dispersal component and construct design matrix
  X2 <- map(kappa_try, ~ kappa_inner_sum(
    y_prev = y_prev,   # or provide fixed values
    wind_matrix = wind_mat,
    dist_matrix = dist_mat,
    d0 = d_0,
    derivative = FALSE,
    kappa = .x
  ))
  
  X_mat <- map(X2, ~cbind(X1 = X1, X2 = .x))
  
  logit_y <- logit(y_cur + 1e-3)
  fits <- map(X_mat, ~ lm(logit_y ~ ., data = as.data.frame(.x)))
  
  phi <- map(fits, function(fit) {

    #total samples
    n <- length(y_prev)
    
    # degrees of freedom
    df <- n - length(coef(fit)) + 1
    
    # predicted values 
    mu <- inv_logit(fitted(fit))

    #SSR of model fit
    ssr <- sum(resid(fit)^2)
    
    #Weight for residual variance
    g_prime <- 1/(mu * (1-mu))
    
    #denominator for sigma
    denom <- df*(g_prime^2)
    
    #compute sigma_hat
    sigma_hat <- ssr/denom
    
    #return phi
    mean((mu*(1-mu))/sigma_hat) - 1
  })
  
  # Extract coefficients
  coefs <- map(fits, coef)

  # Compose result as list of named vectors
  theta_list <- pmap(list(coefs, kappa_try, phi), function(coef_vec, kappa_val, phi_val) {
    c(
      beta  = coef_vec[["(Intercept)"]],
      delta = coef_vec[["X1"]],
      gamma = coef_vec[["X2"]],
      kappa = kappa_val,
      phi   = phi_val
    )
  })

  return(theta_list)
  
}

# --- Main Gradient Descent Function ---
fit_beta_model <- function(y_current, y_prev, wind_matrix, dist_matrix, param_init,
                           alpha = 1e-3, max_iter = 5000, tol = 1e-6, d0 = 0.01) {
  browser()
  n <- length(y_current)
  theta <- param_init
  param_history <- list()
  
  for (iter in 1:max_iter) {
    #Compute dispersal component
    dispersal <- kappa_inner_sum(y_prev, wind_matrix, dist_matrix, d0, theta[['kappa']])
    # --- Compute linear predictor and mu_i ---
    eta <- theta[['beta']]+ theta[['delta']]*y_prev + theta[['gamma']]*dispersal
    mu <- inv_logit(eta)
    
    # --- Compute weights ---
    y_star <- logit(y_current)
    mu_star <- digamma(mu*theta[['phi']]) - digamma((1-mu)*theta[['phi']])
    theta_weight <- theta[['phi']] * (y_star - mu_star) * mu * (1-mu)
    
    # --- Compute gradients (see 00a_ModelingConsiderations.pdf)---
    
    grad_beta <- sum(theta_weight)
    grad_delta <- sum(theta_weight * y_prev)
    grad_gamma <- sum(theta_weight * dispersal)
    grad_kappa <- sum(theta_weight * -theta[['gamma']]*kappa_inner_sum(y_prev, wind_matrix, dist_matrix, d0, theta[['kappa']], derivative = T))
    grad_phi <- sum(
      digamma(theta[['phi']])
      -mu*digamma(mu*theta[['phi']])
      -(1-mu)*digamma(theta[['phi']]*(1-mu))
      + mu*log(y_current)
      +(1-mu)*log(1- y_current)
      )
    
    # --- Update parameters ---
    theta[['beta']]  <- theta[['beta']]  - alpha * grad_beta
    theta[['delta']] <- theta[['delta']] - alpha * grad_delta
    theta[['gamma']] <- theta[['gamma']] - alpha * grad_gamma
    theta[['kappa']]<- theta[['kappa']] - alpha * grad_kappa
    theta[['phi']]   <- theta[['phi']]   - alpha * grad_phi
    
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