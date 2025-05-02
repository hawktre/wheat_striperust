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
  return(rowSums(spread_matrix))
  
}

# --- Initialization ---

initialize_theta <- function(y_cur, y_prev, wind_mat, dist_mat, d_0, kappa_try) {

 #Construct the design matrix
  X1 <- y_prev
  # Get dispersal component and construct design matrix
  X2 <- map(kappa_try, ~kappa_inner_sum(
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
# --- function for log-likelihood ---
loglik_beta <- function(y, mu, phi, sum = T) {
  ll <- lgamma(phi) -
    lgamma(mu * phi) -
    lgamma((1 - mu) * phi) +
    (mu * phi - 1) * log(y) +
    ((1 - mu) * phi - 1) * log(1 - y)
  if (sum) {
    return(sum(ll))
  }else{
    return(ll)
  }
}

# --- function for negative log-likelihood ---
neg_loglik <- function(par, y_current, y_prev, wind_matrix, dist_matrix, d0 = 0.01) {
  beta  <- par["beta"]
  delta <- par["delta"]
  gamma <- par["gamma"]
  kappa <- par["kappa"]
  phi   <- par["phi"]
  
  dispersal <- kappa_inner_sum(y_prev, wind_matrix, dist_matrix, d0, kappa)
  eta <- beta + delta * y_prev + gamma * dispersal
  mu  <- inv_logit(eta)

  -loglik_beta(y_current, mu, phi)
}

# --- function for gradients ---
neg_grad <- function(par, y_current, y_prev, wind_matrix, dist_matrix, d0 = 0.01) {
  beta  <- par["beta"]
  delta <- par["delta"]
  gamma <- par["gamma"]
  kappa <- par["kappa"]
  phi   <- par["phi"]
  
  dispersal <- kappa_inner_sum(y_prev, wind_matrix, dist_matrix, d0, kappa)
  eta <- beta + delta * y_prev + gamma * dispersal
  mu  <- inv_logit(eta)
  
  y_star <- logit(y_current)
  mu_star <- digamma(mu * phi) - digamma((1 - mu) * phi)
  weight <- (y_star - mu_star) * mu * (1 - mu)
  
  d_beta  <-  phi * sum(weight)
  d_delta <-  phi * sum(weight * y_prev)
  d_gamma <-  phi * sum(weight * dispersal)
  d_kappa <-  phi * sum(weight * (-gamma) * kappa_inner_sum(y_prev, wind_matrix, dist_matrix, d0, kappa, derivative = TRUE))
  d_phi   <-  sum(mu * (y_star - mu_star) + log(1 - y_current) - digamma((1 - mu) * phi) + digamma(phi))
  
  # Return negative gradients
  -c(beta = d_beta, delta = d_delta, gamma = d_gamma, kappa = d_kappa, phi = d_phi)
}