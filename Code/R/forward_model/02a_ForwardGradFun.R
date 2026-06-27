## ---------------------------
##
## Script name: 04a_GradDescent-Hurdle_fun.R
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

initialize_theta <- function(y, y_prev, wind_mat, dist_mat, d_0, kappa) {
  
  non_zero <- which(y > 0)
  
  X1 <- y_prev * (1 - y_prev)
  X2 <- kappa_inner_sum(y_prev = y_prev, 
                        wind_matrix = wind_mat,
                        dist_matrix = dist_mat,
                        d0 = d_0,
                        derivative = FALSE,
                        kappa = kappa)
  
  mod_df <- data.frame(y = logit(pmin(pmax(y[non_zero], 1e-5), 1 - 1e-5)),
                       X1 = X1[non_zero],
                       X2 = X2[non_zero])
  
  fit <- lm(y ~ X1 + X2, data = mod_df)
  coef_vec <- coef(fit)
  
  # Total samples for variance approximation
  n <- length(y[non_zero])
  df <- n - length(coef_vec) + 1
  mu <- inv_logit(fitted(fit))
  ssr <- sum(resid(fit)^2)
  g_prime <- 1 / (pmax(mu * (1 - mu), 1e-6))
  sigma_hat <- ssr / (df * (g_prime^2))
  
  phi <- mean((mu * (1 - mu)) / pmax(sigma_hat, 1e-6)) - 1
  phi <- pmin(pmax(phi, 1e-3), 100) # Keep initial variance scale anchored Safely
  
  # DAMP THE INITIAL LAUNCHPAD VALUES
  # Instead of letting OLS pass huge parameters to optim, cap them to a realistic scale
  theta <- c(
    beta  = as.numeric(coef_vec[["(Intercept)"]]),
    delta = pmin(pmax(as.numeric(coef_vec[["X1"]]), -2), 5),   # Bounded launch area
    gamma = pmin(pmax(as.numeric(coef_vec[["X2"]]), -2), 5),   # Bounded launch area
    kappa = kappa,                                             # Starts exactly at grid step
    phi   = log(phi)
  )
  return(theta)
}

# --- function for log-likelihood (see Ospina & Ferrari) ---
loglik_zibeta <- function(y, mu, phi, sum = TRUE, log = TRUE) {
  alpha <- mean(y == 0)
  
  # Calculate log-likelihood for all elements simultaneously
  ll <- ifelse(y == 0, 
               log(alpha), 
               log(1 - alpha) + dbeta(y, shape1 = mu * phi, shape2 = (1 - mu) * phi, log = TRUE))
  
  if (sum) {
    out <- sum(ll)
    if (!log) out <- exp(out)
  } else {
    out <- if (log) ll else exp(ll)
  }
  return(out)
}

# --- function for negative log-likelihood ---
neg_loglik <- function(par, y_current, y_prev, wind_matrix, dist_matrix, d0 = 0.01) {
  # Gather Parameters
  beta  <- par["beta"]
  delta <- par["delta"]
  gamma <- par["gamma"]
  kappa <- par["kappa"]
  log_phi <- par["phi"]

  #Compute fitted values 
  dispersal <- kappa_inner_sum(y_prev, wind_matrix, dist_matrix, d0, kappa)
  auto <- y_prev * (1-y_prev)
  eta <- beta + delta * auto + gamma * dispersal
  mu  <- inv_logit(eta)
  phi <- exp(log_phi)

  -loglik_zibeta(y_current, mu, phi)
}

# --- function for gradients ---
neg_grad <- function(par, y_current, y_prev, wind_matrix, dist_matrix, d0 = 0.01) {
  beta  <- par["beta"]
  delta <- par["delta"]
  gamma <- par["gamma"]
  kappa <- par["kappa"]
  log_phi <- par["phi"]
  phi <- exp(log_phi)

  #Only Non-zero observations
  non_zero <- which(y_current > 0)

  #Compute dispersal terms
  dispersal <- kappa_inner_sum(y_prev, wind_matrix, dist_matrix, d0, kappa)[non_zero]
  dispersal_grad <- kappa_inner_sum(y_prev, wind_matrix, dist_matrix, d0, kappa, derivative = TRUE)[non_zero]
  
  #Compute fitted values (and clamp if necessary)
  y_current <- y_current[non_zero]
  y_prev <- y_prev[non_zero]
  auto <- y_prev * (1-y_prev)
  eta <- beta + delta * auto + gamma * dispersal
  mu  <- inv_logit(eta)
  mu <- pmin(pmax(mu, 1e-6), 1 - 1e-6)

  #Compute weights that show up in all derivatives
  y_star <- logit(y_current)
  mu_star <- digamma(mu * phi) - digamma((1 - mu) * phi)
  weight <- (y_star - mu_star) * mu * (1 - mu)
  
  d_beta  <-  phi * sum(weight)
  d_delta <-  phi * sum(weight * auto)
  d_gamma <-  phi * sum(weight * dispersal)
  d_kappa <-  phi * sum(weight * (-gamma) * dispersal_grad)
  d_phi   <-  sum((digamma(phi)) - mu*(digamma(mu*phi)) - (1 - mu)*digamma((1- mu)*phi) + mu*log(y_current) + (1 - mu)*log(1 - y_current)) * phi
  
  # Return negative gradients
  -c(beta = d_beta, delta = d_delta, gamma = d_gamma, kappa = d_kappa, phi = d_phi)
}
## ---------------------------
## Script name: 02a_ForwardGradFun.R
## Purpose: Functions to run gradient descent on natural scale (except log_phi).
## Author: Trent VanHawkins
## ---------------------------

logit <- function(p) log(p / (1 - p))
inv_logit <- function(x) 1 / (1 + exp(-x))

kappa_inner_sum <- function(y_prev, wind_matrix, dist_matrix, d0, kappa, derivative = FALSE) {
  n <- length(y_prev)
  dist_shifted <- dist_matrix + d0
  log_dist <- log(dist_shifted)
  kernel <- dist_shifted^(-kappa)
  
  y_mat <- matrix(y_prev, nrow = n, ncol = n, byrow = TRUE)
  
  if(derivative){
    spread_matrix <- y_mat * wind_matrix * kernel * log_dist
  } else {
    spread_matrix <- y_mat * wind_matrix * kernel
  }
  diag(spread_matrix) <- 0
  return(rowSums(spread_matrix))
}

# --- Initialization (Back to natural scale) ---
initialize_theta <- function(y, y_prev, wind_mat, dist_mat, d_0, kappa) {
  non_zero <- which(y > 0)
  X1 <- y_prev * (1 - y_prev)
  X2 <- kappa_inner_sum(y_prev = y_prev, wind_matrix = wind_mat, dist_matrix = dist_mat, d0 = d_0, derivative = FALSE, kappa = kappa)
  
  mod_df <- data.frame(y = logit(pmin(pmax(y[non_zero], 1e-5), 1 - 1e-5)), X1 = X1[non_zero], X2 = X2[non_zero])
  fit <- lm(y ~ X1 + X2, data = mod_df)
  
  n <- length(y[non_zero])
  df <- n - length(coef(fit)) + 1
  mu <- inv_logit(fitted(fit))
  ssr <- sum(resid(fit)^2)
  g_prime <- 1 / (pmax(mu * (1 - mu), 1e-6))
  sigma_hat <- ssr / (df * (g_prime^2))
  
  phi <- mean((mu * (1 - mu)) / pmax(sigma_hat, 1e-6)) - 1
  phi <- max(phi, 1e-6)
  coef_vec <- coef(fit)
  
  theta <- c(
    beta  = as.numeric(coef_vec[["(Intercept)"]]),
    delta = as.numeric(coef_vec[["X1"]]),
    gamma = as.numeric(coef_vec[["X2"]]),
    kappa = kappa,
    phi   = log(phi) # Only phi remains on the log scale
  )
  return(theta)
}

# --- Log-likelihood with clipping ---
loglik_zibeta <- function(y, mu, phi, sum = TRUE, log = TRUE) {
  alpha <- mean(y == 0)
  mu <- pmin(pmax(mu, 1e-6), 1 - 1e-6) # Guardrail clipping
  
  ll <- ifelse(y == 0, 
               log(alpha), 
               log(1 - alpha) + dbeta(y, shape1 = mu * phi, shape2 = (1 - mu) * phi, log = TRUE))
  if (sum) {
    out <- sum(ll)
    if (!log) out <- exp(out)
  } else {
    out <- if (log) ll else exp(ll)
  }
  return(out)
}

# --- Negative log-likelihood ---
neg_loglik <- function(par, y_current, y_prev, wind_matrix, dist_matrix, d0 = 0.01) {
  beta  <- par["beta"]
  delta <- par["delta"]
  gamma <- par["gamma"]
  kappa <- par["kappa"]
  log_phi <- par["phi"]

  dispersal <- kappa_inner_sum(y_prev, wind_matrix, dist_matrix, d0, kappa)
  auto <- y_prev * (1 - y_prev)
  eta <- beta + delta * auto + gamma * dispersal
  mu  <- inv_logit(eta)
  phi <- exp(log_phi)

  -loglik_zibeta(y_current, mu, phi)
}

# --- Gradient function (Natural parameters + log-scale phi chain rule) ---
neg_grad <- function(par, y_current, y_prev, wind_matrix, dist_matrix, d0 = 0.01) {
  beta  <- par["beta"]
  delta <- par["delta"]
  gamma <- par["gamma"]
  kappa <- par["kappa"]
  log_phi <- par["phi"]
  phi <- exp(log_phi)

  non_zero <- which(y_current > 0)
  dispersal <- kappa_inner_sum(y_prev, wind_matrix, dist_matrix, d0, kappa)[non_zero]
  dispersal_grad <- kappa_inner_sum(y_prev, wind_matrix, dist_matrix, d0, kappa, derivative = TRUE)[non_zero]
  
  y_current <- y_current[non_zero]
  y_prev <- y_prev[non_zero]
  auto <- y_prev * (1 - y_prev)
  eta <- beta + delta * auto + gamma * dispersal
  mu  <- inv_logit(eta)
  mu  <- pmin(pmax(mu, 1e-6), 1 - 1e-6) # Guardrail clipping

  y_star <- logit(y_current)
  mu_star <- digamma(mu * phi) - digamma((1 - mu) * phi)
  weight <- (y_star - mu_star) * mu * (1 - mu)
  
  d_beta  <-  phi * sum(weight)
  d_delta <-  phi * sum(weight * auto)
  d_gamma <-  phi * sum(weight * dispersal)
  d_kappa <-  phi * sum(weight * (-gamma) * dispersal_grad)
  
  # Natural scale derivative for phi
  d_phi_natural <- sum((digamma(phi)) - mu * (digamma(mu * phi)) - (1 - mu) * digamma((1 - mu) * phi) + mu * log(y_current) + (1 - mu) * log(1 - y_current))
  d_log_phi     <- d_phi_natural * phi # Chain rule wrapper for log_phi

  -c(beta = d_beta, delta = d_delta, gamma = d_gamma, kappa = d_kappa, phi = d_log_phi)
}