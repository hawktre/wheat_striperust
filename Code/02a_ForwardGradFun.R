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

  #Subset non-zero observations
  non_zero <- which(y > 0)
  
  #Construct the design matrix
  X1 <- y_prev * (1- y_prev)
  
  # Get dispersal component and construct design matrix
  X2 <- kappa_inner_sum(y_prev = y_prev, 
                        wind_matrix = wind_mat,
                        dist_matrix = dist_mat,
                        d0 = d_0,
                        derivative = FALSE,
                        kappa = kappa)
  
  mod_df <- data.frame(y = logit(y[non_zero]),
                       X1 = X1[non_zero],
                       X2 = X2[non_zero])
  
  
  fit <- lm(y~X1 + X2, data = mod_df)
  
    
  #total samples
  n <- length(y[non_zero])
  df <- n - length(coef(fit)) + 1
  mu <- inv_logit(fitted(fit))
  ssr <- sum(resid(fit)^2)
  
  #Weight for residual variance
  g_prime <- 1/(mu * (1-mu))
  
  #denominator for sigma
  denom <- df*(g_prime^2)
  
  #compute sigma_hat
  sigma_hat <- ssr/denom
  
  #return phi
  phi <- mean((mu*(1-mu))/sigma_hat) - 1
  phi <- max(phi, 1e-6)
  
  # Extract coefficients
  coef_vec <- coef(fit)
  
  # Compose result as list of named vectors
  theta <- c(
    beta  = coef_vec[["(Intercept)"]],
    delta = coef_vec[["X1"]],
    gamma = coef_vec[["X2"]],
    kappa = kappa,
    phi   = log(phi)
  )
  return(theta)
  
}

# --- function for log-likelihood (see Ospina & Ferrari) ---
loglik_zibeta <- function(y, mu, phi, sum = TRUE, log = TRUE) {
  
  # Gather Fixed Terms
  n <- length(y)

  # Estimate alpha (zero-inflation term)
  alpha <- mean(y == 0)

  # Prepare per-observation log-likelihood vector
  ll <- numeric(n)

  for(i in 1:length(y)){
    if(y[i] == 0){
      ll[i] <- log(alpha)
    } else {
      a <- mu[i] * phi
      b <- (1 - mu[i]) * phi
      ll[i] <- log(1 - alpha) +
      (lgamma(phi) - lgamma(a) - lgamma(b) ) +
      (a - 1) * log(y[i]) +
      (b - 1) * log(1 - y[i])
    }
  }

  if (sum) {
    out <- sum(ll)
    if (!log) out <- exp(out)   # exponentiate after sum
  } else {
    out <- if (log) ll else exp(ll)  # return vector
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
  mu <- pmin(pmax(mu, 1e-6), 1 - 1e-6)
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
