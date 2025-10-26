## ---------------------------
##
## Script name: 05b_EMmod.R
##
## Purpose of script: Fit em mod to get source probabilities
##
## Author: Trent VanHawkins
##
## Date Created: 2025-06-06
##
##
## ---------------------------
## ---------------------------
# --- Logit and Inverse Logit ---
logit <- function(p) log(p / (1 - p))
inv_logit <- function(x) 1 / (1 + exp(-x))

# Dispersal function ------------------------------------------------------
kappa_inner_sum_backward <- function(par_mat, y_prev, wind_mat, dist_mat, d0 = 0.01, derivative = FALSE, group_id) {
  
  n <- length(y_prev)
  K <- length(unique(group_id))

  # Compute constant terms
  dist_shifted <- dist_mat + d0
  log_dist <- log(dist_shifted)
  kappa <- par_mat['kappa',]

  #Initialize matrix for storage
  dispersal_mat <- matrix(0, nrow = n, ncol = K)

  #Assume infection came from group k
  for (k in 1:K){
    group_k <- which(group_id == k)
    y_prev_k <- y_prev[group_k]
    kappa_k <- kappa[k]

    #Compute distance kernel
    dist_kernel_k <- dist_shifted^(-kappa_k)
    if(derivative){
      dist_kernel_k <- -dist_kernel_k * log_dist
    }

    #Compute spread matrix
    wind_dist <- wind_mat * dist_kernel_k
    spread_k <- wind_dist[,group_k, drop = F] %*% y_prev_k

    dispersal_mat[, k] <- spread_k
  }

return(dispersal_mat)  # n x K matrix of dispersal from each group
}

# Mean function -----------------------------------------------------------
get_mu <- function(par, y_prev, dispersal) {

#Compute Autoinfection
auto <- y_prev * (1 - y_prev)
  
#Compute linear predictor
eta <- par[['beta']] + par[['delta']] * auto + par[['gamma']] * dispersal

#invers-logit
mu  <- inv_logit(eta)

return(pmin(pmax(mu, 1e-6), 1 - 1e-6))
}

# Compute Likelihood ---------------------------------------------------------
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

# Function to compute Q (expected log-likelihood)
Q_fun <- function(y, mu_mat, phi, p_mat, pi_vec) {
  #Compute weighted likelihood matrix
  lik_mat <- apply(mu_mat, 2, function(x) loglik_zibeta(y, x, phi, sum = F, log = T))

  #compute Q-val
  Q_val <- sum(p_mat * lik_mat) + sum(p_mat * pi_vec)

  return(-Q_val)
}

# E-step ------------------------------------------------------------------
e_step <- function(y, mu_mat, phi, prior){
  
  # Initialize weighted likelihood matrix
  wl_mat <- prior * apply(mu_mat, 2, function(x) loglik_zibeta(y, x, phi, sum = F, log = F))
  
  #Compute posterior probabilities
  p_mat <- wl_mat / rowSums(wl_mat)
  
  return(p_mat)
}


# M-step ------------------------------------------------------------------
m_step_grad <- function(par, y_current, y_prev, dispersal, dispersal_grad, p_vec, mu_vec) {
  beta  <- par["beta"]
  delta <- par["delta"]
  gamma <- par["gamma"]
  kappa <- par["kappa"]
  log_phi   <- par["phi"]
  phi <- exp(log_phi)

  non_zero <- which(y_current > 0)

  # Create model matrix
  auto_vec <- y_prev * (1 - y_prev)

  #Subset the non-zero part
  y_current <- y_current[non_zero]
  auto_vec <- auto_vec[non_zero]
  dispersal <- dispersal[non_zero]
  dispersal_grad <- dispersal_grad[non_zero]
  p_vec <- p_vec[non_zero]
  mu_vec <- mu_vec[non_zero]
  
  # Compute derivatives
  y_star <- logit(y_current)
  mu_star <- digamma(mu_vec * phi) - digamma((1 - mu_vec) * phi)
  weight <- phi * (y_star - mu_star) * mu_vec * (1 - mu_vec)
  
  # Weight everything by p_mat
  d_beta  <- sum(p_vec * weight)
  d_delta <- sum(p_vec * weight * auto_vec)
  d_gamma <- sum(p_vec * weight * dispersal)
  d_kappa <- sum(p_vec * weight * (-gamma) * dispersal_grad)
  
  
  # phi derivative (as scalar sum over i and k)
  d_phi <- sum(p_vec * ((digamma(phi)) - mu*(digamma(mu*phi)) - (1 - mu)*digamma((1- mu)*phi) + mu*log(y_current) + (1 - mu)*log(1 - y_current) * phi))
  
  -c(beta = d_beta, delta = d_delta, gamma = d_gamma, kappa = d_kappa, phi = d_phi)
}
