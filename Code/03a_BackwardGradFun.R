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
kappa_inner_sum_backward <- function(par, y_prev, wind_mat, dist_mat, d0 = 0.01, derivative = FALSE, group_id) {
  
  n <- length(y_prev)
  groups <- sort(unique(group_id))
  K <- length(groups)
  
  # Compute shifted distance matrix and dispersal kernel
  dist_shifted <- dist_mat + d0
  log_dist <- log(dist_shifted)
  kernel_mat <- dist_shifted^(-par[['kappa']])
  
  # Element-wise dispersal term
  y_mat <- matrix(y_prev, nrow = n, ncol = n, byrow = TRUE)
  spread_mat <- y_mat * wind_mat * kernel_mat
  
  if (derivative) {
    spread_mat <- spread_mat * log_dist
  }
  
  diag(spread_matrix) <- 0  # remove self-contribution
  
  # Initialize result matrix: one column per group s
  dispersal_mat <- matrix(0, nrow = n, ncol = K)
  
  for (k in seq_along(groups)) {
    group_mask_vec <- as.numeric(group_id == groups[k])  # length-n vector
    dispersal_mat[, k] <- spread_matrix %*% group_mask_vec  # matrix-vector product
  }
  
  return(dispersal_mat)  # n x K matrix of dispersal from each group
}
 # Mean function -----------------------------------------------------------
get_mu <- function(par, auto, dispersal) {

  #Compute linear predictor
  eta <- par[['beta']] + par[['delta']] * auto + gamma * dispersal
  
  #invers-logit
  mu  <- inv_logit(eta_mat)
  
  return(pmin(pmax(mu_mat, 1e-8), 1 - 1e-8))
}

# Compute Likelihood ---------------------------------------------------------
loglik_zibeta <- function(y_vec, mu_vec, phi, sum = TRUE, log = TRUE) {
  
  # Gather Fixed Terms
  n <- length(y_vec)

  # Estimate alpha (zero-inflation term)
  alpha <- mean(y_vec == 0)

  # Prepare per-observation log-likelihood vector
  ll <- numeric(n)

  for(i in 1:length(y_vec)){
    if(y_vec[i] == 0){
      ll[i] <- log(alpha)
    }else{
      a <- mu_vec[i] * phi
      b <- (1 - mu_vec[i]) * phi
      ll[i] <- log(1 - alpha) +
      (lgamma(phi) - lgamma(a) - lgamma(b) ) +
      (a - 1) * log(y_vec[i]) +
      (b - 1) * log(1 - y_vec[i])
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
Q_fun <- function(y_vec, mu_mat, phi, p_mat, pi_vec) {
  #Compute weighted likelihood matrix
  lik_mat <- apply(mu_mat, 2, function(x) loglik_zibeta(y_vec, x, phi, sum = F, log = T))
  
  #compute Q-val
  Q_val <- sum(p_mat * lik_mat) + sum(p_mat * pi_vec)

  return(-Q_val)
}

# E-step ------------------------------------------------------------------
e_step <- function(y_vec, mu_mat, phi, prior_vec){
  
  # Initialize weighted likelihood matrix
  wl_mat <- prior_vec * apply(mu_mat, 2, function(x) loglik_zibeta(y_vec, x, phi, sum = F, log = F))
  
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
  phi   <- par["phi"]

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
  d_phi <- sum(p_vec * (digamma(phi)) - 
    mu_vec*(digamma(mu_vec*phi)) - 
    (1 - mu_vec)*digamma((1- mu_vec)*phi) + 
    mu_vec*log(y_current) + 
    (1 - mu_vec)*log(1 - y_current))
  
  -c(beta = d_beta, delta = d_delta, gamma = d_gamma, kappa = d_kappa, phi = d_phi)
}
