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
    }else{
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

Q_fun <- function(par, y_current, y_prev, wind_matrix, dist_matrix, group_id, p_mat, d0 = 0.01) {
  #Compute mu assuming infection from each group
  mu_mat <- get_mu(par = par,
                   y_prev = y_prev,
                   wind_matrix = wind_matrix,
                   dist_matrix = dist_matrix,
                   d0 = d0,
                   group_id = group_id)
  
  # Extract phi from parameter vector
  phi <- par["phi"]

  #Compute weighted likelihood matrix
  lik_mat <- lik_zibeta_backward(y, mu_mat, phi, sum = F, log = T) #When y > 0
  
  #compute Q-val
  Q_val <- sum(p_mat * lik_mat)

  return(-Q_val)
}

# Dispersal function ------------------------------------------------------
kappa_inner_sum_backward <- function(y_prev, wind_matrix, dist_matrix, d0, kappa, 
                            derivative = FALSE, group_id) {
  n <- length(y_prev)
  groups <- sort(unique(group_id))
  S <- length(groups)
  
  # Compute shifted distance matrix and dispersal kernel
  dist_shifted <- dist_matrix + d0
  log_dist <- log(dist_shifted)
  kernel <- dist_shifted^(-kappa)
  
  # Element-wise dispersal term
  y_mat <- matrix(y_prev, nrow = n, ncol = n, byrow = TRUE)
  spread_matrix <- y_mat * wind_matrix * kernel
  
  if (derivative) {
    spread_matrix <- spread_matrix * log_dist
  }
  
  diag(spread_matrix) <- 0  # remove self-contribution
  
  # Initialize result matrix: one column per group s
  dispersal_mat <- matrix(0, nrow = n, ncol = S)
  
  for (s in seq_along(groups)) {
    group_mask_vec <- as.numeric(group_id == groups[s])  # length-n vector
    dispersal_mat[, s] <- spread_matrix %*% group_mask_vec  # matrix-vector product
  }
  
  return(dispersal_mat)  # n x S matrix of dispersal from each group
}


# Mean function -----------------------------------------------------------
get_mu <- function(par, y_prev, wind_matrix, dist_matrix, d0 = 0.01, group_id) {
  beta  <- par["beta"]
  delta <- par["delta"]
  gamma <- par["gamma"]
  kappa <- par["kappa"]
  phi   <- par["phi"]
  
  
  #Compute Covariates
  ## Dispersal 
  dispersal <- kappa_inner_sum_backward(y_prev, wind_matrix, dist_matrix, d0, kappa, group_id = group_id)
  
  n <- length(y_prev)
  S <- ncol(dispersal)
  
  ## Auto-infection
  y_vec <- y_prev * (1 - y_prev)
  y_mat <- matrix(y_vec, nrow = n, ncol = S)  # broadcast across columns
  
  
  eta_mat <- beta + delta * y_mat + gamma * dispersal
  mu_mat  <- inv_logit(eta_mat)
  
  return(pmin(pmax(mu_mat, 1e-6), 1 - 1e-6))  # Clip for stability
}

# E-step ------------------------------------------------------------------
e_step <- function(par, y_current, y_prev, wind_matrix, dist_matrix, d0 = 0.01, group_id, prior){
  browser()
  mu_mat <- get_mu(par = par,
                   y_prev = y_prev,
                   wind_matrix = wind_matrix,
                   dist_matrix = dist_matrix,
                   group_id = group_id,
                   d0 = d0)
  
  n <- length(y_current)
  S <- length(unique(group_id))
  
  # Make the prior a vector of length n
  if(length(prior) == 1){
    prior <- rep(prior, S)
  }
  
  # Initialize weighted likelihood matrix
  wl_mat <- prior * lik_zibeta_backward(y_current, mu_mat, par[['phi']], sum = F, log = F)
  
  #Compute posterior probabilities
  p_mat <- wl_mat / rowSums(wl_mat)
  
  return(p_mat)
}


# M-step ------------------------------------------------------------------
mstep_grad_em <- function(par, y_current, y_prev, wind_matrix, dist_matrix, d0 = 0.01, group_id, p_mat) {
 
  beta  <- par["beta"]
  delta <- par["delta"]
  gamma <- par["gamma"]
  kappa <- par["kappa"]
  phi   <- par["phi"]

  non_zero <- which(y_current > 0)

  # Create model matrices for each term (n x S)
  y_prev_term <- outer(y_prev * (1 - y_prev), rep(1, ncol(mu_mat)))[non_zero, ]
  dispersal_term <- kappa_inner_sum_backward(y_prev, wind_matrix, dist_matrix, d0, kappa, group_id = group_id)[non_zero, ]
  dispersal_term_deriv <- kappa_inner_sum_backward(y_prev, wind_matrix, dist_matrix, d0, kappa, derivative = TRUE, group_id = group_id)[non_zero, ]
  
  # Compute predicted mu matrix (n x S)
  mu_mat <- get_mu(par, y_prev, wind_matrix, dist_matrix, d0, group_id)[non_zero, ]
  
  # Compute derivatives
  y_star <- logit(y_current[non_zero])
  mu_star <- digamma(mu_mat * phi) - digamma((1 - mu_mat) * phi)
  weight <- phi * (y_star - mu_star) * mu_mat * (1 - mu_mat)
  
  
  
  
  # Weight everything by p_mat
  d_beta  <- sum(p_mat * weight)
  d_delta <- sum(p_mat * weight * y_prev_term)
  d_gamma <- sum(p_mat * weight * dispersal_term)
  d_kappa <- sum(p_mat * weight * (-gamma) * dispersal_term_deriv)
  
  
  # phi derivative (as scalar sum over i and k)
  d_phi <- sum(p_mat * (mu_mat * (y_star - mu_star) +
                          log(1 - y_current) -
                          digamma((1 - mu_mat) * phi) +
                          digamma(phi)))
  
  -c(beta = d_beta, delta = d_delta, gamma = d_gamma, kappa = d_kappa, phi = d_phi)
}
