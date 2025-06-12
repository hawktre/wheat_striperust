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

## view outputs in non-scientific notation

options(scipen = 6, digits = 4) 

## ---------------------------
# --- Logit and Inverse Logit ---
logit <- function(p) log(p / (1 - p))
inv_logit <- function(x) 1 / (1 + exp(-x))

# Compute Likelihood ---------------------------------------------------------
lik_beta <- function(y, mu_mat, phi, sum = TRUE, log = TRUE) {
  y <- matrix(y, nrow = nrow(mu_mat), ncol = ncol(mu_mat))  
  phi <- matrix(phi, nrow = nrow(mu_mat), ncol = ncol(mu_mat))  
  
  val <- lgamma(phi) -
    lgamma(mu_mat * phi) -
    lgamma((1 - mu_mat) * phi) +
    (mu_mat * phi - 1) * log(y) +
    ((1 - mu_mat) * phi - 1) * log(1 - y)
  
  if (!log) {
    val <- exp(val)
  }
  
  if (sum) {
    return(rowSums(val))
  } else {
    return(val)
  }
}

Q_fun <- function(y, mu_mat, phi, p_mat) {
  lik_mat <- lik_beta(y, mu_mat, phi, sum = FALSE)
  Q_val <- sum(p_mat * lik_mat)
  return(-Q_val)
}

wrapped_obj <- function(par, y_current, y_prev, wind_matrix, dist_matrix, group_id, p_mat, d0 = 0.01) {
  mu_mat <- get_mu(par = par,
                   y_current = y_current,
                   y_prev = y_prev,
                   wind_matrix = wind_matrix,
                   dist_matrix = dist_matrix,
                   d0 = d0,
                   group_id = group_id)
  phi <- par["phi"]
  Q_fun(y = y_current, mu_mat = mu_mat, phi = phi, p_mat = p_mat)
}



# Dispersal function ------------------------------------------------------
kappa_inner_sum <- function(y_prev, wind_matrix, dist_matrix, d0, kappa, 
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
get_mu <- function(par, y_current, y_prev, wind_matrix, dist_matrix, d0 = 0.01, group_id) {
  beta  <- par["beta"]
  delta <- par["delta"]
  gamma <- par["gamma"]
  kappa <- par["kappa"]
  phi   <- par["phi"]
  
  dispersal <- kappa_inner_sum(y_prev, wind_matrix, dist_matrix, d0, kappa, group_id = group_id)
  n <- length(y_prev)
  S <- ncol(dispersal)
  
  y_vec <- y_prev * (1 - y_prev)
  y_mat <- matrix(y_vec, nrow = n, ncol = S)  # broadcast across columns
  
  eta_mat <- beta + delta * y_mat + gamma * dispersal
  mu_mat  <- inv_logit(eta_mat)
  
  return(pmin(pmax(mu_mat, 1e-6), 1 - 1e-6))  # Clip for stability
}

# E-step ------------------------------------------------------------------
e_step <- function(par, y_current, y_prev, wind_matrix, dist_matrix, d0 = 0.01, group_id, prior){
  
  mu_mat <- get_mu(par, par, y_current, y_prev, wind_matrix, dist_matrix, d0 = 0.01, group_id)
  
  S <- length(unique(group_id))
  
  #If the prior is a constant, make it a vector of length n
  if(length(prior) == 1){
    prior <- rep(prior, S)
  }
  
  #Compute the weighted likelihoods
  wl_mat <- prior * lik_beta(y_current, mu_mat, par[["phi"]], sum = F, log = F)
  
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
  
  phi <- pmax(phi, 1e-6)
  
  # Compute predicted mu matrix (n x S)
  mu_mat <- get_mu(par, y_current, y_prev, wind_matrix, dist_matrix, d0, group_id)
  
  # Compute derivatives
  y_star <- logit(y_current)
  mu_star <- digamma(mu_mat * phi) - digamma((1 - mu_mat) * phi)
  weight <- phi * (y_star - mu_star) * mu_mat * (1 - mu_mat)
  
  
  # Create model matrices for each term (n x S)
  y_prev_term <- outer(y_prev * (1 - y_prev), rep(1, ncol(mu_mat)))
  dispersal_term <- kappa_inner_sum(y_prev, wind_matrix, dist_matrix, d0, kappa, group_id = group_id)
  dispersal_term_deriv <- kappa_inner_sum(y_prev, wind_matrix, dist_matrix, d0, kappa, derivative = TRUE, group_id = group_id)
  
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
