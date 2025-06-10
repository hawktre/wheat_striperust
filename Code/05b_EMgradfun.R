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
lik_beta <- function(y, mu_mat, phi, sum = TRUE, log = TRUE, neg = F) {
  # Expand y and phi to match shape of mu if needed
  y <- matrix(y, nrow = nrow(mu_mat), ncol = ncol(mu_mat))  # Repeat y across columns
  phi <- matrix(phi, nrow = nrow(mu_mat), ncol = ncol(mu_mat))  # Repeat phi across columns (if scalar)
  
  # Compute log-likelihood elementwise
  val <- lgamma(phi) -
    lgamma(mu_mat * phi) -
    lgamma((1 - mu_mat) * phi) +
    (mu_mat * phi - 1) * log(y) +
    ((1 - mu_mat) * phi - 1) * log(1 - y)
  
  # Convert to raw likelihood if needed
  if (!log) {
    val <- exp(val)
  }
  
  # Create variable to return negative if requested
  if(neg){
    neg_val <- -1
  }else{
    neg_val <- 1
  }
  
  # Return row sums or full matrix
  if (sum) {
    return(neg_val * rowSums(val))
  } else {
    return(val)
  }
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
e_step <- function(y, mu_mat, phi, prior, group_id){
  #If the prior is a constant, make it a vector of length n
  if(length(prior) == 1){
    prior <- rep(prior, n)
  }
  
  #Compute the weighted likelihoods
  wl_mat <- prior * lik_beta(y, mu_mat, phi)
  
  #Compute posterior probabilities
  p_mat <- wl_mat / rowSums(wl_mat)
  
  return(p_mat)
}


# M-step ------------------------------------------------------------------
mstep_grad <- function(par, y_current, y_prev, wind_list, dist_list, p_hat_mat, d0 = 0.01) {
  beta  <- par["beta"]
  delta <- par["delta"]
  gamma <- par["gamma"]
  kappa <- par["kappa"]
  phi   <- pmax(par["phi"], 1e-6)
  
  n <- length(y_current)
  S <- length(wind_list)
  
  # Initialize gradients
  d_beta  <- 0
  d_delta <- 0
  d_gamma <- 0
  d_kappa <- 0
  d_phi   <- 0
  
  y_star <- logit(y_current)
  
  for (s in 1:S) {
    wind_matrix <- wind_list[[s]]
    dist_matrix <- dist_list[[s]]
    p_hat <- p_hat_mat[, s]
    
    # Get eta and mu for group s
    dispersal <- kappa_inner_sum(y_prev, wind_matrix, dist_matrix, d0, kappa)
    eta <- beta + delta * (y_prev * (1 - y_prev)) + gamma * dispersal
    mu  <- inv_logit(eta)
    mu  <- pmin(pmax(mu, 1e-6), 1 - 1e-6)
    
    mu_star <- digamma(mu * phi) - digamma((1 - mu) * phi)
    grad <- phi * (y_star - mu_star) * mu * (1 - mu)
    
    # Apply group-weighted contributions
    d_beta  <- d_beta  + sum(p_hat * grad)
    d_delta <- d_delta + sum(p_hat * grad * y_prev * (1 - y_prev))
    d_gamma <- d_gamma + sum(p_hat * grad * dispersal)
    
    dispersal_deriv <- kappa_inner_sum(y_prev, wind_matrix, dist_matrix, d0, kappa, derivative = TRUE)
    d_kappa <- d_kappa + sum(p_hat * grad * (-gamma) * dispersal_deriv)
    
    d_phi <- d_phi + sum(p_hat * (
      mu * (y_star - mu_star) +
        log(1 - y_current) - digamma((1 - mu) * phi) +
        digamma(phi)
    ))
  }
  
  # Return negative gradient vector
  -c(beta = d_beta, delta = d_delta, gamma = d_gamma, kappa = d_kappa, phi = d_phi)
}

