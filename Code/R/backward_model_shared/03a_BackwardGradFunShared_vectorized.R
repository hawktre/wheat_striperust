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
## VECTORIZED VERSION - optimized for speed
##
## ---------------------------
## ---------------------------

# --- Logit and Inverse Logit ---
logit <- function(p) log(p / (1 - p))
inv_logit <- function(x) 1 / (1 + exp(-x))

# Component indicator matrix ----------------------------------------------
create_component_indicator <- function(group_id, components) {
  
  # Creates n x K indicator matrix where C[i,k] = 1 if obs i belongs to component k
  n <- length(group_id)
  K <- ncol(components)
  C <- matrix(0, nrow = n, ncol = K)
  
  for (k in 1:K) {
    C[group_id %in% components[, k], k] <- 1
  }
  return(C)
}

# Vectorized Dispersal function -------------------------------------------
kappa_inner_sum_backward_vectorized <- function(par, y_prev, wind_mat, dist_mat, 
                                                 component_indicator, d0 = 0.01, 
                                                 derivative = FALSE) {
  # Compute distance kernel (same for all components)
  dist_shifted <- dist_mat + d0
  kappa <- par[['kappa']]
  dist_kernel <- dist_shifted^(-kappa)
  
  if (derivative) {
    dist_kernel <- -dist_kernel * log(dist_shifted)
  }
  
  # Compute wind-distance matrix once
  wind_dist <- wind_mat * dist_kernel  # n x n
  
  # Vectorized: weight y_prev by component membership
  # y_prev_weighted is n x K
  y_prev_weighted <- y_prev * component_indicator
  
  # Single matrix multiplication gives all K dispersals at once
  dispersal_all <- wind_dist %*% y_prev_weighted  # n x K
  
  return(dispersal_all)
}

# Vectorized Mean function ------------------------------------------------
get_mu_vectorized <- function(par, y_prev, wind_mat, dist_mat, component_indicator, d0 = 0.01) {
  # Compute dispersal for all components at once (n x K matrix)
  dispersal_all <- kappa_inner_sum_backward_vectorized(
    par = par,
    y_prev = y_prev,
    wind_mat = wind_mat,
    dist_mat = dist_mat,
    component_indicator = component_indicator,
    d0 = d0,
    derivative = FALSE
  )
  
  # Auto-infection term (broadcast to n x K)
  auto <- y_prev * (1 - y_prev)
  
  # Linear predictor for all components (n x K)
  eta <- par["beta"] + par["delta"] * auto + par["gamma"] * dispersal_all
  
  # Inverse logit, clipped
  mu_all <- 1 / (1 + exp(-eta))
  pmin(pmax(mu_all, 1e-6), 1 - 1e-6)
}

# Compute Likelihood ------------------------------------------------------
# Vectorized likelihood that accepts mu as a matrix
loglik_zibeta_vectorized <- function(y, mu_mat, phi, sum = TRUE, log = TRUE) {
  # y: vector of length n
  # mu_mat: n x K matrix (or vector for single component)
  # phi: scalar
  
  n <- length(y)
  
  # Handle both vector and matrix input
  if (is.vector(mu_mat)) {
    mu_mat <- matrix(mu_mat, ncol = 1)
  }
  K <- ncol(mu_mat)
  
  # Compute alpha once
  alpha <- mean(y == 0)
  
  # Initialize log-likelihood matrix (n x K)
  ll_mat <- matrix(0, nrow = n, ncol = K)
  
  # Identify zero observations
  zero_idx <- (y == 0)
  nonzero_idx <- !zero_idx
  
  # Zero observations - broadcast across all K components
  ll_mat[zero_idx, ] <- log(alpha)
  
  # Non-zero observations (vectorized across both n and K)
  if (any(nonzero_idx)) {
    y_nz <- y[nonzero_idx]
    mu_nz <- mu_mat[nonzero_idx, , drop = FALSE]  # subset to nonzero rows, keep all K columns
    
    # Vectorized beta parameters (works element-wise on matrix)
    a <- mu_nz * phi
    b <- (1 - mu_nz) * phi
    
    # Vectorized log-likelihood computation
    ll_mat[nonzero_idx, ] <- log(1 - alpha) +
      lgamma(phi) - lgamma(a) - lgamma(b) +
      (a - 1) * log(y_nz) +  # y_nz broadcasts across columns
      (b - 1) * log(1 - y_nz)
  }
  
  if (sum) {
    out <- colSums(ll_mat)  # sum over observations, return K-length vector
    if (!log) out <- exp(out)
  } else {
    out <- if (log) ll_mat else exp(ll_mat)
  }
  
  return(out)
}

# E-step (VECTORIZED) -----------------------------------------------------
e_step <- function(par, y_current, y_prev, wind, dist, group_id, components, pi_vec){
  
  # Set up variables
  phi <- exp(par[['phi']])
  K <- ncol(components)
  
  # Create component indicator ONCE
  component_indicator <- create_component_indicator(group_id, components)
  
  # Get mu for ALL components at once (n x K matrix)
  mu_all <- get_mu_vectorized(
    par = par,
    y_prev = y_prev,
    wind_mat = wind,
    dist_mat = dist,
    component_indicator = component_indicator
  )
  
  # Vectorized likelihood computation for all components
  lik_mat <- loglik_zibeta_vectorized(y_current, mu_mat = mu_all, phi = phi, sum = FALSE, log = FALSE)
  
  # Multiply each row (observation) element-wise by prior[k] and get weighted posterior probabilities
  wl_mat <- t(t(lik_mat) * pi_vec)
  p_mat <- wl_mat / rowSums(wl_mat)
  
  # Compute observed-data ll
  mix_density <- lik_mat %*% pi_vec
  mix_density <- pmax(mix_density, .Machine$double.eps)
  ll_obs <- sum(log(mix_density))
  
  return(list("p_mat" = p_mat, "ll_obs" = -ll_obs))
}

# M-step (Shared Parameters) ----------------------------------------------
m_step <- function(theta_old, intensity, intensity_prev, wind, dist, group_id, p_mat, components, max_iter = 1000, tol = 1e-4){
  
  fit <- tryCatch(optim(
    par     = theta_old,
    fn      = m_step_obj,
    gr      = m_step_grad,
    method  = "BFGS",
    control = list(maxit = 1000, reltol = 1e-8),
    y_current = intensity,
    y_prev    = intensity_prev,
    wind_mat  = wind,
    dist_mat  = dist,
    group_id  = group_id,
    p_mat = p_mat,
    components = components
  ),
  error = function(e) {
    message(sprintf("optim() error: %s", conditionMessage(e)))
    NULL
  }
  )
  if (is.null(fit) || fit$convergence != 0) {
    theta_new <- theta_old
  } else {
    theta_new <- fit$par
  }
  
  ## Update mixture weights 
  pi <- colSums(p_mat) / nrow(p_mat)
  
  return(list("theta_new" = theta_new, "pi" = pi))
}

# Function to compute Q for M-step optimization (VECTORIZED) -------------
m_step_obj <- function(par, y_current, y_prev, wind_mat, dist_mat, group_id, p_mat, components) {
  
  # Specify number of components
  K <- ncol(p_mat)
  # extract phi (still on log scale)
  phi <- exp(par[["phi"]])
  
  # Create component indicator ONCE
  component_indicator <- create_component_indicator(group_id, components)
  
  # Recompute mu at current parameters - ALL components at once
  mu_all <- get_mu_vectorized(
    par = par,
    y_prev = y_prev,
    wind_mat = wind_mat,
    dist_mat = dist_mat,
    component_indicator = component_indicator
  )
  
  # Compute log-likelihood for each component
  ll_vec <- loglik_zibeta_vectorized(y_current, mu_mat = mu_all, phi = phi, sum = T, log = T)
  
  # Expected complete-data log-likelihood (Q_k)
  ll <- sum(ll_vec)
  
  return(-ll)
}

# M-step gradient (VECTORIZED) --------------------------------------------
m_step_grad <- function(par, y_current, y_prev, wind_mat, dist_mat, group_id, p_mat, components) {
  beta  <- par["beta"]
  delta <- par["delta"]
  gamma <- par["gamma"]
  kappa <- par["kappa"]
  log_phi <- par["phi"]
  phi <- exp(log_phi)
  
  K <- ncol(p_mat)
  non_zero <- which(y_current > 0)
  
  # Create component indicator ONCE
  component_indicator <- create_component_indicator(group_id, components)
  
  # Compute ALL dispersals at once (n x K)
  dispersal <- kappa_inner_sum_backward_vectorized(
    par = par, 
    y_prev = y_prev,
    wind_mat = wind_mat, 
    dist_mat = dist_mat,
    component_indicator = component_indicator,
    derivative = FALSE
  )
  
  # Compute ALL dispersal gradients at once (n x K)
  dispersal_grad <- kappa_inner_sum_backward_vectorized(
    par = par, 
    y_prev = y_prev,
    wind_mat = wind_mat, 
    dist_mat = dist_mat,
    component_indicator = component_indicator,
    derivative = TRUE
  )
  
  # Compute ALL mus at once (n x K)
  mu_mat <- get_mu_vectorized(
    par = par,
    y_prev = y_prev,
    wind_mat = wind_mat,
    dist_mat = dist_mat,
    component_indicator = component_indicator
  )
  
  # Broadcast Y_current and Y_prev to matrices
  y_current_mat <- matrix(y_current, nrow = length(y_current), ncol = K)
  y_prev_mat <- matrix(y_prev, nrow = length(y_prev), ncol = K)
  
  # Subset to nonzero responses
  y_current_nz <- y_current_mat[non_zero, ]
  auto_mat <- y_prev_mat[non_zero, ] * (1 - y_prev_mat[non_zero, ])
  dispersal <- dispersal[non_zero, ]
  dispersal_grad <- dispersal_grad[non_zero, ]
  p_mat <- p_mat[non_zero, ]
  mu_mat <- mu_mat[non_zero, ]
  
  # Core derivatives
  y_star <- logit(y_current_nz)
  mu_star <- digamma(mu_mat * phi) - digamma((1 - mu_mat) * phi)
  weight <- phi * (y_star - mu_star) * mu_mat * (1 - mu_mat)
  
  # Gradients
  d_beta  <- sum(p_mat * weight)
  d_delta <- sum(p_mat * weight * auto_mat)
  d_gamma <- sum(p_mat * weight * dispersal)
  d_kappa <- sum(p_mat * weight * (-gamma) * dispersal_grad)
  
  # φ gradient (log-scale)
  d_phi_raw <- p_mat * (
    digamma(phi) -
      mu_mat * digamma(mu_mat * phi) -
      (1 - mu_mat) * digamma((1 - mu_mat) * phi) +
      mu_mat * log(y_current_nz) +
      (1 - mu_mat) * log(1 - y_current_nz)
  )
  
  d_phi <- sum(phi * d_phi_raw)  # chain rule for log(φ)
  
  -c(beta = d_beta, delta = d_delta, gamma = d_gamma, kappa = d_kappa, phi = d_phi)
}