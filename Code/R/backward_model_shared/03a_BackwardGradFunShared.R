## ---------------------------
##
## Script name: 03a_BackwardGradFunShared_Vectorized.R
##
## Purpose of script: Fit em mod to get source probabilities
##                    Updated for log-scale parameters (delta, gamma, kappa, phi)
##                    and zero boundary clipping removals.
##
## Author: Trent VanHawkins
##
## Date Created: 2025-06-06
##
## VECTORIZED VERSION - optimized for speed
##
## ---------------------------

# --- Logit and Inverse Logit ---
logit <- function(p) log(p / (1 - p))
inv_logit <- function(x) 1 / (1 + exp(-x))

# Component indicator matrix ----------------------------------------------
create_component_indicator <- function(group_id, components) {
  n <- length(group_id)
  K <- ncol(components)
  C <- matrix(0, nrow = n, ncol = K)

  for (k in 1:K) {
    C[group_id %in% components[, k], k] <- 1
  }
  return(C)
}

# Vectorized Dispersal function -------------------------------------------
kappa_inner_sum_backward_vectorized <- function(
  par,
  y_prev,
  wind_mat,
  dist_mat,
  component_indicator,
  d0 = 0.01,
  derivative = FALSE,
  return_both = FALSE
) {
  dist_shifted <- dist_mat + d0
  # Extract kappa from the log-scale parameter
  kappa <- exp(par[['log_kappa']])
  dist_kernel <- dist_shifted^(-kappa)

  # Wind-distance matrix
  wind_dist <- wind_mat * dist_kernel

  # Vectorized computation
  y_prev_weighted <- y_prev * component_indicator
  dispersal <- wind_dist %*% y_prev_weighted

  if (derivative || return_both) {
    log_dist <- log(dist_shifted)
    dist_kernel_grad <- -dist_kernel * log_dist
    wind_dist_grad <- wind_mat * dist_kernel_grad
    dispersal_grad <- wind_dist_grad %*% y_prev_weighted

    if (return_both) {
      return(list(value = dispersal, gradient = dispersal_grad))
    } else {
      return(dispersal_grad)
    }
  }

  return(dispersal)
}

# Vectorized Mean function ------------------------------------------------
get_mu_vectorized <- function(
  par,
  y_prev,
  wind_mat,
  dist_mat,
  component_indicator,
  d0 = 0.01
) {
  # par maps log-parameters
  beta  <- par["beta"]
  delta <- exp(par["log_delta"])
  gamma <- exp(par["log_gamma"])

  dispersal_all <- kappa_inner_sum_backward_vectorized(
    par = par,
    y_prev = y_prev,
    wind_mat = wind_mat,
    dist_mat = dist_mat,
    component_indicator = component_indicator,
    d0 = d0,
    derivative = FALSE
  )

  auto <- y_prev * (1 - y_prev)
  eta <- beta + delta * auto + gamma * dispersal_all
  
  # Return pure continuous probabilities without clipping thresholds
  inv_logit(eta)
}

# Compute Likelihood ------------------------------------------------------
loglik_zibeta_vectorized <- function(y, mu_mat, phi, sum = TRUE, log = TRUE) {
  n <- length(y)

  if (is.vector(mu_mat)) {
    mu_mat <- matrix(mu_mat, ncol = 1)
  }
  K <- ncol(mu_mat)
  alpha <- mean(y == 0)

  ll_mat <- matrix(0, nrow = n, ncol = K)
  zero_idx <- (y == 0)
  nonzero_idx <- !zero_idx

  ll_mat[zero_idx, ] <- log(alpha)

  if (any(nonzero_idx)) {
    y_nz <- y[nonzero_idx]
    mu_nz <- mu_mat[nonzero_idx, , drop = FALSE]

    a <- mu_nz * phi
    b <- (1 - mu_nz) * phi

    ll_mat[nonzero_idx, ] <- log(1 - alpha) +
      lgamma(phi) -
      lgamma(a) -
      lgamma(b) +
      (a - 1) * log(y_nz) + 
      (b - 1) * log(1 - y_nz)
  }

  if (sum) {
    out <- colSums(ll_mat)
    if (!log) out <- exp(out)
  } else {
    out <- if (log) ll_mat else exp(ll_mat)
  }

  return(out)
}

# E-step (VECTORIZED) -----------------------------------------------------
e_step <- function(
  par,
  y_current,
  y_prev,
  wind,
  dist,
  group_id,
  components,
  pi_vec,
  component_indicator = NULL
) {
  phi <- exp(par[['log_phi']])
  K <- ncol(components)

  if (is.null(component_indicator)) {
    component_indicator <- create_component_indicator(group_id, components)
  }

  mu_all <- get_mu_vectorized(
    par = par,
    y_prev = y_prev,
    wind_mat = wind,
    dist_mat = dist,
    component_indicator = component_indicator
  )

  lik_mat <- loglik_zibeta_vectorized(
    y_current,
    mu_mat = mu_all,
    phi = phi,
    sum = FALSE,
    log = FALSE
  )

  wl_mat <- t(t(lik_mat) * pi_vec)
  p_mat <- wl_mat / rowSums(wl_mat)

  mix_density <- lik_mat %*% pi_vec
  mix_density <- pmax(mix_density, .Machine$double.eps)
  ll_obs <- sum(log(mix_density))

  return(list("p_mat" = p_mat, "ll_obs" = -ll_obs))
}

# M-step (Shared Parameters) ----------------------------------------------
m_step <- function(
  theta_old,
  intensity,
  intensity_prev,
  wind,
  dist,
  group_id,
  p_mat,
  components,
  component_indicator = NULL,
  max_iter = 1000,
  tol = 1e-4
) {
  # Assumes incoming theta_old maps parameters matching the new log layout
  fit <- tryCatch(
    optim(
      par = theta_old,
      fn = m_step_obj,
      gr = m_step_grad,
      method = "BFGS",
      control = list(maxit = 1000, reltol = 1e-8),
      y_current = intensity,
      y_prev = intensity_prev,
      wind_mat = wind,
      dist_mat = dist,
      group_id = group_id,
      p_mat = p_mat,
      components = components,
      component_indicator = component_indicator
    ),
    error = function(e) {
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
m_step_obj <- function(
  par,
  y_current,
  y_prev,
  wind_mat,
  dist_mat,
  group_id,
  p_mat,
  components,
  component_indicator = NULL
) {
  K <- ncol(p_mat)
  phi <- exp(par[["log_phi"]])

  if (is.null(component_indicator)) {
    component_indicator <- create_component_indicator(group_id, components)
  }

  mu_all <- get_mu_vectorized(
    par = par,
    y_prev = y_prev,
    wind_mat = wind_mat,
    dist_mat = dist_mat,
    component_indicator = component_indicator
  )

  ll_mat <- loglik_zibeta_vectorized(
    y_current,
    mu_mat = mu_all,
    phi = phi,
    sum = FALSE,
    log = TRUE
  )
  ll <- sum(p_mat * ll_mat)

  return(-ll)
}

# M-step gradient (VECTORIZED) --------------------------------------------
m_step_grad <- function(
  par,
  y_current,
  y_prev,
  wind_mat,
  dist_mat,
  group_id,
  p_mat,
  components,
  component_indicator = NULL
) {
  # Pull to natural scale to preserve the structure of your original calculus
  beta  <- par["beta"]
  delta <- exp(par["log_delta"])
  gamma <- exp(par["log_gamma"])
  kappa <- exp(par["log_kappa"])
  phi   <- exp(par["log_phi"])

  K <- ncol(p_mat)
  non_zero <- which(y_current > 0)

  if (is.null(component_indicator)) {
    component_indicator <- create_component_indicator(group_id, components)
  }

  dispersal_combined <- kappa_inner_sum_backward_vectorized(
    par = par,
    y_prev = y_prev,
    wind_mat = wind_mat,
    dist_mat = dist_mat,
    component_indicator = component_indicator,
    return_both = TRUE
  )

  dispersal <- dispersal_combined$value
  dispersal_grad <- dispersal_combined$gradient

  mu_mat <- get_mu_vectorized(
    par = par,
    y_prev = y_prev,
    wind_mat = wind_mat,
    dist_mat = dist_mat,
    component_indicator = component_indicator
  )

  # Broadcast vectors to matrices
  y_current_mat <- matrix(y_current, nrow = length(y_current), ncol = K)
  y_prev_mat    <- matrix(y_prev, nrow = length(y_prev), ncol = K)

  # Subset to nonzero responses
  y_current_nz   <- y_current_mat[non_zero, , drop = FALSE]
  auto_mat       <- y_prev_mat[non_zero, , drop = FALSE] * (1 - y_prev_mat[non_zero, , drop = FALSE])
  dispersal      <- dispersal[non_zero, , drop = FALSE]
  dispersal_grad <- dispersal_grad[non_zero, , drop = FALSE]
  p_mat          <- p_mat[non_zero, , drop = FALSE]
  mu_mat         <- mu_mat[non_zero, , drop = FALSE]

  # Core derivatives
  y_star <- logit(y_current_nz)
  mu_star <- digamma(mu_mat * phi) - digamma((1 - mu_mat) * phi)
  weight <- phi * (y_star - mu_star) * mu_mat * (1 - mu_mat)

  # Natural scale gradients
  d_beta  <- sum(p_mat * weight)
  d_delta <- sum(p_mat * weight * auto_mat)
  d_gamma <- sum(p_mat * weight * dispersal)
  d_kappa <- sum(p_mat * weight * (-gamma) * dispersal_grad)
  
  d_phi_natural <- sum(p_mat *
    (digamma(phi) -
      mu_mat * digamma(mu_mat * phi) -
      (1 - mu_mat) * digamma((1 - mu_mat) * phi) +
      mu_mat * log(y_current_nz) +
      (1 - mu_mat) * log(1 - y_current_nz)))

  # Apply structural chain rule conversions to log-scale variables
  d_log_delta <- d_delta * delta
  d_log_gamma <- d_gamma * gamma
  d_log_kappa <- d_kappa * kappa
  d_log_phi   <- d_phi_natural * phi

  -c(
    beta      = d_beta,
    log_delta = d_log_delta,
    log_gamma = d_log_gamma,
    log_kappa = d_log_kappa,
    log_phi   = d_log_phi
  )
}