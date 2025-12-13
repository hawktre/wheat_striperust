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
kappa_inner_sum_backward <- function(par, y_prev, wind_mat, dist_mat, d0 = 0.01, derivative = FALSE, group_id, component) {

  # Compute constant terms
  dist_shifted <- dist_mat + d0
  log_dist <- log(dist_shifted)
  kappa <- par[['kappa']]

  #Assume infection came from group k
  group_k <- which(group_id %in% component)
  y_prev_k <- y_prev[group_k]

  #Compute distance kernel
  dist_kernel_k <- dist_shifted^(-kappa)
  if(derivative){
    dist_kernel_k <- -dist_kernel_k * log_dist
  }

  #Compute spread matrix
  wind_dist <- wind_mat * dist_kernel_k
  dispersal_k <- wind_dist[,group_k, drop = F] %*% y_prev_k

  return(dispersal_k)  # n x K matrix of dispersal from each group
}

# Mean function -----------------------------------------------------------
get_mu <- function(par, y_prev, wind_mat, dist_mat, group_id, component, d0 = 0.01) {
  # recompute dispersal at current kappa
  dispersal <- kappa_inner_sum_backward(
    par = par,
    y_prev = y_prev,
    wind_mat = wind_mat,
    dist_mat = dist_mat,
    d0 = d0,
    derivative = FALSE,
    group_id = group_id,
    component = component
  )

  # autoinfection term
  auto <- y_prev * (1 - y_prev)

  # linear predictor
  eta <- par["beta"] + par["delta"] * auto + par["gamma"] * dispersal

  # inverse logit, clipped for numerical stability
  mu <- 1 / (1 + exp(-eta))
  pmin(pmax(mu, 1e-6), 1 - 1e-6)
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
  K <- ncol(mu_mat)
  phi <- exp(phi)
  
  lik_list <-   lapply(seq_len(K), function(k) {
      loglik_zibeta(y, mu = mu_mat[, k], phi = phi[k], sum = F, log = T)
    })
  
  lik_mat <- do.call(cbind, lik_list)

  inner_sum <- t(t(lik_mat) + log(pi_vec))

  #compute Q-val
  Q_val <- sum(p_mat * inner_sum)

  return(-Q_val)
}

loglik_obs <- function(y, mu_mat, phi, pi_vec) {
  K <- ncol(mu_mat)
  phi <- exp(phi)
  
  # Compute log-likelihood contributions for each component
  lik_list <- lapply(seq_len(K), function(k) {
    loglik_zibeta(y, mu = mu_mat[, k], phi = phi[k], sum = FALSE, log = FALSE)
  })
  
  # Combine into n × K matrix of component densities (not logs)
  lik_mat <- do.call(cbind, lik_list)
  
  # Mixture density for each observation: sum_k pi_k * f_k(y_i)
  mix_density <- lik_mat %*% pi_vec
  
  # Avoid log(0)
  mix_density <- pmax(mix_density, .Machine$double.eps)
  
  # Observed-data log-likelihood
  loglik_val <- sum(log(mix_density))
  
  return(loglik_val)
}



# E-step ------------------------------------------------------------------
e_step <- function(y, mu_mat, phi, prior){
  phi <- exp(phi)
  K <- ncol(mu_mat)

  lik_list <-   lapply(seq_len(K), function(k) {
      loglik_zibeta(y, mu = mu_mat[, k], phi = phi[k], sum = F, log = F)
    })
  
  lik_mat <- do.call(cbind, lik_list)

  # Multiply each row (observation) element-wise by prior[k]
  wl_mat <- t(t(lik_mat) * prior)

  # Posterior probabilities for each observation and component
  p_mat <- wl_mat / rowSums(wl_mat)
  
  return(p_mat)
}

# Function to compute Q for M-step optimization
m_step_obj <- function(par, y_current, y_prev, wind_mat, dist_mat, group_id, component, p_vec) {
  # extract phi (still on log scale)
  phi <- exp(par[["phi"]])

  # recompute mu at current parameters
  mu_vec <- get_mu(
    par       = par,
    y_prev    = y_prev,
    wind_mat  = wind_mat,
    dist_mat  = dist_mat,
    group_id  = group_id,
    component = component
  )

  # compute log-likelihood for each observation
  ll_vec <- loglik_zibeta(y_current, mu_vec, phi, sum = FALSE, log = TRUE)

  # expected complete-data log-likelihood (Q_k)
  Q_k <- sum(p_vec * ll_vec)

  return(-Q_k)
}


# M-step (Individual Parameters (for eack K)) ------------------------------------------------------------------
m_step_grad <- function(par, y_current, y_prev, wind_mat, dist_mat, group_id, component, p_vec) {
  beta  <- par["beta"]
  delta <- par["delta"]
  gamma <- par["gamma"]
  kappa <- par["kappa"]
  log_phi <- par["phi"]
  phi <- exp(log_phi)

  non_zero <- which(y_current > 0)

  # Recompute dispersal and its derivative (with respect to kappa)
  dispersal <- kappa_inner_sum_backward(
    par = par, 
    y_prev = y_prev,
    wind_mat = wind_mat, 
    dist_mat = dist_mat,
     derivative = FALSE,
    group_id = group_id, 
    component = component
  )

  dispersal_grad <- kappa_inner_sum_backward(
    par = par, 
    y_prev = y_prev,
    wind_mat = wind_mat, 
    dist_mat = dist_mat,
     derivative = T,
    group_id = group_id, 
    component = component
  )

  # Recompute mu
  mu_vec <- get_mu(
    par = par,
    y_prev = y_prev,
    wind_mat = wind_mat,
    dist_mat = dist_mat,
    group_id = group_id,
    component = component
  )

  # Subset to nonzero responses
  y_current <- y_current[non_zero]
  auto_vec <- y_prev[non_zero] * (1 - y_prev[non_zero])
  dispersal <- dispersal[non_zero]
  dispersal_grad <- dispersal_grad[non_zero]
  p_vec <- p_vec[non_zero]
  mu_vec <- mu_vec[non_zero]

  # Core derivatives
  y_star <- logit(y_current)
  mu_star <- digamma(mu_vec * phi) - digamma((1 - mu_vec) * phi)
  weight <- phi * (y_star - mu_star) * mu_vec * (1 - mu_vec)

  # Gradients
  d_beta  <- sum(p_vec * weight)
  d_delta <- sum(p_vec * weight * auto_vec)
  d_gamma <- sum(p_vec * weight * dispersal)
  d_kappa <- sum(p_vec * weight * (-gamma) * dispersal_grad)

  # φ gradient (log-scale)
  d_phi_raw <- sum(p_vec * (
    digamma(phi) -
      mu_vec * digamma(mu_vec * phi) -
      (1 - mu_vec) * digamma((1 - mu_vec) * phi) +
      mu_vec * log(y_current) +
      (1 - mu_vec) * log(1 - y_current)
  ))
  d_phi <- phi * d_phi_raw  # chain rule for log(φ)

  -c(beta = d_beta, delta = d_delta, gamma = d_gamma, kappa = d_kappa, phi = d_phi)
}

