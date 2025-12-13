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
      loglik_zibeta(y, mu = mu_mat[, k], phi = phi, sum = F, log = T)
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
    loglik_zibeta(y, mu = mu_mat[, k], phi = phi, sum = FALSE, log = FALSE)
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
e_step <- function(par, y_current, y_prev, wind, dist, group_id, components, pi_vec){
  #Set up variables
  phi <- exp(par[['phi']])
  K <- ncol(components)

  #Compute mu_mat using the current parameters
  mu_mat <- do.call(cbind, lapply(seq_len(K), function(k) get_mu(
    par = par,
    y_prev = y_prev,
    wind_mat = wind,
    dist_mat = dist,
    group_id = group_id,
    component = components[,k]
  )))

  #Compute the likelihood of each data point
  lik_list <- lapply(seq_len(K), function(k) {
      loglik_zibeta(y_current, mu = mu_mat[, k], phi = phi, sum = F, log = F)
    })
  
  lik_mat <- do.call(cbind, lik_list)

  # Multiply each row (observation) element-wise by prior[k]
  wl_mat <- t(t(lik_mat) * pi_vec)

  # Posterior probabilities for each observation and component
  p_mat <- wl_mat / rowSums(wl_mat)
  
  return(p_mat)
}

# M-step (Shared Parameters) ------------------------------------------------------------------
m_step <- function(theta_old, intensity, intensity_prev, wind, dist, group_id, p_mat, components, max_iter = 1000, tol = 1e-4){
  fit <- tryCatch(optim(
          par     = theta_old,
          fn      = m_step_obj,
          gr      = m_step_grad,
          method  = "BFGS",
          control = list(maxit = 1000, reltol = 1e-4),
          y_current = intensity,
          y_prev    = intensity_prev,
          wind_mat  = wind,
          dist_mat  = dist,
          group_id  = group_id,
          p_mat = p_mat,
          components = components
        ),
        error = function(e) {
          message(sprintf("optim() error for component %s: %s", k, conditionMessage(e)))
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
    
  return(list("theta_new" = theta_new,
  "pi" = pi))
}

# Function to compute Q for M-step optimization
m_step_obj <- function(par, y_current, y_prev, wind_mat, dist_mat, group_id, p_mat, components) {
  #Specify number of components
  K <- ncol(p_mat)
  # extract phi (still on log scale)
  phi <- exp(par[["phi"]])

  # recompute mu at current parameters
  mu_mat <- do.call(cbind, lapply(seq_len(K), function(k) get_mu(
    par = par,
    y_prev = y_prev,
    wind_mat = wind_mat,
    dist_mat = dist_mat,
    group_id = group_id,
    component = components[,k]
  )))

  # compute log-likelihood for each observation
  ll_mat <- do.call(cbind, lapply(seq_len(K), function(k) loglik_zibeta(y_current, 
    mu_mat[,k], 
    phi, 
    sum = FALSE, 
    log = TRUE)))

  # expected complete-data log-likelihood (Q_k)
  Q_k <- sum(p_mat * ll_mat)

  return(-Q_k)
}



m_step_grad <- function(par, y_current, y_prev, wind_mat, dist_mat, group_id, p_mat, components) {
  beta  <- par["beta"]
  delta <- par["delta"]
  gamma <- par["gamma"]
  kappa <- par["kappa"]
  log_phi <- par["phi"]
  phi <- exp(log_phi)


  K <- ncol(p_mat)
  non_zero <- which(y_current > 0)

  # Recompute dispersal and its derivative (with respect to kappa) (n x K)
  dispersal <- do.call(cbind, lapply(seq_len(K), function(k) kappa_inner_sum_backward(
    par = par, 
    y_prev = y_prev,
    wind_mat = wind_mat, 
    dist_mat = dist_mat,
     derivative = FALSE,
    group_id = group_id, 
    component = components[,k]
  ))) 

  dispersal_grad <- do.call(cbind, lapply(seq_len(K), function(k) kappa_inner_sum_backward(
    par = par, 
    y_prev = y_prev,
    wind_mat = wind_mat, 
    dist_mat = dist_mat,
     derivative = TRUE,
    group_id = group_id, 
    component = components[,k]
  ))) 

  # Recompute mu (n x K)
  mu_mat <- do.call(cbind, lapply(seq_len(K), function(k) get_mu(
    par = par,
    y_prev = y_prev,
    wind_mat = wind_mat,
    dist_mat = dist_mat,
    group_id = group_id,
    component = components[,k]
  )))

  #Broadcast Y_current and Y_prev to matrices
  y_current <- matrix(y_current, nrow = length(y_current), ncol = K)
  y_prev <- matrix(y_prev, nrow = length(y_prev), ncol = K)

  # Subset to nonzero responses
  y_current <- y_current[non_zero, ]
  auto_mat <- y_prev[non_zero, ] * (1 - y_prev[non_zero, ])
  dispersal <- dispersal[non_zero,]
  dispersal_grad <- dispersal_grad[non_zero,]
  p_mat <- p_mat[non_zero,]
  mu_mat <- mu_mat[non_zero,]

  # Core derivatives
  y_star <- logit(y_current)
  mu_star <- digamma(mu_mat * phi) - digamma((1 - mu_mat) * phi)
  weight <- phi * (y_star - mu_star) * mu_mat * (1 - mu_mat)

  # Gradients
  d_beta  <- sum(p_mat * weight)
  d_delta <- sum(p_mat * weight * auto_vec)
  d_gamma <- sum(p_mat * weight * dispersal)
  d_kappa <- sum(p_mat * weight * (-gamma) * dispersal_grad)

  # φ gradient (log-scale)
  d_phi_raw <- p_mat * (
    digamma(phi) -
      mu_mat * digamma(mu_mat * phi) -
      (1 - mu_mat) * digamma((1 - mu_mat) * phi) +
      mu_mat * log(y_current) +
      (1 - mu_mat) * log(1 - y_current)
  )
  
  d_phi <- sum(phi * d_phi_raw)  # chain rule for log(φ)

  -c(beta = d_beta, delta = d_delta, gamma = d_gamma, kappa = d_kappa, phi = d_phi)
}


