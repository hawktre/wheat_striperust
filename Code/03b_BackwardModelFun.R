## ---------------------------
##
## Script name: 05c_RunEM.R
##
## Purpose of script: Run the EM algorithm using functions from 05b
##
## Author: Trent VanHawkins
##
## Date Created: 2025-06-11
##
##
## ---------------------------

## ---------------------------

## load up the packages we will need:  (uncomment as required)
library(here)
library(parallel)
library(data.table)

## Read in the needed functions
source(here("Code/03a_BackwardGradFun.R"))
source(here("Code/02a_ForwardGradFun.R"))

# Backward Fit ------------------------------------------------------------
backward_fit <- function(config, blk, trt, vst, mod_dat, inits, max_iter = 100, tol = 1e-2) {
  
  ## Extract needed data
  intensity <- mod_dat$intensity[,blk,trt,vst]
  intensity_prev <- mod_dat$intensity[,blk,trt,as.numeric(vst)-1]
  wind <- mod_dat$wind[,,blk,trt,vst]
  dist <- mod_dat$dist
  group_id <- mod_dat$groups[, config]
  
  # Initializations
  K <- length(unique(group_id)) #Number of mixture components
  pi <- rep(1/K, K) #Uniform prior
  theta_init <- inits
  theta_mat <- matrix(theta_init, nrow = length(theta_init), ncol = K) #Create matrix of parameters for each component
  rownames(theta_mat) <- names(theta_init)
  
  #Initialize mu_mat
  # Combine results into an n × K matrix
  mu_list <- lapply(seq_len(K), function(k) {
  get_mu(
    par        = theta_mat[, k],     # parameters for component k
    y_prev     = intensity_prev,     # previous intensities
    wind_mat   = wind,               # full wind matrix
    dist_mat   = dist,               # full distance matrix
    group_id   = group_id,           # vector of group assignments
    component  = k                 # current component index
  )
})

  # Combine results into an n × K matrix
  mu_mat <- do.call(cbind, mu_list)
  ## Set up convergence criteria and algorithm tracking
  max_em_iter <- max_iter
  q_track <- numeric(max_em_iter + 1)
  numCores <- detectCores() - 1
  # Initial E-step
  p_mat <- e_step(y = intensity, mu_mat = mu_mat, phi = theta_mat['phi',], prior = pi)

  ## Compute Q
  q_track[1] <- Q_fun(y = intensity, mu_mat = mu_mat, phi = theta_mat['phi',], pi_vec = pi, p_mat = p_mat)
  # EM loop
  for (iter in 1:max_em_iter) {
    # M-step (optimize theta for each component k) ---------------------------------------------------------------
    
    theta_mat <- do.call(cbind, mclapply(seq_len(K), function(k) {
      fit <- tryCatch(
        optim(
          par     = theta_mat[, k],
          fn      = m_step_obj,
          gr      = m_step_grad,
          method  = "BFGS",
          control = list(maxit = 1000, reltol = 1e-4),
          y_current = intensity,
          y_prev    = intensity_prev,
          wind_mat  = wind,
          dist_mat  = dist,
          group_id  = group_id,
          component = k,
          p_vec     = p_mat[, k]
        ),
        error = function(e) {
          message(sprintf("optim() error for component %s: %s", k, conditionMessage(e)))
          NULL
        }
      )
      if (is.null(fit) || fit$convergence != 0) {
        theta_mat[, k]
      } else {
        fit$par
      }
    }, mc.cores = numCores, mc.preschedule = FALSE))
    
    ## Update mixture weights 
    pi <- colSums(p_mat)/nrow(p_mat)

    # (E-step will go here next)

    ## Update mu_mat
    mu_mat <- do.call(cbind, lapply(seq_len(K), function(k) {
      get_mu(
        par        = theta_mat[, k],     # parameters for component k
        y_prev     = intensity_prev,     # previous intensities
        wind_mat   = wind,               # full wind matrix
        dist_mat   = dist,               # full distance matrix
        group_id   = group_id,           # vector of group assignments
        component  = k                  # current component index
      )}))
    
    p_mat <- e_step(y = intensity, mu_mat = mu_mat, phi = theta_mat['phi',], prior = pi)
    
    ## Compute Q
    q_track[iter + 1] <- Q_fun(y = intensity, mu_mat = mu_mat, phi = theta_mat['phi',], pi_vec = pi, p_mat = p_mat)

    if (!is.finite(q_track[iter + 1]) || (iter + 1) == max_em_iter) {
      return(data.table(
        config = config,
        block = blk,
        treat = trt,
        visit = as.numeric(vst),
        em_iters = iter + 1,
        converged = FALSE,
        Q_track = list(q_track[1:iter + 1]),
        Q_final = q_track[iter+1],
        theta = list(theta_mat),
        p_mat = list(p_mat),
        pi = list(pi)
      ))
    }
    
    if ((iter > 1 && abs(q_track[iter + 1] - q_track[iter]) < tol * (1 + abs(q_track[iter])))) {
      return(data.table(
        config = config,
        block = blk,
        treat = trt,
        visit = as.numeric(vst),
        em_iters = iter + 1,
        converged = T,
        Q_track = list(q_track[1:iter +1]),
        Q_final = q_track[iter+1],
        theta = list(theta_mat),
        p_mat = list(p_mat),
        pi = list(pi)))
    }


  }
}

source_pred <- function(config, blk, trt, vst, p_mat, mod_dat) {
  
  # Ensure groups is a factor with levels = 1:n_groups
  groups <- as.factor(mod_dat$groups[,config])
  levels(groups) <- sort(unique(groups))
  
  # 1. Compute p̄ matrix (average p_mat rows by destination group)
  group_ids <- sort(unique(groups))
  n_groups <- length(group_ids)
  
  p_bar <- matrix(NA, n_groups, n_groups)
  rownames(p_bar) <- colnames(p_bar) <- as.character(group_ids)
  
  for (g in group_ids) {
    rows_in_group <- which(groups == g)
    p_bar[as.character(g), ] <- colMeans(p_mat[rows_in_group, , drop = FALSE])
  }
  
  # 2. Predict source group for each destination group (row-wise argmax)
  predicted_source <- which(p_bar == max(p_bar), arr.ind = TRUE)[[2]]
  
  # 3. Compare to ground truth
  true_source <- mod_dat$truth[blk, config]  # vector of length n_groups
  
  # 4. Compute distance-weighted accuracy
  dist <- mod_dat$grid_dist[[config]]
  error <- dist[predicted_source, true_source]
  dist_acc <- 1 - (error / max(dist[predicted_source,]))  # normalized distance accuracy
  
  result <- list(
    config = config,
    block = blk,
    treat = trt,
    visit = vst,
    p_bar = list(p_bar),
    predicted_source = predicted_source,
    true_source = true_source,
    dist_error = error,
    dist_acc = dist_acc
  )
  # 5. Compute other metrics 
  return(as.data.table(result))
}

