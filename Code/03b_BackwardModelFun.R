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
library(tidyverse)

## Read in the needed functions
source(here("Code/03a_BackwardGradFun.R"))
source(here("Code/02a_ForwardGradFun.R"))

# Backward Fit ------------------------------------------------------------
backward_fit <- function(config, blk, trt, vst, mod_dat, forward_fits, kappa = NULL) {

  # Setup
  intensity <- mod_dat$intensity[, blk, trt, vst]
  intensity_prev <- mod_dat$intensity[, blk, trt, as.numeric(vst) - 1]
  wind <- mod_dat$wind[, , blk, trt, vst]
  dist <- mod_dat$dist
  groups <- mod_dat$groups[, config]
  n_groups <- parse_number(config)
  prior <- rep(1 / n_groups, n_groups)
  
  # Starting values from forward fit (grid search or MLE approach)
  if(!is.null(kappa)){
    theta <- initialize_theta(y_cur = intensity, 
                                   y_prev = intensity_prev, 
                                   wind_mat = wind, 
                                   dist_mat = dist, 
                                   kappa_try = as.numeric(kappa),
                                   d_0 = 0.01)
  }else{
    theta <- forward_fits[block == blk & treat == as.numeric(trt) & visit == as.numeric(vst)
  ][["theta"]][[1]]
  }
  
  
  max_em_iter <- 1000
  tol <- 1e-6
  q_track <- numeric(max_em_iter)
  
  # E-step

  ## Compute mu_matrix
  mu_mat <- get_mu(par = theta, y_prev = intensity_prev, wind_matrix = wind, dist_matrix = dist, group_id = groups)

  ## Compute conditional probability matrix
  p_mat <- e_step(y_vec = intensity, mu_mat = mu_mat, phi = theta[["phi"]], prior_vec = prior)
  
  # EM loop
  for (em_iter in seq_len(max_em_iter)) {
    fit <- tryCatch(
      optim(
        par = theta,
        fn = wrapped_obj,
        gr = mstep_grad_em,
        method = "BFGS",
        control = list(maxit = 1000, reltol = tol),
        p_mat = p_mat,
        y_current = intensity,
        y_prev = intensity_prev,
        wind_matrix = wind,
        dist_matrix = dist,
        group_id = groups
      ),
      error = function(e) {
        message(sprintf(
          "EM step %d failed in backward_fit [config=%s, block=%s, treat=%s, visit=%s]: %s",
          em_iter, config, blk, trt, vst, conditionMessage(e)
        ))
        
        stop(e)
      }
    )
    
    if (is.null(fit) || fit$convergence != 0 || !is.finite(fit$value)) {
      return(data.table(
        config = config,
        block = blk,
        treat = as.numeric(trt),
        visit = as.numeric(vst),
        kappa = kappa,
        em_iters = em_iter,
        converged = FALSE,
        neg_loglik = NA_real_,
        final_theta = list(theta),
        p_mat = list(p_mat)
      ))
    }
    
    #Update parameters
    theta_new <- fit$par
    
    # Update E-step
    p_mat <- e_step(
      par = theta_new,
      y_current = intensity,
      y_prev = intensity_prev,
      wind_matrix = wind,
      dist_matrix = dist,
      group_id = groups,
      prior = prior
    )
    
    # Compute Q
    mu_mat <- get_mu(
      par = theta_new,
      y_prev = intensity_prev,
      wind_matrix = wind,
      dist_matrix = dist,
      group_id = groups
    )
    
    q_val <- Q_fun(
      y = intensity,
      mu_mat = mu_mat,
      phi = theta_new[["phi"]],
      p_mat = p_mat
    )
    
    if (!is.finite(q_val)) {
      return(data.table(
        config = config,
        block = blk,
        treat = as.numeric(trt),
        visit = as.numeric(vst),
        kappa = kappa,
        em_iters = em_iter,
        converged = FALSE,
        neg_loglik = NA_real_,
        final_theta = list(theta),
        p_mat = list(p_mat)
      ))
    }
    q_track[em_iter] <- q_val
    
    
    if (em_iter > 1 && abs(q_track[em_iter] - q_track[em_iter - 1]) < tol) {
      break
    }
    
    theta <- theta_new
  }
  
  # Return a single-row data.table
  result <- data.table(
    config = config,
    block = blk,
    treat = trt,
    visit = vst,
    kappa = kappa,
    em_iters = em_iter,
    converged = em_iter < max_em_iter,
    neg_loglik = q_track[em_iter],
    final_theta = list(theta_new),
    p_mat = list(p_mat)
  )
  
  return(result)
}

# Backward_Preds ----------------------------------------------------------
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
    predicted_source = predicted_source,
    true_source = true_source,
    dist_error = error,
    dist_acc = dist_acc
  )
  # 5. Compute other metrics 
  return(as.data.table(result))
}