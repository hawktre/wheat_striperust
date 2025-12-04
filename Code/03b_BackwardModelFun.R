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
backward_fit <- function(config, blk, trt, vst, mod_dat, inits, max_iter = 100, tol = 1e-4) {

  ## Extract needed data
  intensity <- mod_dat$intensity[, blk, trt, vst]
  intensity_prev <- mod_dat$intensity[, blk, trt, as.numeric(vst) - 1]
  wind <- mod_dat$wind[, , blk, trt, vst]
  dist <- mod_dat$dist
  group_id <- mod_dat$groups[, config]
  
  # Initializations
  S <- as.numeric(trt)
  combos <- combn(sort(unique(group_id)), S)
  K <- ncol(combos)
  pi <- rep(1 / K, K)
  theta_init <- inits
  theta_old <- matrix(theta_init, nrow = length(theta_init), ncol = K)
  rownames(theta_old) <- names(theta_init)
  
  # Initialize mu_mat
  mu_list <- lapply(seq_len(K), function(k) {
    get_mu(
      par        = theta_old[, k],
      y_prev     = intensity_prev,
      wind_mat   = wind,
      dist_mat   = dist,
      group_id   = group_id,
      component  = combos[,k]
    )
  })
  mu_mat <- do.call(cbind, mu_list)
  
  ## Tracking
  max_em_iter <- max_iter
  q_track <- numeric(max_em_iter)
  observed_ll <- numeric(max_em_iter)
  numCores <- detectCores() - 1
  
  # Initial E-step
  p_mat <- e_step(y = intensity, mu_mat = mu_mat, phi = theta_old['phi',], prior = pi)
  
  ## Initial Q and log-likelihood
  q_track[1] <- Q_fun(y = intensity, mu_mat = mu_mat, phi = theta_old['phi',], pi_vec = pi, p_mat = p_mat)
  observed_ll[1] <- loglik_obs(y = intensity, mu_mat = mu_mat, phi = theta_old['phi',], pi_vec = pi)
  
  # EM loop
  for (iter in 2:max_em_iter) {
    
    # M-step
    theta_new <- do.call(cbind, mclapply(seq_len(K), function(k) {
      fit <- tryCatch(
        optim(
          par     = theta_old[, k],
          fn      = m_step_obj,
          gr      = m_step_grad,
          method  = "BFGS",
          control = list(maxit = 1000, reltol = 1e-4),
          y_current = intensity,
          y_prev    = intensity_prev,
          wind_mat  = wind,
          dist_mat  = dist,
          group_id  = group_id,
          component = combos[,k],
          p_vec     = p_mat[, k]
        ),
        error = function(e) {
          message(sprintf("optim() error for component %s: %s", k, conditionMessage(e)))
          NULL
        }
      )
      if (is.null(fit) || fit$convergence != 0) {
        theta_old[, k]
      } else {
        fit$par
      }
    }, mc.cores = numCores, mc.preschedule = FALSE))
    
    ## Update mixture weights 
    pi <- colSums(p_mat) / nrow(p_mat)
    
    ## Update mu_mat
    mu_mat <- do.call(cbind, lapply(seq_len(K), function(k) {
      get_mu(
        par        = theta_new[, k],
        y_prev     = intensity_prev,
        wind_mat   = wind,
        dist_mat   = dist,
        group_id   = group_id,
        component  = combos[,k]
      )
    }))
    
    ## E-step: recompute responsibilities
    p_mat <- e_step(y = intensity, mu_mat = mu_mat, phi = theta_new['phi',], prior = pi)
    
    ## Compute observed log-likelihood and Q
    observed_ll[iter] <- loglik_obs(y = intensity, mu_mat = mu_mat, phi = theta_new['phi',], pi_vec = pi)
    q_track[iter] <- Q_fun(y = intensity, mu_mat = mu_mat, phi = theta_new['phi',], pi_vec = pi, p_mat = p_mat)
    
    ## Convergence checks
    if (!is.finite(observed_ll[iter]) || iter == max_em_iter) {
      return(data.table(
        config = config, block = blk, treat = trt, visit = as.numeric(vst),
        em_iters = iter, converged = FALSE,
        Q_track = list(q_track[1:iter]), Q_final = q_track[iter],
        theta = list(theta_new), p_mat = list(p_mat), pi = list(pi)
      ))
    }
    
    theta_diff <- max(abs(theta_new - theta_old))
    observed_ll_diff <- abs(observed_ll[iter] - observed_ll[iter - 1])
    if (( observed_ll_diff < tol * (1 + abs(observed_ll[iter - 1]))) ||
        theta_diff < tol) {
      return(data.table(
        config = config, block = blk, treat = trt, visit = as.numeric(vst),
        em_iters = iter, converged = TRUE,
        Q_track = list(q_track[1:iter]), Q_final = q_track[iter],
        theta = list(theta_new), p_mat = list(p_mat), pi = list(pi)
      ))
    }
    
    theta_old <- theta_new
  }
}


source_pred <- function(config, blk, trt, vst, p_mat, mod_dat) {
  browser()
  # Ensure groups is a factor with levels = 1:n_groups
  groups <- as.factor(mod_dat$groups[,config])
  levels(groups) <- sort(unique(groups))
  
  # 1. Compute p̄ matrix (average p_mat rows by destination group)
  group_ids <- sort(unique(groups))
  n_groups <- length(group_ids)

  S <- as.numeric(trt)
  combos <- combn(sort(unique(group_ids)), S)
  K <- ncol(combos)
  
  p_bar <- matrix(NA, n_groups, K)
    
  for (g in 1:n_groups) {
    rows_in_group <- which(groups %in% group_ids[g])
    p_bar[g, ] <- colMeans(p_mat[rows_in_group, , drop = FALSE])
  }

  # 2. Predict source group for each destination group (row-wise argmax)
  predicted_source <- combos[,which(p_bar == max(p_bar), arr.ind = T)[[2]]]
  
  # 3. Compare to ground truth
  true_source <- as.numeric(mod_dat$truth[blk,trt,,config][1:trt])  
  
  # 4. Compute error (closest that is not already assigned)
  dist <- mod_dat$grid_dist[[config]]
  error <- dist[predicted_source, true_source]
  test <- apply(error, 1, function(x) min(x))
  error[!correct,]
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


error <- matrix(c(
  0.0, 15.0, 20.0, 25.0,  # Prediction 1: clearly closest to source 1
  18.0,  0, 22.0, 30.0,  # Prediction 2: clearly closest to source 2
  10.0, 12.0, 10.0, 35.0,  # Prediction 3: tied (10.0) between sources 1 and 3
  14.0, 16.0, 10.0, 40.0   # Prediction 4: also tied (10.0) for source 3
), nrow = 4, byrow = TRUE)

rownames(error) <- c("1", "2", "3", "4")
colnames(error) <- c("1", "2", "3", "4")

correct <- apply(error, 1, function(x) 0 %in% x)
incorrect_dist <- error[!correct, !correct]
wrong_dist <- apply(incorrect_dist, 1, min)

if(any(duplicated(wrong_dist))){
  
}
