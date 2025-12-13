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
library(RcppHungarian)
## Read in the needed functions
source(here("Code/R/backward_model/03a_BackwardGradFun.R"))
source(here("Code/R/forward_model//02a_ForwardGradFun.R"))

# Backward Fit ------------------------------------------------------------
backward_fit <- function(config, blk, trt, vst, mod_dat, inits, max_iter = 100, tol = 1e-4) {
  
  ## Extract needed data
  intensity <- mod_dat$intensity[, blk, trt, vst]
  intensity_prev <- mod_dat$intensity[, blk, trt, as.numeric(vst) - 1]
  wind <- mod_dat$wind[, , blk, trt, vst]
  dist <- mod_dat$dist
  group_id <- mod_dat$groups[, config]
  
  # Initializations
  n_src <- length(unique(na.omit(mod_dat$truth[blk, trt,,config])))
  S <- n_src
  combos <- combn(sort(unique(group_id)), S)
  K <- ncol(combos)
  prior <- rep(1 / K, K)
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
  numCores <- min(K, detectCores())
  
  # Initial E-step
  p_mat <- e_step(y = intensity, mu_mat = mu_mat, phi = theta_old['phi',], pi = prior)
  
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
        config = config, block = blk, treat = trt, visit = as.numeric(vst), n_src = n_src,
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
        config = config, block = blk, treat = trt, visit = as.numeric(vst), n_src = n_src,
        em_iters = iter, converged = TRUE,
        Q_track = list(q_track[1:iter]), Q_final = q_track[iter],
        theta = list(theta_new), p_mat = list(p_mat), pi = list(pi)
      ))
    }
    
    theta_old <- theta_new
  }
}


# Create function for distance-weighted accuracy -------------------------
dist_acc <- function(error_mat){
  
  #Just return the metric if length 1
  if(length(error_mat) == 1){
    d_pred <- error_mat
    return(d_pred)
  }
  
  #Make sure n_preds is the same as n_srcs
  if(ncol(error_mat) != nrow(error_mat)){
    print("Error: Not a square Matrix")
    stop()
  }
  
  #Assign the number of sources
  S <- ncol(error_mat)

  #Create a vector to store distances
  d_pred <- numeric(S)
  names(d_pred) <- colnames(error_mat)

  #Check if any are correct and assign them
  correct_source <- apply(error_mat, 2, function(x) 0 %in% x)
  correct_pred <- apply(error_mat, 1, function(x) 0 %in% x)
  
  if(sum(correct_source) > 0){
    d_pred[correct_source] <- 0
    
    #If they are all correct, just returnt all zeros
    if(all(correct_source)){return(d_pred)}
  }

  #Assign the remaining predictions and get their errors
  incorrect_dist <- error_mat[!correct_pred,!correct_source]
  
  #If there is only one incorrect, assign it to the remaining source
  if(length(incorrect_dist) == 1){
    d_pred[!correct_source] <- incorrect_dist
  }else{
  assignment <- HungarianSolver(incorrect_dist)
  d_pred[!correct_source] <- apply(assignment$pairs, 1, function(x) incorrect_dist[x[1], x[2]])
  }
  #Compute error metric
  return(d_pred)
}

source_pred <- function(config, blk, trt, vst, n_src, p_mat, mod_dat) {
 
  # Ensure groups is a factor with levels = 1:n_groups
  groups <- as.factor(mod_dat$groups[,config])
  levels(groups) <- sort(unique(groups))
  
  # 1. Compute p̄ matrix (average p_mat rows by destination group)
  group_ids <- sort(unique(groups))
  n_groups <- length(group_ids)

  S <- as.numeric(n_src)
  combos <- combn(sort(unique(group_ids)), S)
  K <- ncol(combos)
  
  p_bar <- matrix(NA, n_groups, K)
    
  for (g in 1:n_groups) {
    rows_in_group <- which(groups %in% group_ids[g])
    p_bar[g, ] <- colMeans(p_mat[rows_in_group, , drop = FALSE])
  }

  # 2. Predict source group for each destination group (row-wise argmax)
  predicted_source <- combos[,which(p_bar == max(p_bar, na.rm = T), arr.ind = T)[[2]]]
  
  # 3. Compare to ground truth
  true_source <- sort(unique(as.numeric(mod_dat$truth[blk,trt,,config][1:trt])))
  
  # 4. Compute error (closest that is not already assigned)
  dist <- mod_dat$grid_dist[[config]]
  error <- dist[predicted_source, true_source]
  d_pred <- dist_acc(error)
  weighted_acc <- 1 - (d_pred/max(dist))
  acc <- mean(predicted_source %in% true_source)
  n_correct <- sum(predicted_source %in% true_source)
  result <- list(
    config = config,
    block = blk,
    treat = trt,
    visit = vst,
    n_src = n_src,
    p_bar = list(p_bar),
    mean_error = mean(d_pred),
    n_correct = n_correct,
    acc = acc,
    dist_acc = mean(weighted_acc),
    component_dist_acc = list(weighted_acc)
  )
  # 5. Compute other metrics 
  return(as.data.table(result))
}
