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
source(here("Code/R/backward_model_shared/03a_BackwardGradFunShared_vectorized.R"))
source(here("Code/R/forward_model/02a_ForwardGradFun.R"))

# Backward Fit ------------------------------------------------------------
backward_fit <- function(config, blk, trt, vst, mod_dat, inits, max_iter = 100, tol = 1e-4) {
  
  ## Extract needed data
  intensity <- mod_dat$intensity[, blk, trt, vst]
  intensity_prev <- mod_dat$intensity[, blk, trt, as.numeric(vst) - 1]
  wind <- mod_dat$wind[, , blk, trt, vst]
  dist <- mod_dat$dist
  group_id <- mod_dat$groups[, config]
  
  # Initializations
  S <- length(unique(na.omit(mod_dat$truth[blk, trt,,config])))
  combos <- combn(sort(unique(group_id)), S)
  K <- ncol(combos)
  prior <- rep(1 / K, K)
  
  ## Tracking
  max_em_iter <- max_iter
  observed_ll <- numeric(max_em_iter)
  
  # Initial E-step
  theta_old <- inits
  E <- e_step(par = theta_old, 
                    y_current = intensity, 
                    y_prev = intensity_prev, 
                    wind = wind, 
                    dist = dist, 
                    group_id = group_id, 
                    components = combos,
                  pi_vec = prior)
  
  ## Initial Q and log-likelihood
  if(!is.na(E$ll_obs)){
    observed_ll[1] <- E$ll_obs
  }


  # EM loop
  for (iter in 2:max_em_iter) {

    # M-step
    M <- m_step(theta_old = theta_old, 
                    intensity = intensity,
                    intensity_prev = intensity_prev, 
                    wind = wind, 
                    dist = dist, 
                    group_id = group_id, 
                    p_mat = E$p_mat, 
                    components = combos)
    
    ## E-step: recompute responsibilities
    E <- e_step(par = M$theta_new, 
                    y_current = intensity, 
                    y_prev = intensity_prev, 
                    wind = wind, 
                    dist = dist, 
                    group_id = group_id, 
                    components = combos,
                  pi_vec = M$pi)
    
    ## Compute observed log-likelihood and Q
    observed_ll[iter] <- E$ll_obs
    
    
    ## Convergence checks
    if (!is.finite(observed_ll[iter]) || iter == max_em_iter) {
      return(data.table(
        config = config, block = blk, treat = trt, visit = as.numeric(vst), n_src = S,
        em_iters = iter, converged = FALSE, Q_final = observed_ll[iter],
        theta = list(M$theta_new), p_mat = list(E$p_mat)
      ))
    }
    
    theta_diff <- max(abs(M$theta_new - theta_old))
    observed_ll_diff <- abs(observed_ll[iter] - observed_ll[iter - 1])
    if (observed_ll_diff < tol || theta_diff < tol) {
      return(data.table(
        config = config, block = blk, treat = trt, visit = as.numeric(vst), n_src = S,
        em_iters = iter, converged = TRUE, Q_final = observed_ll[iter],
        theta = list(M$theta_new), p_mat = list(E$p_mat)
      ))
    }
    
    theta_old <- M$theta_new
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
    mean_error = mean(d_pred),
    n_correct = n_correct,
    acc = acc,
    dist_acc = mean(weighted_acc),
    component_dist_acc = list(weighted_acc)
  )
  # 5. Compute other metrics 
  return(as.data.table(result))
}
