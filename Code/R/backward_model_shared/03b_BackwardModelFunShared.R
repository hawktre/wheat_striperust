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
source(here(
  "Code/R/backward_model_shared/03a_BackwardGradFunShared.R"
))
source(here("Code/R/forward_model/02a_ForwardGradFun.R"))

# Backward Fit ------------------------------------------------------------
backward_fit <- function(
  config,
  blk,
  trt,
  vst,
  n_src,
  mod_dat,
  inits,
  max_iter = 100,
  tol = 1e-4
) {
  ## Extract needed data
  intensity <- mod_dat$intensity[, blk, trt, vst]
  intensity_prev <- mod_dat$intensity[, blk, trt, as.numeric(vst) - 1]
  wind <- mod_dat$wind[,, blk, trt, vst]
  dist <- mod_dat$dist
  group_id <- mod_dat$groups[, config]

  # Initializations
  S <- n_src
  combos <- combn(sort(unique(group_id)), S)
  K <- ncol(combos)
  prior <- rep(1 / K, K)

  # Pre-compute component indicator matrix ONCE
  component_indicator <- create_component_indicator(group_id, combos)

  ## Tracking
  max_em_iter <- max_iter
  observed_ll <- numeric(max_em_iter)

  # Initial E-step
  theta_old <- inits
  E <- e_step(
    par = theta_old,
    y_current = intensity,
    y_prev = intensity_prev,
    wind = wind,
    dist = dist,
    group_id = group_id,
    components = combos,
    pi_vec = prior,
    component_indicator = component_indicator
  )

  ## Initial Q and log-likelihood
  if (!is.na(E$ll_obs)) {
    observed_ll[1] <- E$ll_obs
  }

  # EM loop
  for (iter in 2:max_em_iter) {
    # M-step
    M <- m_step(
      theta_old = theta_old,
      intensity = intensity,
      intensity_prev = intensity_prev,
      wind = wind,
      dist = dist,
      group_id = group_id,
      p_mat = E$p_mat,
      components = combos,
      component_indicator = component_indicator
    )

    ## E-step: recompute responsibilities
    E <- e_step(
      par = M$theta_new,
      y_current = intensity,
      y_prev = intensity_prev,
      wind = wind,
      dist = dist,
      group_id = group_id,
      components = combos,
      pi_vec = M$pi,
      component_indicator = component_indicator
    )

    ## Compute observed log-likelihood and Q
    observed_ll[iter] <- E$ll_obs

    # Compute BIC
    n_obs <- length(intensity)
    n_params <- 5 + K - 1 # 5 shared parameters + K mixing weights
    bic <- -2 * (-observed_ll[iter]) + n_params * log(n_obs)

    ## Convergence checks
    if (!is.finite(observed_ll[iter]) || iter == max_em_iter) {
      return(data.table(
        config = config,
        block = blk,
        treat = trt,
        visit = as.numeric(vst),
        n_src = S,
        em_iters = iter,
        converged = FALSE,
        loglik = observed_ll[iter],
        bic = bic,
        n_params = n_params,
        theta = list(M$theta_new),
        p_mat = list(E$p_mat)
      ))
    }

    theta_diff <- max(abs(M$theta_new - theta_old))
    observed_ll_diff <- abs(observed_ll[iter] - observed_ll[iter - 1])
    if (observed_ll_diff < tol || theta_diff < tol) {
      return(data.table(
        config = config,
        block = blk,
        treat = trt,
        visit = as.numeric(vst),
        n_src = S,
        em_iters = iter,
        converged = TRUE,
        loglik = observed_ll[iter],
        bic = bic,
        n_params = n_params,
        theta = list(M$theta_new),
        p_mat = list(E$p_mat)
      ))
    }

    theta_old <- M$theta_new
  }
}


# Create function for distance-weighted accuracy -------------------------
dist_acc <- function(error_mat) {
  #Just return the metric if length 1
  if (length(error_mat) == 1) {
    d_pred <- error_mat
    return(d_pred)
  }

  #Make sure n_preds is the same as n_srcs
  if (ncol(error_mat) != nrow(error_mat)) {
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

  if (sum(correct_source) > 0) {
    d_pred[correct_source] <- 0

    #If they are all correct, just returnt all zeros
    if (all(correct_source)) {
      return(d_pred)
    }
  }

  #Assign the remaining predictions and get their errors
  incorrect_dist <- error_mat[!correct_pred, !correct_source]

  #If there is only one incorrect, assign it to the remaining source
  if (length(incorrect_dist) == 1) {
    d_pred[!correct_source] <- incorrect_dist
  } else {
    assignment <- HungarianSolver(incorrect_dist)
    d_pred[!correct_source] <- apply(assignment$pairs, 1, function(x) {
      incorrect_dist[x[1], x[2]]
    })
  }
  #Compute error metric
  return(d_pred)
}

source_pred <- function(config, blk, trt, vst, n_src, p_mat, mod_dat) {
  intensity <- mod_dat$intensity[, blk, trt, vst]

  # Ensure groups is a factor with levels = 1:n_groups
  groups <- as.factor(mod_dat$groups[, config])
  levels(groups) <- sort(unique(groups))

  # 1. Compute p̄ matrix (average p_mat rows by destination group)
  group_ids <- sort(unique(groups))
  n_groups <- length(group_ids)

  #Set up number of sources to predict
  S <- as.numeric(n_src)
  combos <- combn(sort(unique(group_ids)), S)
  K <- ncol(combos)

  p_bar <- rowsum(p_mat, groups) / tabulate(groups)

  # 2. Predict source group for each destination group (row-wise argmax)
  predicted_source <- combos[, which(
    p_bar == max(p_bar, na.rm = T),
    arr.ind = T
  )[[2]]]

  # 2b. Naive prediction: top S groups by mean intensity
  group_max_intensity <- as.numeric(tapply(intensity, groups, max))
  naive_source <- as.numeric(group_ids[order(
    group_max_intensity,
    decreasing = TRUE
  )[1:S]])

  # 3. Compare to ground truth
  true_source <- sort(unique(as.numeric(mod_dat$truth[blk, trt, , config][
    1:trt
  ])))

  # 4. Helper to compute accuracy metrics given a predicted source
  compute_metrics <- function(predicted) {
    if (length(predicted) == length(true_source)) {
      dist <- mod_dat$grid_dist[[config]]
      error <- dist[predicted, true_source]
      d_pred <- dist_acc(error)
      weighted_acc <- 1 - (d_pred / max(dist))
      acc <- mean(predicted %in% true_source)
      n_correct <- sum(predicted %in% true_source)
    } else {
      d_pred <- NA
      weighted_acc <- NA
      acc <- NA
      n_correct <- NA
    }
    list(
      mean_error = mean(d_pred),
      n_correct = n_correct,
      acc = acc,
      dist_acc = mean(weighted_acc),
      component_dist_acc = list(weighted_acc),
      predicted_source = list(predicted)
    )
  }

  model_metrics <- compute_metrics(predicted_source)
  naive_metrics <- compute_metrics(naive_source)

  result <- data.table(
    config = config,
    block = blk,
    treat = trt,
    visit = as.numeric(vst),
    n_src = n_src,
    true_source = list(true_source),
    # Model predictions
    predicted_source = model_metrics$predicted_source,
    mean_error = model_metrics$mean_error,
    n_correct = model_metrics$n_correct,
    acc = model_metrics$acc,
    dist_acc = model_metrics$dist_acc,
    component_dist_acc = model_metrics$component_dist_acc,
    # Naive predictions
    naive_source = naive_metrics$predicted_source,
    naive_mean_error = naive_metrics$mean_error,
    naive_n_correct = naive_metrics$n_correct,
    naive_acc = naive_metrics$acc,
    naive_dist_acc = naive_metrics$dist_acc,
    naive_component_dist_acc = naive_metrics$component_dist_acc
  )

  return(result)
}
