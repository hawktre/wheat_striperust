library(here)
library(data.table)
library(dplyr)
library(purrr)
library(parallel)  # Use base parallel package
source(here("Code/02a_GradDescentFun.R"))
source(here("Code/03a_EMgradfun.R"))
source(here("Code/04a_SimFunc.R"))

# Read in data
forward_fits <- readRDS(here("DataProcessed/results/forward_model/forward_fits.rds"))
mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat_arrays.rds"))

## ---------------------------
##
## Script name: 04a_SimFunc.R
##
## Purpose of script: Compile all necessary functions for simulation.
##
## Author: Trent VanHawkins
##
## Date Created: 2025-07-22
##
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## ---------------------------

## view outputs in non-scientific notation

options(scipen = 6, digits = 4) 

## ---------------------------

## load up the packages we will need:  (uncomment as required)

# Function to simulate the data -------------------------------------------
disease_sim <- function(pars, mu, pi) {
  
  phi <- pars[["phi"]]
  
  # Convert to shape parameters
  alpha <- mu * phi
  beta <- phi * (1 - mu)
  
  # Minimum perceptible value
  min_detectable <- 0.0001
  
  # Draw from beta and clip
  sim_dat <- map2_dbl(alpha, beta, ~max(min_detectable, rbeta(1, shape1 = .x, shape2 = .y)))
  
  # Apply zero inflation
  zeros <- rbinom(length(sim_dat), size = 1, prob = 1 - pi)
  
  return(sim_dat * zeros)
}

# Forward Fit -------------------------------------------------------------
forward_fit <- function(blk, trt, vst, mod_dat, dist, kappa_try) {
  browser()
  ## Extract Needed Model Data
  intensity <- mod_dat$intensity[, blk, trt, vst]
  intensity_prev <- mod_dat$intensity[, blk, trt, as.numeric(vst) - 1]
  wind <- mod_dat$wind[, , blk, trt, vst]
  pi <- mean(intensity == 0)
  non_zero <- which(intensity != 0)
  
  ## Set up strage for temproary results
  results_list <- list()
  
  for (kappa in kappa_try) {
    #Generate initial values for current kappa in grid search
    init_theta <- initialize_theta(y_cur = intensity[non_zero], 
                                   y_prev = intensity_prev[non_zero], 
                                   wind_mat = wind[non_zero, non_zero], 
                                   dist_mat = dist[non_zero, non_zero], 
                                   kappa_try = as.numeric(kappa),
                                   d_0 = 0.01)
    
     fit <- tryCatch(
      optim(
        par = init_theta,
        fn = neg_loglik,
        gr = neg_grad,
        method = "BFGS",
        control = list(maxit = 5000, reltol = 1e-8),
        y_current = intensity[non_zero],
        y_prev = intensity_prev[non_zero],
        wind_matrix = wind[non_zero, non_zero],
        dist_matrix = dist[non_zero, non_zero]
      ),
      error = function(e) {
        message(sprintf("forward_fit failed: %s [blk=%s, trt=%s, vst=%s, kappa=%s]",
                        conditionMessage(e), blk, trt, vst, kappa))
        return(NULL)
      }
    )
    
    if (!is.null(fit)) {
      params <- fit$par
      
      #Create covariates with the data
      X1 <- intensity_prev * (1 - intensity_prev)
      X2 <- kappa_inner_sum(y_prev = intensity_prev,
                            wind_matrix = wind,
                            dist_matrix = dist,
                            d0 = 0.01,
                            kappa = params[['kappa']],
                            derivative = F)
      
      # Compute predictions
      eta <- params[['beta']] + params[['delta']]*X1 + params[['gamma']]*X2
      mu_hat <- inv_logit(eta)
      fitted <- (1-pi)*mu_hat
      
      #Compute MSE
      mse <- sum((intensity-fitted)^2)/length(intensity)
      
      results_list[[length(results_list) + 1]] <- list(
        block = blk,
        treat = as.numeric(trt),
        visit = as.numeric(vst),
        iters = fit$counts[2],
        init_kappa = as.numeric(kappa),
        neg_loglik = fit$value,
        MSE = mse,
        converged = fit$convergence == 0,
        theta = list(fit$par),
        pi = pi
      )
    }
    
  }
  
  if (length(results_list) > 0) {
    visit_dt <- rbindlist(results_list)
    visit_dt <- visit_dt[neg_loglik > -1e5 & converged == T]
    if (nrow(visit_dt) > 0) {
      return(visit_dt[which.min(neg_loglik)])
    }
  }
  
  return(NULL)
}


# Backward Fit ------------------------------------------------------------
backward_fit <- function(config, blk, trt, vst, mod_dat, forward_fits) {
  browser()
  # Setup
  intensity <- mod_dat$intensity[, blk, trt, vst]
  intensity_prev <- mod_dat$intensity[, blk, trt, as.numeric(vst) - 1]
  wind <- mod_dat$wind[, , blk, trt, vst]
  dist <- mod_dat$dist
  groups <- mod_dat$groups[, config]
  non_zero <- which(intensity != 0)
  n_groups <- readr::parse_number(config)
  prior <- 1 / n_groups
  
  # Starting values from forward fit
  theta <- forward_fits[
    block == blk & treat == as.numeric(trt) & visit == as.numeric(vst)
  ][["theta"]][[1]]
  
  max_em_iter <- 1000
  tol <- 1e-8
  q_track <- numeric(max_em_iter)
  
  # Initial E-step
  p_mat <- e_step(
    par = theta,
    y_current = intensity,
    y_prev = intensity_prev,
    wind_matrix = wind,
    dist_matrix = dist,
    group_id = groups,
    prior = prior
  )
  
  # EM loop
  for (em_iter in seq_len(max_em_iter)) {
    fit <- tryCatch(
      optim(
        par = theta,
        fn = wrapped_obj,
        gr = mstep_grad_em,
        method = "BFGS",
        control = list(maxit = 1000, reltol = 1e-8),
        p_mat = p_mat[non_zero, ],
        y_current = intensity[non_zero],
        y_prev = intensity_prev[non_zero],
        wind_matrix = wind[non_zero, non_zero],
        dist_matrix = dist[non_zero, non_zero],
        group_id = groups[non_zero]
      ),
      error = function(e) {
        msg <- sprintf(
          "EM step %d failed in backward_fit [config=%s, block=%s, treat=%s, visit=%s]: %s",
          em_iter, config, blk, trt, vst, conditionMessage(e)
        )
        warning(msg)
        
        return(NULL)
      }
    )
    
    if (is.null(fit) || fit$convergence != 0 || !is.finite(fit$value)) {
      return(data.table(
        config = config,
        block = blk,
        treat = as.numeric(trt),
        visit = as.numeric(vst),
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
    
    q_track[em_iter] <- q_val
    
    if (em_iter > 1) {
      diff <- abs(q_track[em_iter] - q_track[em_iter - 1])
      if(diff < tol){break}
    }
    
    theta <- theta_new
  }
  
  # Return a single-row data.table
  result <- data.table(
    config = config,
    block = blk,
    treat = as.numeric(trt),
    visit = as.numeric(vst),
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


# Wrapper Function for the simulation -------------------------------------
single_sim <- function(sim_id, dat, forward_mod, output_dir = here("DataProcessed/results/simulation")) {

  base_seed <- 404
  set.seed(base_seed + sim_id)
    # Set up indices
    blocks <- dimnames(dat$intensity)[["block"]]
    treats <- dimnames(dat$intensity)[["treat"]]
    visits <- dimnames(dat$intensity)[["visit"]]
    kappa_try <- dimnames(dat$inits)[["kappa_try"]]
    configs <- dimnames(dat$groups)[["config"]]
    
    # 1. Simulate new intensity values
    intensity_sim <- dat$intensity
    for (blk in blocks) {
      for (trt in treats) {
        for (vst in visits[-1]) {
          fit <- forward_mod[block == blk & treat == as.numeric(trt) & visit == as.numeric(vst)]
          pars <- fit[["theta"]][[1]]
          pi <- fit[["pi"]]
          fitted <- fit[["fitted"]][[1]]
          intensity_sim[, blk, trt, vst] <- disease_sim(pars, fitted, pi)
        }
      }
    }
    dat$intensity <- intensity_sim
    
    # 3. Fit Forward Model
    forward_fit("A", "1", "2", dat, dat$dist, kappa_try)
    combos_forward <- expand.grid(block = blocks, treat = treats, visit = visits[-1], stringsAsFactors = FALSE)
    forward <- pmap(combos_forward, ~forward_fit(..1, ..2, ..3, dat, dat$dist, kappa_try)) %>% rbindlist()
    
    # 2. Fit backward model
    combos_backward <- expand.grid(config = configs, bk = blocks, trt = treats, vst = visits[-1], stringsAsFactors = FALSE)
    backward <- pmap(combos_backward, ~backward_fit(..1, ..2, ..3, ..4, dat, forward)) %>% rbindlist()
    
    # 4. Source prediction for treat == 1
    backward_t1 <- backward[treat == 1]
    sources_predicted <- pmap(backward_t1, ~source_pred(config = ..1,
                                                        blk = ..2,
                                                        trt = ..3,
                                                        vst = ..4,
                                                        p_mat = ..9,
                                                        mod_dat = dat)) %>% 
      rbindlist()
    
    results_merge <- left_join(backward, forward, by = c("block", "treat", "visit"), suffix = c(".backward", ".forward")) %>% 
      dplyr::select(-c(final_theta, p_mat)) %>% 
      left_join(sources_predicted, by = c("config", "block", "treat", "visit")) %>% 
      mutate(sim = sim_id) %>% 
      dplyr::select(sim, everything())
    # 5. Return results
    return(results_merge)
}


test <- single_sim(1739, mod_dat, forward_fits)
