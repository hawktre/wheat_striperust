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

library(tidyverse)
library(here)

# Read in the data --------------------------------------------------------
em_dat <- readRDS(here("DataProcessed/experimental/em_dat.rds"))

## Read in the needed functions
source(here("Code/05b_EMgradfun.R"))


em_summary <- list()
# Run all models in a for-loop --------------------------------------------
for (config in names(em_dat)) {
  for(plt in names(em_dat[[config]])){
    for(visit in names(em_dat[[config]][[plt]])){
    cat("Config:", config, "Plot:", plt, "Visit:", visit)
      # if(config == "stripe_8h" && plt == "D2" && visit == "visit4"){
      #   browser()
      # }
      # Specifications for the model -------------------------------------------------
      dat <- em_dat[[config]][[plt]][[visit]]
      n_groups <- length(unique(dat$group_id))
    
      ## Figure out which ones are zero (for hurdle-model)
      non_zero <- which(dat$y_cur != 0)
      
      ## set up Parameters
      max_em_iter <- 1000
      tol <- 1e-8
      theta <- dat$inits
      q_track <- numeric(max_em_iter)
      theta_track <- list()
      prior <- 1/n_groups
      
      # Initial E-step
      p_mat <- e_step(par = theta,
                      y_current = dat$y_cur,
                      y_prev = dat$y_prev,
                      wind_matrix = dat$wind_mat,
                      dist_matrix = dat$dist_mat,
                      group_id = dat$group_id,
                      prior = prior)
      
      # Loop
      for (em_iter in seq_len(max_em_iter)) {
        
        # M-step
        fit <- try(
          optim(
            par = theta,
            fn = wrapped_obj,
            gr = mstep_grad_em,
            method = "BFGS",
            control = list(maxit = 1000, reltol = 1e-8),
            p_mat = p_mat[non_zero,],
            y_current = dat$y_cur[non_zero],
            y_prev = dat$y_prev[non_zero],
            wind_matrix = dat$wind_mat[non_zero, non_zero],
            dist_matrix = dat$dist_mat[non_zero, non_zero],
            group_id = dat$group_id[non_zero]
          ),
          silent = TRUE
        )
        
        if (inherits(fit, "try-error")) {
          warning(paste("EM step", em_iter, "failed"))
          break
        }
        
        theta_new <- fit$par
        theta_track[[em_iter]] <- theta_new
        
        # Update E-step
        p_mat <- e_step(par = theta_new,
                        y_current = dat$y_cur,
                        y_prev = dat$y_prev,
                        wind_matrix = dat$wind_mat,
                        dist_matrix = dat$dist_mat,
                        group_id = dat$group_id,
                        prior = prior)
        
        # Compute Q
        mu_mat <- get_mu(par = theta_new,
                         y_prev = dat$y_prev,
                         wind_matrix = dat$wind_mat,
                         dist_matrix = dat$dist_mat,
                         group_id = dat$group_id)
        
        q_val <- Q_fun(y = dat$y_cur, mu_mat = mu_mat, phi = theta_new[['phi']], p_mat = p_mat)
        q_track[em_iter] <- q_val
        
        # Convergence check
        if (em_iter > 1 && abs(q_track[em_iter] - q_track[em_iter - 1]) < tol) {
          message("EM converged at iteration ", em_iter)
          break
        }
        
        theta <- theta_new
      }
      em_summary[[config]][[plt]][[visit]] <- list(em_iters = em_iter,
                                                   converged = em_iter < max_em_iter,
                                                   final_neg_loglik = q_track[em_iter],
                                                   final_theta = theta_new,
                                                   p_mat = p_mat)
    }
  }
}

saveRDS(em_summary, file = "DataProcessed/results/backward_model/backwards_mle_fits.rds")



# Sensitivity Analysis ----------------------------------------------------
init.length <- length(em_dat$stripe_4$A1$visit2$og_inits)
em_sensitivity <- list()
# Run all models in a for-loop --------------------------------------------
for (inits in 1:init.length){
  for (config in names(em_dat)){
    for(plt in names(em_dat[[config]])){
      for(visit in names(em_dat[[config]][[plt]])){
        
        cat("Config:", config, "Plot:", plt, "Visit:", visit, "Kappa:", inits)
        dat <- em_dat[[config]][[plt]][[visit]]
        n_groups <- length(unique(dat$group_id))
        
        ## Figure out which ones are zero (for hurdle-model)
        non_zero <- which(dat$y_cur != 0)
        
        ## set up Parameters
        max_em_iter <- 1000
        tol <- 1e-8
        theta <- dat[["og_inits"]][[inits]]
        kappa_try <- theta[["kappa"]]
        q_track <- numeric(max_em_iter)
        theta_track <- list()
        prior <- 1/n_groups
        
        # Initial E-step
        p_mat <- e_step(par = theta,
                        y_current = dat$y_cur,
                        y_prev = dat$y_prev,
                        wind_matrix = dat$wind_mat,
                        dist_matrix = dat$dist_mat,
                        group_id = dat$group_id,
                        prior = prior)
        
        # Loop
        for (em_iter in seq_len(max_em_iter)) {
          
          # M-step
          fit <- try(
            optim(
              par = theta,
              fn = wrapped_obj,
              gr = mstep_grad_em,
              method = "BFGS",
              control = list(maxit = 1000, reltol = 1e-8),
              p_mat = p_mat[non_zero,],
              y_current = dat$y_cur[non_zero],
              y_prev = dat$y_prev[non_zero],
              wind_matrix = dat$wind_mat[non_zero, non_zero],
              dist_matrix = dat$dist_mat[non_zero, non_zero],
              group_id = dat$group_id[non_zero]
            ),
            silent = TRUE
          )
          
          if (inherits(fit, "try-error")) {
            warning(paste("EM step", em_iter, "failed"))
            break
          }
          
          theta_new <- fit$par
          theta_track[[em_iter]] <- theta_new
          
          # Update E-step
          p_mat <- e_step(par = theta_new,
                          y_current = dat$y_cur,
                          y_prev = dat$y_prev,
                          wind_matrix = dat$wind_mat,
                          dist_matrix = dat$dist_mat,
                          group_id = dat$group_id,
                          prior = prior)
          
          # Compute Q
          mu_mat <- get_mu(par = theta_new,
                           y_prev = dat$y_prev,
                           wind_matrix = dat$wind_mat,
                           dist_matrix = dat$dist_mat,
                           group_id = dat$group_id)
          
          q_val <- Q_fun(y = dat$y_cur, mu_mat = mu_mat, phi = theta_new[['phi']], p_mat = p_mat)
          q_track[em_iter] <- q_val
          
          # Convergence check
          if (em_iter > 1 && abs(q_track[em_iter] - q_track[em_iter - 1]) < tol) {
            message("EM converged at iteration ", em_iter)
            break
          }
          
          theta <- theta_new
        }
        em_sensitivity[[config]][[plt]][[visit]][[paste0("kappa_", kappa_try)]] <- list(em_iters = em_iter,
                                                     converged = em_iter < max_em_iter,
                                                     final_neg_loglik = q_track[em_iter],
                                                     final_theta = theta_new,
                                                     p_mat = p_mat)
      }
    }
  }
}

saveRDS(em_sensitivity, file = "DataProcessed/results/backward_model/backwards_sensitivity_fits.rds")

