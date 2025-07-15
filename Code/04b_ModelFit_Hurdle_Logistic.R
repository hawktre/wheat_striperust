## ---------------------------
##
## Script name: 04_ModelFit_Hurdle_Logistic.R
##
## Purpose of script: Fit model with logistic autoinfection term
##
## Author: Trent VanHawkins
##
## Date Created: 2025-05-16
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

library(tidyverse)
library(here)

## ---------------------------
##
## Script name: 03_ModelFit_ZIB.R
##
## Purpose of script: Fit a zero-inflated beta model 
##
## Author: Trent VanHawkins
##
## Date Created: 2025-05-13
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

library(tidyverse)
library(here)

# Read in the data --------------------------------------------------------
mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat.rds"))
source(here("Code/01a_DataFormat_Fun.R"))
source(here("Code/04a_GradDescent-Hurdle_fun.R"))


# Fit the model -----------------------------------------------------------

free_fits <- list()

for (plot_id in names(mod_dat)) {
  plot_visits <- mod_dat[[plot_id]]
  
  for (visit_name in names(plot_visits)) {
    dat <- plot_visits[[visit_name]]
    results <- list()
    
    for (i in seq_along(dat$theta_init)) {
      
      
      init_theta <- dat$theta_init[[i]]
      
      pi <- mean(dat$y_cur == 0)
      non_zero <- which(dat$y_cur != 0)
      
      fit <- try(
        optim(
          par = init_theta,
          fn = neg_loglik,
          gr = neg_grad,
          method = "BFGS",
          control = list(maxit = 5000, reltol = 1e-8),
          y_current = dat$y_cur[non_zero],
          y_prev = dat$y_prev[non_zero],
          wind_matrix = dat$wind_mat[non_zero, non_zero],
          dist_matrix = dat$dist_mat[non_zero, non_zero]
        ),
        silent = FALSE  
      )
      
      if (!inherits(fit, "try-error")) {
        results[[i]] <- tibble(
          plot_id = plot_id,
          visit = visit_name,
          iters = fit$counts,
          init_kappa = init_theta["kappa"],
          neg_loglik = fit$value,
          converged = fit$convergence == 0,
          theta = list(fit$par),
          pi = pi
        )
      }
    }
    
    if (length(results) > 0) {
      visit_df <- bind_rows(results) %>% filter(neg_loglik > -1e5) %>% arrange(neg_loglik) 
      best_fit <- visit_df %>% slice(1)
      free_fits[[paste(plot_id, visit_name, sep = "_")]] <- best_fit
    }
  }
}





# Fit the model (constrained) ------------------------------------------------------
# Define gamma upper bounds to try
gamma_max_vals <- c(100, 200, 500, 1000, 2000, 5000)

# Store results for each gamma max
constrained_fits <- list()

for (gamma_max in gamma_max_vals) {
  message("Fitting models with gamma_max = ", gamma_max)
  
  # Define bounds (allowing unconstrained except for gamma and phi)
  lower_bounds <- c(-Inf, -Inf, -Inf, 0.01, 0.01)
  upper_bounds <- c(Inf, Inf, gamma_max, Inf, Inf)
  
  fits <- list()
  
  for (plot_id in names(mod_dat)) {
    plot_visits <- mod_dat[[plot_id]]
    
    for (visit_name in names(plot_visits)) {
      dat <- plot_visits[[visit_name]]
      results <- list()
      
      for (i in seq_along(dat$theta_init)) {
        init_theta <- dat$theta_init[[i]]
        
        # Hurdle-specific preprocessing
        pi <- mean(dat$y_cur == 0)
        non_zero <- which(dat$y_cur != 0)
        
        fit <- tryCatch({
          optim(
            par = init_theta,
            fn = neg_loglik,
            gr = neg_grad,
            method = "L-BFGS-B",
            lower = lower_bounds,
            upper = upper_bounds,
            control = list(maxit = 5000),
            y_current = dat$y_cur[non_zero],
            y_prev = dat$y_prev[non_zero],
            wind_matrix = dat$wind_mat[non_zero, non_zero],
            dist_matrix = dat$dist_mat[non_zero, non_zero]
          )
        }, error = function(e) NULL)
        
        if (!is.null(fit)) {
          results[[i]] <- tibble(
            plot_id = plot_id,
            visit = visit_name,
            gamma_max = gamma_max,
            iters = fit$counts,
            init_kappa = init_theta["kappa"],
            neg_loglik = fit$value,
            converged = fit$convergence == 0,
            theta = list(fit$par),
            pi = pi
          )
        }
      }
      
      if (length(results) > 0) {
        visit_df <- bind_rows(results) %>% filter(neg_loglik > -1e5) %>% arrange(neg_loglik)
        best_fit <- visit_df %>% slice(1)
        fits[[paste(plot_id, visit_name, sep = "_")]] <- best_fit
      }
    }
  }
  
  # Save best fit results for this gamma_max value
  constrained_fits[[paste0("gamma_", gamma_max)]] <- fits
}

all_fits <- list("free" = free_fits, "constrained" = constrained_fits)
saveRDS(all_fits, here("DataProcessed/results/forward_model/all_fits_hurdle_logistic.rds"))

