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
source(here("Code/01b_GradDescent_fun.R"))

# Fit the model -----------------------------------------------------------

fits <- list()

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
        silent = FALSE  # <-- This prints the error
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
      visit_df <- bind_rows(results) %>% arrange(neg_loglik)
      best_fit <- visit_df %>% slice(1)
      fits[[paste(plot_id, visit_name, sep = "_")]] <- best_fit
    }
  }
}

saveRDS(fits, here("DataProcessed/results/all_fits_hurdle.rds"))
