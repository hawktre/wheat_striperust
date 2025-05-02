## ---------------------------
##
## Script name: 02_ModelFit.R
##
## Purpose of script: Model fitting for forward model
##
## Author: Trent VanHawkins
##
## Date Created: 2025-04-22
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
source(here("Code/01a_DataFormat_Fun.R"))
source(here("Code/01b_GradDescent_fun.R"))


# Read in the data --------------------------------------------------------
stripe <- readRDS(here("DataProcessed/experimental/stripe_clean.rds"))
wind <- readRDS(here("DataProcessed/wind/wind_clean.rds"))


# Create lagged data frame ------------------------------------------------
stripe.lag <- stripe %>% 
  arrange(plotID_new, plant_num, visit) %>% 
  group_by(plant_id) %>%
  mutate(intensity = intensity/100,
         intensity_prev = lag(intensity),
         date_prev = lag(date)) %>% 
  ungroup() %>% 
  filter(visit != "visit1")


# Create data object for modeling -----------------------------------------
## Create a df of survey periods
survey_periods <- stripe.lag %>% 
  select(plotID_new, visit, date, date_prev) %>% 
  distinct()

mod_dat <- list()

for (i in 1:nrow(survey_periods)){

  #subset the data to the current plot-visit combo
  cur.plot <- survey_periods$plotID_new[i]
  cur.visit <- survey_periods$visit[i]
  
  plot_dat <- stripe.lag %>% 
    filter(plotID_new == cur.plot,
           visit == cur.visit)
  
  #Subset y_prev
  y_prev <-  plot_dat$intensity_prev
  
  #Subset y_cur
  y_cur <- plot_dat$intensity
  
  #Get Plant ID
  plant_id <- plot_dat$plant_id
  
  #distance and direction matrices
  dist_dir <- get_dist_dir(cur.plot, plot_dat)
  
  #wind matrix
  wind_mat <- get_wind_mat(first_day = survey_periods$date_prev[i],
                           last_day = survey_periods$date[i],
                           dir.mat = dist_dir$dir,
                           wind = wind)
  #Define values of kappa to try
  kappa_try <- seq(1,5,0.5)
  
  #Initialize theta
  theta_init <- initialize_theta(y_cur = y_cur, y_prev = y_prev, wind_mat, dist_dir$dist, 
                                 d_0 = 0.01, kappa_try)
  
  mod_dat[[cur.plot]][[cur.visit]] <- list("y_prev" = y_prev,
                                       "y_cur" = y_cur,
                                       "plant_id" = plant_id,
                                       "dist_mat" = dist_dir$dist,
                                       "dir_mat" = dist_dir$dir,
                                       "wind_mat" = wind_mat,
                                       "theta_init" = theta_init)
}

saveRDS(mod_dat, here("DataProcessed/experimental/mod_dat.rds"))


# Fit the model (unconstrained) -------------------------------------------

all_fits <- list()

for (plot_id in names(mod_dat)) {
  plot_visits <- mod_dat[[plot_id]]
  
  for (visit_name in names(plot_visits)) {
    dat <- plot_visits[[visit_name]]
    results <- list()
    
    for (i in seq_along(dat$theta_init)) {
      init_theta <- dat$theta_init[[i]]
      
      fit <- tryCatch({
        optim(
          par = init_theta,
          fn = neg_loglik,
          gr = neg_grad,
          method = "BFGS",
          control = list(maxit = 5000, reltol = 1e-8),
          y_current = dat$y_cur,
          y_prev = dat$y_prev,
          wind_matrix = dat$wind_mat,
          dist_matrix = dat$dist_mat
        )
      }, error = function(e) NULL)
      
      if (!is.null(fit)) {
        results[[i]] <- tibble(
          plot_id = plot_id,
          visit = visit_name,
          iters = fit$counts,
          init_kappa = init_theta["kappa"],
          neg_loglik = fit$value,
          converged = fit$convergence == 0,
          theta = list(fit$par)
        )
      }
    }
    
    if (length(results) > 0) {
      visit_df <- bind_rows(results) %>% arrange(neg_loglik)
      best_fit <- visit_df %>% slice(1)
      all_fits[[paste(plot_id, visit_name, sep = "_")]] <- best_fit
    }
  }
}

saveRDS(all_fits, here("DataProcessed/results/fits_free.rds"))

# Fit the model (constrained) ------------------------------------------------------
# Define gamma upper bounds to try
gamma_max_vals <- c(100, 200, 500, 1000)

# Store results for each gamma max
all_fit_results <- list()

# Loop over gamma_max values
for (gamma_max in gamma_max_vals) {
  
  message("Fitting models with gamma_max = ", gamma_max)
  
  # Define lower and upper bounds (phi capped at 200)
  lower_bounds <- c(-Inf, -Inf, -Inf, -Inf, 0.01)
  upper_bounds <- c(Inf, Inf, gamma_max, Inf, Inf)
  
  all_fits <- list()
  
  for (plot_id in names(mod_dat)) {
    plot_visits <- mod_dat[[plot_id]]
    
    for (visit_name in names(plot_visits)) {
      dat <- plot_visits[[visit_name]]
      results <- list()
      
      for (i in seq_along(dat$theta_init)) {
        init_theta <- dat$theta_init[[i]]
        
        fit <- tryCatch({
          optim(
            par = init_theta,
            fn = neg_loglik,
            gr = neg_grad,
            method = "L-BFGS-B",
            lower = lower_bounds,
            upper = upper_bounds,
            control = list(maxit = 5000),
            y_current = dat$y_cur,
            y_prev = dat$y_prev,
            wind_matrix = dat$wind_mat,
            dist_matrix = dat$dist_mat
          )
        }, error = function(e) NULL)
        
        if (!is.null(fit)) {
          results[[i]] <- tibble(
            plot_id = plot_id,
            visit = visit_name,
            gamma_max = gamma_max,  # track which gamma cap was used
            iters = fit$counts,
            init_kappa = init_theta["kappa"],
            neg_loglik = fit$value,
            converged = fit$convergence == 0,
            theta = list(fit$par)
          )
        }
      }
      
      if (length(results) > 0) {
        visit_df <- bind_rows(results) %>% arrange(neg_loglik)
        best_fit <- visit_df %>% slice(1)
        all_fits[[paste(plot_id, visit_name, sep = "_")]] <- best_fit
      }
    }
  }
  
  # Store results for this gamma_max
  all_fit_results[[as.character(gamma_max)]] <- all_fits
}

# Save all fits from all gamma_max values
write_rds(all_fit_results, here("DataProcessed/results/fits_by_gamma_max.rds"))

greek_cols <- wesanderson::wes_palette("Darjeeling1", 5, type = "discrete")


# Plot constrained fits ---------------------------------------------------

