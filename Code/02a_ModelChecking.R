## ---------------------------
##
## Script name: 02a_ModelChecking.R
##
## Purpose of script: Get predictions, deviance residuals, etc.
##
## Author: Trent VanHawkins
##
## Date Created: 2025-04-30
## ---------------------------

## view outputs in non-scientific notation

options(scipen = 6, digits = 4) 

## ---------------------------

## load up the packages we will need:  (uncomment as required)

library(tidyverse)
library(here)


# Read in results & Data ---------------------------------------------------------
fits <- readRDS(here("DataProcessed/results/fits_by_gamma_max.rds"))

mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat.rds"))

## Read in necessary functions
source(here("Code/01b_GradDescent_fun.R"))


# Define function for residuals -------------------------------------------
compute_deviance_resid <- function(y, mu_hat, phi_hat) {
  # Saturated model log-likelihood: mu = y
  ll_sat <- loglik_beta(y, y, phi_hat, sum = F)
  # Fitted model log-likelihood
  ll_fit <- loglik_beta(y, mu_hat, phi_hat, sum = F)
  
  sqrt_term <- pmax(0, 2 * (ll_sat - ll_fit))
  # Deviance residuals
  sign(y - mu_hat) * sqrt(sqrt_term)
}

# Get predictions & residuals ---------------------------------------------------------
fits_200 <- fits$`200`
for (plt_visit in names(fits_200)) {
    #Subset visit and data
    fit <- fits_200[[plt_visit]]
    dat <- mod_dat[[fit$plot_id]][[fit$visit]]
    
    #Subset the parameters
    params <- unlist(fit$theta)
    
    #Create covariates
    X1 <- dat$y_prev
    X2 <- kappa_inner_sum(y_prev = dat$y_prev,
                          wind_matrix = dat$wind_mat,
                          dist_matrix = dat$dist_mat,
                          d0 = 0.01,
                          kappa = params[['kappa']],
                          derivative = F)
    
    mod_dat[[fit$plot_id]][[fit$visit]][['X2']] <- X2
    
    # Compute predictions
    y_pred <- params[['beta']] + params[['delta']]*X1 + params[['gamma']]*X2
    dev_resid <- compute_deviance_resid(y = dat$y_cur, mu_hat = inv_logit(y_pred), phi_hat = params[['phi']])
    dev <- sum(dev_resid^2)
    # Assign predictions and residuals
    fits_200[[plt_visit]]$y_pred <- list(y_pred)
    fits_200[[plt_visit]]$resid <- list(dev_resid) 
    fits_200[[plt_visit]]$deviance <- dev 
    
    #Assign th the data object for plotting
    mod_dat[[fit$plot_id]][[fit$visit]][['y_pred']] <- inv_logit(y_pred)
    mod_dat[[fit$plot_id]][[fit$visit]][['dev_resid']] <- dev_resid
}

saveRDS(fits_200, here("DataProcessed/results/fits_final.rds"))
saveRDS(mod_dat, here("DataProcessed/results/mod_dat_appended.rds"))
