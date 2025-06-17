## ---------------------------
##
## Script name: 03a_ModelChecking_Hurdle.R
##
## Purpose of script:
##
## Author: Trent VanHawkins
##
## Date Created: 2025-05-13
##
##
## ---------------------------

## view outputs in non-scientific notation

options(scipen = 6, digits = 4) 

## ---------------------------

## load up the packages we will need:  (uncomment as required)

library(tidyverse)
library(here)

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
fits <- readRDS(here("DataProcessed/results/all_fits_hurdle_logistic.rds"))

mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat.rds"))

## Read in necessary functions
source(here("Code/04a_GradDescent-Hurdle_fun.R"))


# Define function for residuals -------------------------------------------
compute_deviance_resid <- function(y, mu_hat, pi, phi_hat) {
  is_zero <- y == 0
  is_pos <- !is_zero
  
  ll_fit <- numeric(length(y))
  ll_sat <- numeric(length(y))
  
  # Zero values
  ll_fit[is_zero] <- log(pi)
  ll_sat[is_zero] <- log(pi)  # same since model fits perfectly
  
  # Positive values
  ll_fit[is_pos] <- log(1 - pi) + loglik_beta(y[is_pos], mu_hat[is_pos], phi_hat, sum = FALSE)
  ll_sat[is_pos] <- log(1 - pi) + loglik_beta(y[is_pos], y[is_pos], phi_hat, sum = FALSE)
  
  sqrt_term <- pmax(0, 2 * (ll_sat - ll_fit))
  sign(y - (1 - pi) * mu_hat) * sqrt(sqrt_term)
}

# Get predictions & residuals ---------------------------------------------------------
fit.summarise <- function(fits, mod_dat){
  for (plt_visit in names(fits)) {

    #Subset visit and data
    fit <- fits[[plt_visit]]
    dat <- mod_dat[[fit$plot_id]][[fit$visit]]
    
    #Subset the parameters
    params <- unlist(fit$theta)
    pi <- fit$pi
    
    #Create covariates
    X1 <- dat$y_prev * (1 - dat$y_prev)
    X2 <- kappa_inner_sum(y_prev = dat$y_prev,
                          wind_matrix = dat$wind_mat,
                          dist_matrix = dat$dist_mat,
                          d0 = 0.01,
                          kappa = params[['kappa']],
                          derivative = F)
    
    mod_dat[[fit$plot_id]][[fit$visit]][['X1']] <- X1
    mod_dat[[fit$plot_id]][[fit$visit]][['X2']] <- X2
    
    # Compute predictions
    eta <- params[['beta']] + params[['delta']]*X1 + params[['gamma']]*X2
    mu_hat <- inv_logit(eta)
    fitted <- (1-pi)*mu_hat
    dev_resid <- compute_deviance_resid(y = dat$y_cur, mu_hat = mu_hat, phi_hat = params[['phi']], pi = pi)
    dev <- sum(dev_resid^2)
    # Assign predictions and residuals
    fits[[plt_visit]]$y_pred <- list(fitted)
    fits[[plt_visit]]$resid <- list(dev_resid)
    fits[[plt_visit]]$deviance <- dev 
    
    #Assign th the data object for plotting
    mod_dat[[fit$plot_id]][[fit$visit]][['y_pred']] <- fitted
    mod_dat[[fit$plot_id]][[fit$visit]][['dev_resid']] <- dev_resid
    
    
  }
  return(list(fits = fits,
              mod_dat = mod_dat))
}

# Apply fit.summarise to free model
out_free <- fit.summarise(fits$free, mod_dat)

out_constrained<- list()
for (gamma_max in names(fits$constrained)) {
  fit <- fits$constrained[[gamma_max]]
  
  out_constrained[[gamma_max]] <- fit.summarise(fit, mod_dat)
} 

out_all <- list("free" = out_free, "constrained" = out_constrained)
saveRDS(out_all, here("DataProcessed/results/fits_appended_hurdle_logistic.rds"))

