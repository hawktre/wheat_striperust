## ---------------------------
##
## Script name: 06a_StripeSim.R
##
## Purpose of script: Run simulation for full model 
##
## Author: Trent VanHawkins
##
## Date Created: 2025-07-17
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

# Read in requried data ---------------------------------------------------
true_mod <- readRDS(here("DataProcessed/results/forward_model/fits_appended_hurdle_logistic.rds")) #Forward Model Fits
mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat.rds")) #Covariate data
raw_data <- readRDS(here("DataProcessed/experimental/stripe_clean.rds"))

# Function to simulate the data -------------------------------------------
sim_dat <- function(mu, phi){
  #convert back to shape and scale 
  alpha <- mu * phi
  beta <- phi*(1-mu)
  
  #Draw from beta distribution
  rbeta(length(mu), alpha, beta)
}


# Set up array for initial obs and simulated data ------------------------
n_plants <- max(raw_data$plant_num)
n_treat <- length(unique(raw_data$inoculum_total))
n_blocks <- length(unique(raw_data$block))
n_visits <- length(unique(raw_data$visit))

sim_array <- array(NA_real_,
                   dim = c(n_plants, n_blocks, n_treat, n_visits),
                   dimnames = list(
                     plant = paste0(1:n_plants),  # no plant names
                     block = LETTERS[1:n_blocks],
                     treat = paste0(sort(unique(raw_data$inoculum_total))),
                     visit = paste0(2:n_visits)
                   ))


# Fill in the array -------------------------------------------------------

for (plt in dimnames(sim_array)$block) {
  for (trt in dimnames(sim_array)$treat) {
    for(visit in dimnames(sim_array)$visit){
      name <- paste0(plt,trt,"_visit",visit)
      
      cur.fit <- true_mod$free$fits[[name]]
      
      mu_vec <- unlist(cur.fit$y_pred)
      phi <- unlist(cur.fit$theta)[['phi']]
      
      if(visit == "1"){
        sim_array[,plt,trt,visit] <- sim_dat(mu_vec, phi)
      }
      sim_array[,plt,trt,visit] <- sim_dat(mu_vec, phi)
    }
  }
}
