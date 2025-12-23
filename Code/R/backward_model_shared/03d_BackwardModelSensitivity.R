## ---------------------------
##
## Script name: 03c_BackwardModelFit.R
##
## Purpose of script: Fit the backward source-prediction model
##
## Author: Trent VanHawkins
##
## Date Created: 2025-08-20
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

library(here)
library(tidyverse)
library(data.table)

## Read in necessary functions
source(here("Code/R/backward_model_shared/03b_BackwardModelFunShared.R"))

## Read in forward fits
forward <- readRDS(here("DataProcessed/results/forward_model/forward_fits.rds"))
mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat_arrays.rds"))

# Set up Indices ----------------------------------------------------------
blocks <- dimnames(mod_dat$intensity)[["block"]]
treats <- dimnames(mod_dat$intensity)[["treat"]]
visits <- dimnames(mod_dat$intensity)[["visit"]]
configs <- dimnames(mod_dat$groups)[["config"]]
kappa_try <- c(0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0)


# Get the inits for each value of kappa -----------------------------------

# Fit the model -----------------------------------------------------------
# 2. Fit backward model
combos_backward <- expand.grid(config = configs, blk = blocks, trt = treats, vst = visits[-1], kappa = kappa_try, stringsAsFactors = FALSE)
combos_backward$init <- pmap(combos_backward %>% select(-config), function(blk, trt, vst, kappa) {
  intensity <- mod_dat$intensity[, blk, trt, vst]
  intensity_prev <- mod_dat$intensity[, blk, trt, as.numeric(vst) - 1]
  wind <- mod_dat$wind[,, blk, trt, vst]
  dist <- mod_dat$dist
  
  return(initialize_theta(
    y = intensity,
    y_prev = intensity_prev,
    wind_mat = wind,
    dist_mat = dist,
    kappa = kappa,
    d_0 = 0.01
  ))
})



start <- Sys.time()
backward <- lapply(seq_len(1), function(k) {
  #Subset the current combination
  combo <- combos_backward[k,]
  #Run the fitting
  tmp <- backward_fit(config = combo$config,
  blk = combo$blk,
  trt = combo$trt,
  vst = combo$vst,
  inits = combo$init[[1]],
  mod_dat = mod_dat, 
  max_iter = 2000,
  tol = 1e-4)
  #Annotate which initial value of kappa was used
  tmp <- tmp |> mutate(kappa = combo$kappa)
  return(tmp)
}) |> rbindlist()
end <- Sys.time()
runtime <- difftime(end, start, units = "mins")  # could be "mins", "hours", etc.
message("Runtime = ", round(runtime, 2), " minutes")

backward_trim <- backward %>% 
  select(config, block, treat, visit, kappa, everything()) %>% 
  group_by(config, block, treat, visit) %>%
  slice_min(Q_final) %>% 
  ungroup() %>% 
  as.data.table()
# 4. Source prediction 
sources_predicted <- backward_trim %>% 
  select(config, block, treat, visit, n_src, p_mat) %>% 
  pmap(~source_pred(config = ..1,
                    blk = ..2,
                    trt = ..3,
                    vst = ..4,
                    n_src = ..5,
                    p_mat = ..6,

                    mod_dat = mod_dat)) %>% 
  rbindlist()

results <- left_join(backward_trim, sources_predicted, by = c("config", "block", "treat", "visit", "n_src"))
saveRDS(results, here("DataProcessed/results/backward_model/backward_fits_sensitivity.rds"))