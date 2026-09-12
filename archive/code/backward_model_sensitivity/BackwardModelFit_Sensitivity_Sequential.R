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
library(dplyr)
library(purrr)
library(data.table)

## Read in necessary functions
source(here("Code/R/backward_model_shared/03b_BackwardModelFunShared.R"))

## Read in model data
mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat_arrays.rds"))

# Set up Indices ----------------------------------------------------------
blocks <- dimnames(mod_dat$intensity)[["block"]]
treats <- dimnames(mod_dat$intensity)[["treat"]]
visits <- dimnames(mod_dat$intensity)[["visit"]]
configs <- dimnames(mod_dat$groups)[["config"]]
kappa_try <- c(0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0)


# Get the inits for each value of kappa -----------------------------------
combos_backward <- expand.grid(config = configs, blk = blocks, trt = treats, vst = visits[-1], kappa = kappa_try, stringsAsFactors = FALSE) |>
  filter(!(config == "64" & trt == "4"))

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
backward_results <- lapply(seq_len(nrow(combos_backward)), function(row) {
  combo <- combos_backward[row,]
  backward_fit(config = combo$config,
                                blk = combo$blk,
                                trt = combo$trt,
                                vst = combo$vst,
                                inits = combo$init[[1]],
                                mod_dat = mod_dat,
                                tol = 1e-4,
                                max_iter = 200)
}) |> rbindlist()
end <- Sys.time()
runtime <- difftime(end, start, units = "mins")
message("Runtime = ", round(runtime, 2), " minutes")

saveRDS(backward_results, here("DataProcessed/results/backward_model/backward_fits_sensitivity_sequential.rds"))