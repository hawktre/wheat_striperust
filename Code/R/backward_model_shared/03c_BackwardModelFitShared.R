## ---------------------------
##
## Script name: 03c_BackwardModelFit.R
##
## Purpose of script: Fit the backward source-prediction model (array version)
##
## Author: Trent VanHawkins
##
## Date Created: 2025-08-20
##
##
## ---------------------------

options(scipen = 6, digits = 4)

library(dplyr)
library(here)
library(parallel)
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

# Create full combinations
combos_backward <- expand.grid(config = configs, blk = blocks, trt = treats, vst = visits[-1], stringsAsFactors = FALSE)
combos_backward <- left_join(combos_backward, forward %>% select(block, treat, visit, theta), 
                             by = c("blk"="block", "trt"="treat", "vst"="visit")) 



start <- Sys.time()
backward <- mclapply(seq_len(nrow(combos_backward)), function(row) {
  combo <- combos_backward[row,]
  backward_fit(combo$config, combo$blk, combo$trt, combo$vst, mod_dat, combo$theta[[1]], max_iter = 200, tol = 1e-8)
}, mc.cores = detectCores()-1, mc.preschedule = F) |> rbindlist()
end <- Sys.time()
runtime <- difftime(end, start, units = "mins")
message("Runtime = ", round(runtime, 2), " minutes")

# Save individual result
saveRDS(backward, here("DataProcessed/results/backward_model/backward_fits_shared.rds"))