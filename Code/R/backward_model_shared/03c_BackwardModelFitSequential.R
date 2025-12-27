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
                             by = c("blk"="block", "trt"="treat", "vst"="visit")) |> 
  filter(!(config == "64" & trt == "4"))

start <- Sys.time()
backward_results <- lapply(seq_len(nrow(combos_backward)), function(row) {
  combo <- combos_backward[row,]
  message("Config: ", combo$config, ", Block: ", combo$blk, ", Treat: ", combo$trt, ", Visit: ", combo$vst)
  backward_fit(config = combo$config,
                                blk = combo$blk,
                                trt = combo$trt,
                                vst = combo$vst,
                                inits = combo$theta[[1]],
                                mod_dat = mod_dat,
                                tol = 1e-4,
                                max_iter = 200)
}) |> rbindlist()
end <- Sys.time()
runtime <- difftime(end, start, units = "mins")
message("Runtime = ", round(runtime, 2), " minutes")