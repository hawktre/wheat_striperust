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

library(tidyverse)
library(furrr)
library(future)
library(here)
library(data.table)

## Read in necessary functions
source(here("Code/03b_BackwardModelFun.R"))

## Read in forward fits
forward <- readRDS(here("DataProcessed/results/forward_model/forward_fits.rds"))
mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat_arrays.rds"))

# Set up Indices ----------------------------------------------------------
blocks <- dimnames(mod_dat$intensity)[["block"]]
treats <- dimnames(mod_dat$intensity)[["treat"]]
visits <- dimnames(mod_dat$intensity)[["visit"]]
configs <- dimnames(mod_dat$groups)[["config"]]

# Fit the model -----------------------------------------------------------
# 2. Fit backward model
combos_backward <- expand.grid(config = configs, blk = blocks, trt = treats, vst = visits[-1], stringsAsFactors = FALSE)

start <- Sys.time()
backward <- pmap(combos_backward, ~backward_fit(..1, ..2, ..3, ..4, mod_dat = mod_dat, forward), .progress = T) %>% rbindlist()
end <- Sys.time()
runtime <- difftime(end, start, units = "mins")  # could be "mins", "hours", etc.
message("Runtime = ", round(runtime, 2), " minutes")

# 4. Source prediction for treat == 1
backward_t1 <- backward[treat == 1]
sources_predicted <- pmap(backward_t1, ~source_pred(config = ..1,
                                                    blk = ..2,
                                                    trt = ..3,
                                                    vst = ..4,
                                                    p_mat = ..10,
                                                    mod_dat = mod_dat)) %>% 
  rbindlist()

results_merge <- left_join(backward |> mutate(treat = as.factor(treat), 
visit = as.factor(visit)), forward |> mutate(treat = as.factor(treat), 
visit = as.factor(visit)), by = c("block", "treat", "visit"), suffix = c(".backward", ".forward")) %>% 
  left_join(sources_predicted |> mutate(treat = as.factor(treat), 
visit = as.factor(visit)), by = c("config", "block", "treat", "visit")) 

saveRDS(results_merge, here("DataProcessed/results/backward_model/backward_fits.rds"))

sources_predicted |> 
  group_by(config, visit) |> 
  summarise(mean_acc = mean(dist_acc))
