## ---------------------------
##
## Script name: 03d_CombineBackwardResults.R
##
## Purpose of script: Combine backward model array results
##
## ---------------------------

library(dplyr)
library(purrr)
library(here)
library(data.table)

source(here("Code/R/backward_model/03b_BackwardModelFun.R"))

forward <- readRDS(here("DataProcessed/results/forward_model/forward_fits.rds"))
mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat_arrays.rds"))

# Read all individual results
array_dir <- here("DataProcessed/results/backward_model/array_results_shared")
result_files <- list.files(array_dir, pattern = "^backward_.*\\.rds$", full.names = TRUE)

message("Found ", length(result_files), " result files")

backward <- map(result_files, readRDS) %>% rbindlist()

saveRDS(backward, here("DataProcessed/results/backward_model/backward_fits_shared.rds"))

message("Successfully combined all results")