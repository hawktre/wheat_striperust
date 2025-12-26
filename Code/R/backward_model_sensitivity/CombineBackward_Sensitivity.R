## ---------------------------
##
## Script name: 03d_CombineBackwardResults.R
##
## Purpose of script: Combine backward model array results
##
## ---------------------------

library(here)
library(purrr)
library(dplyr)
library(data.table)


# Read all individual results
array_dir <- here("DataProcessed/results/backward_model/array_results_sensitivity")
result_files <- list.files(array_dir, pattern = "^backward_.*\\.rds$", full.names = TRUE)

message("Found ", length(result_files), " result files")

backward <- rbindlist(map(result_files, readRDS))

saveRDS(backward, here("DataProcessed/results/backward_model/backward_fits_sensitivity.rds"))

message("Successfully combined all results")