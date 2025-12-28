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

source(here("Code/R/backward_model_shared/03b_BackwardModelFunShared.R"))

# Read all individual results
array_dir <- here("DataProcessed/results/simulation/batch_results/")
result_files <- list.files(array_dir, pattern = "^simulation_.*\\.rds$", full.names = TRUE)

message("Found ", length(result_files), " result files")

sims <- map(result_files, readRDS) %>% rbindlist()

saveRDS(sims, here("DataProcessed/results/simulation/sims.rds"))

message("Successfully combined all results")