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

# Create full combinations
combos_backward <- expand.grid(config = configs, blk = blocks, trt = treats, vst = visits[-1], stringsAsFactors = FALSE)
combos_backward <- left_join(combos_backward, forward %>% select(block, treat, visit, theta), 
                             by = c("blk"="block", "trt"="treat", "vst"="visit")) |> 
  filter(!(config == "64" & trt == 4))

# Get array task ID (which row to process)
#task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
task_id <- 2
if (is.na(task_id)) {
  stop("SLURM_ARRAY_TASK_ID not set. This script must be run as a SLURM array job.")
}

# Process only this task's row
combo <- combos_backward[task_id, ]

message("Processing task ", task_id, " of ", nrow(combos_backward))
message("Config: ", combo$config, ", Block: ", combo$blk, ", Treat: ", combo$trt, ", Visit: ", combo$vst)

start <- Sys.time()
backward_result <- backward_fit(config = combo$config,
                                blk = combo$blk,
                                trt = combo$trt,
                                vst = combo$vst,
                                inits = combo$theta[[1]],
                                mod_dat = mod_dat,
                                tol = 1e-4,
                                max_iter = 1000)
end <- Sys.time()
runtime <- difftime(end, start, units = "mins")
message("Runtime = ", round(runtime, 2), " minutes")

# Save individual result
output_dir <- here("DataProcessed/results/backward_model/array_results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(backward_result, file.path(output_dir, paste0("backward_", task_id, ".rds")))

message("Task ", task_id, " completed successfully")