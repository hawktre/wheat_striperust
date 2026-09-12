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


combos_backward_t1 <- combos_backward |> 
  filter(trt == 1) |> 
  arrange(config, blk, trt, vst, kappa) |> 
  group_by(config, blk, trt, vst) |> 
  mutate(batch = cur_group_id()) |> 
  ungroup()
# Fit the model -----------------------------------------------------------

# Get array task ID (which row to process)
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

if (is.na(task_id)) {
  stop("SLURM_ARRAY_TASK_ID not set. This script must be run as a SLURM array job.")
}

# Process only this task's row
combos <- combos_backward_t1 |> filter(batch == task_id)

message("Processing task ", task_id, " of ", max(combos_backward_t1$batch))
message("Config: ", unique(combos$config), ", Block: ", unique(combos$blk), ", Treat: ", unique(combos$trt), ", Visit: ", unique(combos$vst))

start <- Sys.time()
backward_result <- lapply(seq_len(nrow(combos)), function(i)
  backward_fit(config = combos$config[i],
                                blk = combos$blk[i],
                                trt = combos$trt[i],
                                vst = combos$vst[i],
                                inits = combos$init[i][[1]],
                                mod_dat = mod_dat,
                                tol = 1e-4,
                                max_iter = 200)) |> rbindlist()

backward_result$kappa <- unique(combos$kappa)
end <- Sys.time()
runtime <- difftime(end, start, units = "mins")
message("Runtime = ", round(runtime, 2), " minutes")
backward_result$runtime <- runtime

# Save individual result
output_dir <- here("DataProcessed/results/backward_model/array_results_sensitivity")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(backward_result, file.path(output_dir, paste0("backward_blk", unique(combos$blk),"_trt",unique(combos$trt),"_vst",unique(combos$vst),"_config",unique(combos$config),"_sensitivity.rds")))

message("Task ", task_id, " completed successfully")
  