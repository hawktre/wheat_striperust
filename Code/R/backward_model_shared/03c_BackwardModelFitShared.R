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
source(here(
  "Code/R/backward_model_shared/03a_BackwardGradFunShared_vectorized.R"
))
source(here("Code/R/backward_model_shared/03b_BackwardModelFunShared.R"))
## Read in forward fits
forward <- readRDS(here("DataProcessed/results/forward_model/forward_fits.rds"))
mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat_arrays.rds"))

# Set up Indices ----------------------------------------------------------
blocks <- dimnames(mod_dat$intensity)[["block"]]
treats <- dimnames(mod_dat$intensity)[["treat"]]
visits <- dimnames(mod_dat$intensity)[["visit"]]
configs <- dimnames(mod_dat$groups)[["config"]]
n_src <- c(1, 2, 3, 4)

# Create full combinations
combos_backward <- expand.grid(
  config = configs,
  blk = blocks,
  trt = treats,
  n_src = n_src,
  vst = visits[-1],
  stringsAsFactors = FALSE
) |>
  filter(!(config == "4" & n_src == 4))
combos_backward <- left_join(
  combos_backward,
  forward %>% select(block, treat, visit, theta),
  by = c("blk" = "block", "trt" = "treat", "vst" = "visit")
)

# Get array task ID (which row to process)
# task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# if (is.na(task_id)) {
#   stop("SLURM_ARRAY_TASK_ID not set. This script must be run as a SLURM array job.")
# }

# Process only this task's row
combo <- combos_backward[956, ]
# message("Processing task ", task_id, " of ", nrow(combos_backward))
# message("Config: ", combo$config, ", Block: ", combo$blk, ", Treat: ", combo$trt,", N_SRC: ", combo$n_src, ", Visit: ", combo$vst)

start <- Sys.time()

# Fit model
backward_result <- backward_fit(
  config = combo$config,
  blk = combo$blk,
  trt = combo$trt,
  vst = combo$vst,
  n_src = combo$n_src,
  inits = combo$theta[[1]],
  mod_dat = mod_dat,
  tol = 1e-4,
  max_iter = 200
)

# Get predictions
if (backward_result$converged && !is.null(backward_result$p_mat)) {
  predictions <- source_pred(
    config = combo$config,
    blk = combo$blk,
    trt = combo$trt,
    vst = combo$vst,
    n_src = combo$n_src,
    p_mat = backward_result$p_mat[[1]],
    mod_dat = mod_dat
  )
} else {
  predictions <- data.table(
    config = combo$config,
    block = combo$blk,
    treat = combo$trt,
    visit = combo$vst,
    n_src = combo$n_src,
    mean_error = NA,
    n_correct = NA,
    acc = NA,
    dist_acc = NA,
    component_dist_acc = list(NA),
    predicted_source = list(NA),
    true_source = list(NA)
  )
}

# Drop large objects and merge
backward_result <- backward_result[, !c("p_mat", "pi")]
final_result <- merge(
  backward_result,
  predictions,
  by = c("config", "block", "treat", "visit", "n_src")
)

runtime <- difftime(Sys.time(), start, units = "mins")
message("Runtime = ", round(runtime, 2), " minutes")

# Save result
output_dir <- here("DataProcessed/results/backward_model/array_results_shared")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(
  final_result,
  file.path(
    output_dir,
    paste0(
      "backward_blk",
      combo$blk,
      "_trt",
      combo$trt,
      "_vst",
      combo$vst,
      "_config",
      combo$config,
      "_nsrc",
      combo$n_src,
      "_shared.rds"
    )
  )
)

message("Task ", task_id, " completed successfully")
