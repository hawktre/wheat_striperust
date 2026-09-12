## ---------------------------
##
## Script name: 04b_RunSim.R
##
## Purpose of script: Run the simulation
##
## Author: Trent VanHawkins
##
## Date Created: 2025-07-23
##
##
## ---------------------------

library(here)
library(data.table)
library(dplyr)
library(purrr)
library(parallel)

source(here("Code/R/forward_model/02b_ForwardModelFun.R"))
source(here("Code/R/backward_model_shared/03a_BackwardGradFunShared.R"))
source(here("Code/R/backward_model_shared/03b_BackwardModelFunShared.R"))
source(here("Code/R/simulation/04a_SimFunc.R"))

# Read in data
forward_fits <- readRDS(here(
  "DataProcessed/results/forward_model/forward_fits.rds"
))
mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat_arrays.rds"))
kappa_try <- exp(seq(log(0.5), log(4.0), length.out = 8))

# Take args from command line
args <- commandArgs(trailingOnly = TRUE)

nsim <- 100
if (length(args) >= 1) {
  nsim <- as.numeric(args[1])
}

# Get array task ID
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(task_id)) {
  stop("SLURM_ARRAY_TASK_ID not set.")
}
cat("Running", nsim, "simulations for task", task_id, "\n")

# Use SLURM_CPUS_PER_TASK if available
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(ncores) || ncores <= 0) {
  ncores <- parallel::detectCores(logical = FALSE)
}
ncores <- min(nsim, ncores)
cat("Using", ncores, "cores\n")

# Run the simulations
sim_list <- mclapply(
  seq_len(nsim),
  function(i) {
    t0 <- Sys.time()
    cur.sim <- (task_id * nsim) + i

    result <- single_sim(
      sim_id = cur.sim,
      dat = mod_dat,
      forward_mod = forward_fits,
      kappa_try = kappa_try,
      output_dir = here("DataProcessed/results/simulation/errors")
    )

    result <- result |>
      group_by(config, treat, visit) |>
      mutate(bic_selected = bic == min(bic)) |>
      filter(
        bic_selected |
          lengths(predicted_source) == lengths(true_source)
      ) |>
      ungroup()

    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(
      sprintf(
        "Sim %05d done in %.1f mins at %s\n",
        cur.sim,
        elapsed,
        format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      ),
      file = file.path(
        here("DataProcessed/results/simulation/logs"),
        "sim_progress.log"
      ),
      append = TRUE
    )

    result
  },
  mc.cores = ncores,
  mc.preschedule = FALSE
)

# Combine and save
sims <- rbindlist(sim_list)

output_dir <- here("DataProcessed/results/simulation/batch_results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(
  sims,
  file.path(output_dir, paste0("simulation_batch", task_id, ".rds"))
)

message("Task ", task_id, " completed successfully")
