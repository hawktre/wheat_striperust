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
library(parallel)  # Use base parallel package
source(here("Code/R/forward_model/02b_ForwardModelFun.R"))
source(here("Code/R/backward_model_shared/03b_BackwardModelFunSimulation.R"))
source(here("Code/R/simulation/04a_SimFunc.R"))

# Read in data
forward_fits <- readRDS(here("DataProcessed/results/forward_model/forward_fits.rds"))
mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat_arrays.rds"))
kappa_try <- c(0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0)

# Take args from command line
args <- commandArgs(trailingOnly = TRUE)

# Set default value of simulations
nsim <- 100
# Override with argument if provided
if (length(args) >= 1) {
  nsim <- as.numeric(args[1])
}

# Get array task ID (which row to process)
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

cat("Running", nsim, "simulations\n")

# Use SLURM_CPUS_PER_TASK if available
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))

if (is.na(ncores) || ncores <= 0) {
  ncores <- parallel::detectCores(logical = FALSE)
}
ncores <- min(nsim, ncores)
cat("Using", ncores, "cores for simulations\n")

sim_list <- mclapply(seq_len(nsim), function(i) {
  t0 <- Sys.time()
  cur.sim <- (task_id * nsim) + i
  result <- single_sim(cur.sim, mod_dat, forward_fits, kappa_try = kappa_try, output_dir = here("DataProcessed/results/simulation/errors"))

  t1 <- Sys.time()
  elapsed <- as.numeric(difftime(t1, t0, units = "mins"))

  log_msg <- sprintf("Sim %05d done in %.1f mins at %s\n", cur.sim, elapsed, format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  cat(log_msg, file = file.path(here("DataProcessed/results/simulation/logs"), paste0("sim_progress.log")), append = TRUE)

  result
}, mc.cores = ncores, mc.preschedule = FALSE)


# Combine results into a data.frame
sims <- rbindlist(sim_list)

# Save individual result
output_dir <- here("DataProcessed/results/simulation/batch_results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(sims, file.path(output_dir, paste0("simulation_batch", task_id,".rds")))

message("Task ", task_id, " completed successfully")