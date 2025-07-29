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

user_lib <- "~/R/x86_64-pc-linux-gnu-library/4.4"
.libPaths(user_lib)

library(here)
library(data.table)
library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(purrr)
library(parallel)  # Use base parallel package
source(here("Code/02a_GradDescentFun.R"))
source(here("Code/03a_EMgradfun.R"))
source(here("Code/04a_SimFunc.R"))

# Read in data
forward_fits <- readRDS(here("DataProcessed/results/forward_model/forward_fits.rds"))
mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat_arrays.rds"))

# Take args from command line
args <- commandArgs(trailingOnly = TRUE)

# Set default value of simulations
nsim <- 20
# Override with argument if provided
if (length(args) >= 1) {
  nsim <- as.numeric(args[1])
}

cat("Running", nsim, "simulations\n")

ncores <- detectCores(logical = FALSE)  # Physical cores only (optional)
ncores <- min(nsim, ncores)

set.seed(404)

sim_list <- mclapply(1:nsim, function(i) {
  t0 <- Sys.time()

  result <- single_sim(i, mod_dat, forward_fits, output_dir = here("DataProcessed/results/simulation/errors"))

  t1 <- Sys.time()
  elapsed <- as.numeric(difftime(t1, t0, units = "mins"))

  log_msg <- sprintf("Sim %05d done in %.1f mins at %s\n", i, elapsed, format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  cat(log_msg, file = file.path(here("DataProcessed/results/simulation"), "sim_progress.log"), append = TRUE)

  result
}, mc.cores = ncores, mc.preschedule = FALSE)


# Combine results into a data.frame
sims <- rbindlist(sim_list)

# Save
saveRDS(sims, here("DataProcessed/results/simulation/sims.rds"))
