## ---------------------------
##
## Script name: 04b_SimTest.R
##
## Purpose of script: Test simulation on a few random seeds locally
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
library(purrr)
library(parallel)
library(data.table)

## Read in necessary functions
source(here("Code/R/forward_model/02b_ForwardModelFun.R"))
source(here("Code/R/backward_model_shared/03b_BackwardModelFunShared.R"))
source(here("Code/R/simulation/04a_SimFunc.R"))

## Read in data
forward_mod <- readRDS(here(
  "DataProcessed/results/forward_model/forward_fits.rds"
))
mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat_arrays.rds"))
kappa_try <- exp(seq(log(0.5), log(4.0), length.out = 8))

# Test seeds --------------------------------------------------------------
test_seeds <- sample(1:1000, 1)
message("Testing seeds: ", paste(test_seeds, collapse = ", "))

start <- Sys.time()

results <- mclapply(
  test_seeds,
  function(sid) {
    message("Running sim_id = ", sid)
    single_sim(
      sim_id = sid,
      dat = mod_dat,
      forward_mod = forward_mod,
      kappa_try = kappa_try
    )
  },
  mc.cores = parallel::detectCores() - 1
)

runtime <- difftime(Sys.time(), start, units = "mins")
message("Total runtime: ", round(runtime, 2), " minutes")

# Report
failed <- which(sapply(results, is.null))
if (length(failed) > 0) {
  message("Failed seeds: ", paste(test_seeds[failed], collapse = ", "))
} else {
  message("All seeds completed successfully.")
}

final <- rbindlist(results[!sapply(results, is.null)])
print(final)
