## ---------------------------
## Script name: 04b_RunSim_debug.R
## Purpose: Debug individual failing simulation seeds locally
## ---------------------------

library(here)
library(data.table)
library(dplyr)
library(purrr)

source(here("Code/R/forward_model/02b_ForwardModelFun.R"))
source(here("Code/R/backward_model_shared/03b_BackwardModelFunShared.R"))
source(here("Code/R/simulation/04a_SimFunc.R"))

# Read in data
forward_fits <- readRDS(here(
  "DataProcessed/results/forward_model/forward_fits.rds"
))
mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat_arrays.rds"))
kappa_try <- c(0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0)

# --- Set your seeds to debug here ---
seeds_to_debug <- c(1)

sim_list <- vector("list", length(seeds_to_debug))

for (i in seq_along(seeds_to_debug)) {
  seed <- seeds_to_debug[i]
  cat("Running seed", seed, "...\n")

  result <- tryCatch(
    {
      single_sim(
        seed,
        mod_dat,
        forward_fits,
        kappa_try = kappa_try,
        output_dir = here("DataProcessed/results/simulation/errors")
      )
    },
    error = function(e) {
      cat("ERROR on seed", seed, ":", conditionMessage(e), "\n")
      NULL
    }
  )

  if (!is.null(result)) {
    result <- result |>
      group_by(config, treat, visit) |>
      mutate(bic_selected = bic == min(bic)) |>
      filter(
        bic_selected | lengths(predicted_source) == lengths(true_source)
      ) |>
      ungroup()

    cat("Seed", seed, "completed —", nrow(result), "rows\n")
  }

  sim_list[[i]] <- result
}

# Check each element before binding
for (i in seq_along(sim_list)) {
  cat(
    "Seed",
    seeds_to_debug[i],
    "-> class:",
    class(sim_list[[i]]),
    "| is NULL:",
    is.null(sim_list[[i]]),
    "\n"
  )
  if (!is.null(sim_list[[i]])) {
    cat("  visit class:", class(sim_list[[i]]$visit), "\n")
  }
}

# Attempt rbindlist
sims <- tryCatch(
  {
    rbindlist(sim_list)
  },
  error = function(e) {
    cat("rbindlist ERROR:", conditionMessage(e), "\n")
    NULL
  }
)
