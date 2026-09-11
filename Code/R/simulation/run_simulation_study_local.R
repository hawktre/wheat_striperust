# Small, checkpointed local runner for the whole-epidemic simulation study.

library(dplyr)
library(here)
library(parallel)
library(readr)

source(here("Code", "R", "simulation", "simulation_study_functions.R"))

mod_dat <- readRDS(here("DataProcessed", "experimental", "mod_dat_arrays.rds"))
generating_fits <- readRDS(here(
  "DataProcessed", "results", "source_detection_main", "forward_fits.rds"
))

local_simulations <- 5L
blocks <- dimnames(mod_dat$intensity)[["block"]]
treatments <- dimnames(mod_dat$intensity)[["treat"]]
scenario_grid <- expand.grid(
  simulation = seq_len(local_simulations),
  block = blocks,
  treat = treatments,
  stringsAsFactors = FALSE
)

output_directory <- here(
  "DataProcessed", "results", "simulation", "whole_epidemic_study", "local"
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
checkpoint_file <- file.path(output_directory, "local_checkpoint.rds")

requested_cores <- suppressWarnings(as.integer(Sys.getenv(
  "SOURCE_DETECTION_CORES", unset = "4"
)))
n_cores <- max(1L, min(requested_cores, length(blocks) * length(treatments)))

if (file.exists(checkpoint_file)) {
  results <- readRDS(checkpoint_file)
  if (length(results) != nrow(scenario_grid)) {
    results <- vector("list", nrow(scenario_grid))
  }
} else {
  results <- vector("list", nrow(scenario_grid))
}

for (simulation_id in seq_len(local_simulations)) {
  indices <- which(scenario_grid$simulation == simulation_id)
  indices <- indices[vapply(results[indices], is.null, logical(1))]
  if (!length(indices)) next

  message(
    "Running simulation replicate ", simulation_id,
    " (", length(indices), " scenarios) on ", n_cores, " cores"
  )
  batch_results <- parallel::mclapply(indices, function(i) {
    spec <- scenario_grid[i, ]
    message(
      "Simulation ", spec$simulation, ", block ", spec$block,
      ", treatment ", spec$treat
    )
    run_simulation_scenario(
      simulation = spec$simulation,
      block = spec$block,
      treat = spec$treat,
      mod_dat = mod_dat,
      generating_fit = generating_fits[[paste(
        spec$block, spec$treat, sep = "_"
      )]],
      save_fits = FALSE
    )
  }, mc.cores = n_cores, mc.preschedule = FALSE)
  results[indices] <- batch_results
  saveRDS(results, checkpoint_file)
}

combined <- combine_simulation_results(results)
for (name in names(combined)) {
  write_csv(combined[[name]], file.path(output_directory, paste0(name, ".csv")))
}
message("Local simulation study complete: ", output_directory)
