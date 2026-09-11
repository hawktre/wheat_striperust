# Run a small local pilot of recursively simulated whole epidemics.

library(dplyr)
library(here)
library(readr)

source(here("Code", "R", "source_detection", "source_detection_functions.R"))
source(here("Code", "R", "simulation", "whole_epidemic_pilot_functions.R"))

mod_dat <- readRDS(here("DataProcessed", "experimental", "mod_dat_arrays.rds"))
forward_fits <- readRDS(here(
  "DataProcessed", "results", "source_detection_main", "forward_fits.rds"
))
stripe <- readRDS(here("DataProcessed", "experimental", "stripe_clean.rds"))
coordinates <- stripe |>
  select(plant_id, east, north) |>
  distinct() |>
  arrange(plant_id)

# Kept intentionally small enough for local experimentation.
pilot_blocks <- c("A", "B")
pilot_treatments <- c("1", "2", "4")
pilot_replicates <- 1:4
mechanisms <- c("secondary_spread", "persistent_sources")

epidemics <- list()
epidemic_summaries <- list()
plant_summaries <- list()
oracle_diagnostics <- list()
index <- 0L

for (block in pilot_blocks) {
  for (treat in pilot_treatments) {
    forward_fit <- forward_fits[[paste(block, treat, sep = "_")]]
    configs <- c("4", "8h", "8v", "16")
    if (treat == "1") configs <- c(configs, "64")

    observed <- mod_dat$intensity[, block, treat, ]
    source_plants <- get_true_source_groups(mod_dat, block, treat, "64")
    epidemic_summaries[[length(epidemic_summaries) + 1L]] <- epidemic_summary(
      observed, block, treat, 0L, "observed", "observed"
    )
    plant_summaries[[length(plant_summaries) + 1L]] <- plant_level_summary(
      observed, coordinates, source_plants, block, treat, 0L, "observed", "observed"
    )

    for (mechanism in mechanisms) {
      for (replicate in pilot_replicates) {
        index <- index + 1L
        seed <- 10000L + index
        message(
          "Simulating block ", block, ", treatment ", treat,
          ", mechanism ", mechanism, ", replicate ", replicate
        )
        simulated <- simulate_whole_epidemic(
          mod_dat, forward_fit, block, treat, mechanism, seed
        )
        simulated_mod_dat <- replace_experiment_intensity(
          mod_dat, block, treat, simulated
        )
        epidemic_summaries[[length(epidemic_summaries) + 1L]] <- epidemic_summary(
          simulated, block, treat, replicate, mechanism, "simulated"
        )
        plant_summaries[[length(plant_summaries) + 1L]] <- plant_level_summary(
          simulated, coordinates, source_plants, block, treat, replicate,
          mechanism, "simulated"
        )
        for (config in configs) {
          oracle_diagnostics[[length(oracle_diagnostics) + 1L]] <-
            oracle_source_diagnostics(
              simulated_mod_dat, forward_fit, block, treat, config,
              replicate, mechanism
            )
        }
        epidemics[[paste(block, treat, mechanism, replicate, sep = "_")]] <- simulated
      }
    }
  }
}

output_directory <- here(
  "DataProcessed", "results", "simulation", "whole_epidemic_pilot"
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
write_csv(bind_rows(epidemic_summaries), file.path(output_directory, "epidemic_summary.csv"))
write_csv(bind_rows(plant_summaries), file.path(output_directory, "plant_severity.csv"))
write_csv(bind_rows(oracle_diagnostics), file.path(output_directory, "oracle_source_diagnostics.csv"))
saveRDS(epidemics, file.path(output_directory, "simulated_epidemics.rds"))

source(here("Code", "R", "simulation", "plot_whole_epidemic_pilot.R"))
message("Whole-epidemic pilot complete: ", output_directory)
