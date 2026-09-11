# Readable, sequential reference version of the primary analysis pipeline.
#
# This script is documentation and quality-control code. The production
# analysis remains Code/R/source_detection/run_main_analysis.R.

library(dplyr)
library(here)
library(readr)
library(RcppHungarian)

source(here(
  "Code", "R", "source_detection", "readable_reference",
  "source_detection_functions_readable.R"
))

model_data <- readRDS(here(
  "DataProcessed", "experimental", "mod_dat_arrays.rds"
))

blocks <- dimnames(model_data$intensity)[["block"]]
treatments <- dimnames(model_data$intensity)[["treat"]]
study_visits <- dimnames(model_data$intensity)[["visit"]]
transition_visits <- study_visits[-1L]
available_resolutions <- dimnames(model_data$groups)[["config"]]

resolutions_for_single_source_treatment <- intersect(
  c("4", "8h", "8v", "16", "64"),
  available_resolutions
)
resolutions_for_multiple_source_treatments <- intersect(
  c("4", "8h", "8v", "16"),
  available_resolutions
)
initial_kappa_grid <- exp(seq(log(0.5), log(2), length.out = 5))

# Fit one full-forward initialization for each observed experiment.
forward_results <- list()
for (block in blocks) {
  for (treatment in treatments) {
    experiment_name <- paste(block, treatment, sep = "_")
    forward_results[[experiment_name]] <- fit_forward_model_readable(
      block = block,
      treatment = treatment,
      model_data = model_data,
      initial_kappa_grid = initial_kappa_grid
    )
  }
}

# Fit the known source count under every primary temporal scope.
all_em_results <- list()
result_number <- 0L

for (block in blocks) {
  for (treatment in treatments) {
    if (as.integer(treatment) == 1L) {
      resolutions <- resolutions_for_single_source_treatment
    } else {
      resolutions <- resolutions_for_multiple_source_treatments
    }

    experiment_name <- paste(block, treatment, sep = "_")
    initial_theta <- forward_results[[experiment_name]]

    for (resolution in resolutions) {
      number_of_sources <- length(get_true_sources_readable(
        model_data, block, treatment, resolution
      ))

      for (endpoint_visit in transition_visits) {
          # Independent transition fit.
          result_number <- result_number + 1L
          all_em_results[[result_number]] <- fit_source_detection_readable(
            block = block,
            treatment = treatment,
            resolution = resolution,
            number_of_sources = number_of_sources,
            model_data = model_data,
            initial_theta = initial_theta,
            single_transition = endpoint_visit
          )

          # Cumulative fit. The final endpoint is the whole-study fit.
          result_number <- result_number + 1L
          all_em_results[[result_number]] <- fit_source_detection_readable(
            block = block,
            treatment = treatment,
            resolution = resolution,
            number_of_sources = number_of_sources,
            model_data = model_data,
            initial_theta = initial_theta,
            through_visit = endpoint_visit
          )
      }
    }
  }
}

# Stack the readable outputs.
candidate_results <- data.frame()
parameter_results <- data.frame()
accuracy_results <- data.frame()

for (fit in all_em_results) {
  candidate_results <- bind_rows(
    candidate_results,
    fit$candidate_summary
  )
  parameter_results <- bind_rows(
    parameter_results,
    fit$parameter_summary
  )
  accuracy_results <- bind_rows(
    accuracy_results,
    fit$accuracy_summary
  )
}

# These files are intentionally separate from the production output folder.
reference_output_directory <- here(
  "DataProcessed", "results", "source_detection_readable_reference"
)
dir.create(reference_output_directory, recursive = TRUE, showWarnings = FALSE)

write_csv(
  candidate_results,
  file.path(reference_output_directory, "candidate_posteriors.csv")
)
write_csv(
  parameter_results,
  file.path(reference_output_directory, "parameter_estimates.csv")
)
write_csv(
  accuracy_results,
  file.path(reference_output_directory, "source_accuracy.csv")
)

message(
  "Readable reference analysis complete: ", reference_output_directory
)
