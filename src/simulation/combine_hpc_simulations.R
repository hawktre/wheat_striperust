# Combine only complete whole-simulation HPC result files.

library(dplyr)
library(here)
library(readr)

arguments <- commandArgs(trailingOnly = TRUE)
expected_simulations <- if (length(arguments)) as.integer(arguments[1]) else NA_integer_

study_directory <- here(
  "output", "simulation", "whole_epidemic_study", "hpc"
)
complete_directory <- file.path(study_directory, "complete")
failed_directory <- file.path(study_directory, "failed")
combined_directory <- file.path(study_directory, "combined")
dir.create(combined_directory, recursive = TRUE, showWarnings = FALSE)

files <- sort(list.files(
  complete_directory, pattern = "^simulation_[0-9]+\\.rds$", full.names = TRUE
))
if (!length(files)) stop("No complete HPC simulation files were found.")

objects <- lapply(files, readRDS)
valid <- vapply(objects, function(x) {
  isTRUE(x$status$simulation_complete) && !is.null(x$results)
}, logical(1))
if (!all(valid)) stop("At least one file in the complete directory is invalid.")

simulation_ids <- vapply(objects, function(x) x$status$simulation, integer(1))
if (anyDuplicated(simulation_ids)) stop("Duplicate simulation IDs were found.")
if (!is.na(expected_simulations) && length(objects) != expected_simulations) {
  stop(
    "Found ", length(objects), " complete simulations; expected ",
    expected_simulations, "."
  )
}

status <- bind_rows(lapply(objects, `[[`, "status")) |>
  arrange(simulation)
component_names <- names(objects[[1]]$results)
combined <- setNames(lapply(component_names, function(component) {
  bind_rows(lapply(objects, function(x) x$results[[component]]))
}), component_names)

write_csv(status, file.path(combined_directory, "simulation_status.csv"))

failed_files <- sort(list.files(
  failed_directory, pattern = "^simulation_[0-9]+\\.rds$", full.names = TRUE
))
if (length(failed_files)) {
  failed_status <- bind_rows(lapply(failed_files, function(path) {
    readRDS(path)$status
  })) |>
    arrange(simulation)
} else {
  failed_status <- tibble()
}
attempt_status <- status |> mutate(result_set = "included")
if (nrow(failed_status)) {
  attempt_status <- bind_rows(
    attempt_status,
    failed_status |> mutate(result_set = "failed")
  )
}
write_csv(
  attempt_status,
  file.path(combined_directory, "all_attempt_status.csv")
)

for (component in names(combined)) {
  write_csv(
    combined[[component]],
    file.path(combined_directory, paste0(component, ".csv"))
  )
}
saveRDS(
  list(status = status, results = combined),
  file.path(combined_directory, "combined_simulations.rds")
)
message("Combined ", length(objects), " complete simulation replicates")
