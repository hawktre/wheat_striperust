# Run one complete Monte Carlo replicate on an HPC array task.
#
# One replicate contains every block x treatment scenario. It is marked usable
# only if all required forward and backward fits complete successfully.

library(dplyr)
library(here)
library(parallel)

source(here("src", "simulation", "simulation_study_functions.R"))

arguments <- commandArgs(trailingOnly = TRUE)
simulation <- if (length(arguments)) {
  as.integer(arguments[1])
} else {
  as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA_character_))
}
if (is.na(simulation) || simulation < 1L) {
  stop("A positive simulation ID is required.")
}

requested_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
if (is.na(requested_cores) || requested_cores < 1L) requested_cores <- 1L

mod_dat <- readRDS(here("data", "processed", "experimental", "mod_dat_arrays.rds"))
generating_fits <- readRDS(here(
  "output", "source_detection", "forward_fits.rds"
))
scenario_grid <- expand.grid(
  block = dimnames(mod_dat$intensity)[["block"]],
  treat = dimnames(mod_dat$intensity)[["treat"]],
  stringsAsFactors = FALSE
)
n_cores <- min(requested_cores, nrow(scenario_grid))

transition_visits <- dimnames(mod_dat$intensity)[["visit"]][-1L]
available_configs <- dimnames(mod_dat$groups)[["config"]]
expected_forward_transitions <- nrow(scenario_grid) * length(transition_visits)
expected_backward_fits <- sum(vapply(
  scenario_grid$treat,
  function(treat) {
    2L * length(transition_visits) *
      length(simulation_configs(treat, available_configs))
  },
  integer(1)
))

study_directory <- here(
  "output", "simulation", "whole_epidemic_study", "hpc"
)
complete_directory <- file.path(study_directory, "complete")
failed_directory <- file.path(study_directory, "failed")
dir.create(complete_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(failed_directory, recursive = TRUE, showWarnings = FALSE)
file_stem <- sprintf("simulation_%06d", simulation)
complete_file <- file.path(complete_directory, paste0(file_stem, ".rds"))
failed_file <- file.path(failed_directory, paste0(file_stem, ".rds"))

if (file.exists(complete_file)) {
  message("Simulation ", simulation, " already has a complete result; skipping")
  quit(save = "no", status = 0L)
}

message(
  "Running complete simulation ", simulation, " with ", n_cores,
  " block-treatment workers"
)
started <- Sys.time()
scenario_results <- parallel::mclapply(
  seq_len(nrow(scenario_grid)),
  function(i) {
    spec <- scenario_grid[i, ]
    tryCatch(
      run_simulation_scenario(
        simulation = simulation,
        block = spec$block,
        treat = spec$treat,
        mod_dat = mod_dat,
        generating_fit = generating_fits[[paste(
          spec$block, spec$treat, sep = "_"
        )]],
        save_fits = FALSE
      ),
      error = function(e) {
        list(
          scenario_status = data.frame(
            simulation = simulation,
            block = spec$block,
            treat = as.character(spec$treat),
            forward_success = FALSE,
            forward_elapsed_seconds = NA_real_,
            error = paste("Unhandled scenario error:", conditionMessage(e))
          )
        )
      }
    )
  },
  mc.cores = n_cores,
  mc.preschedule = FALSE
)

# A worker-level crash can be returned by mclapply as a try-error rather than
# passing through the tryCatch inside the worker. Convert that case to the same
# scenario-status structure used for ordinary fitting errors.
for (i in seq_along(scenario_results)) {
  result <- scenario_results[[i]]
  if (is.list(result) && !is.null(result$scenario_status)) next

  spec <- scenario_grid[i, ]
  scenario_results[[i]] <- list(
    scenario_status = data.frame(
      simulation = simulation,
      block = spec$block,
      treat = as.character(spec$treat),
      forward_success = FALSE,
      forward_elapsed_seconds = NA_real_,
      error = paste("Worker failure:", as.character(result)[1])
    )
  )
}

combined <- combine_simulation_results(scenario_results)
all_scenarios_returned <- nrow(combined$scenario_status) == nrow(scenario_grid)
all_forward_fits_succeeded <- all_scenarios_returned &&
  isTRUE(all(combined$scenario_status$forward_success))
all_forward_optimizers_converged <- all_forward_fits_succeeded &&
  nrow(combined$forward_metrics) == expected_forward_transitions &&
  isTRUE(all(combined$forward_metrics$convergence == 0L))
all_backward_fits_returned <- all_forward_fits_succeeded &&
  nrow(combined$backward_fit_status) == expected_backward_fits
all_backward_fits_succeeded <- all_backward_fits_returned &&
  isTRUE(all(combined$backward_fit_status$success))
all_backward_fits_converged <- all_backward_fits_succeeded &&
  nrow(combined$backward_diagnostics) == expected_backward_fits &&
  isTRUE(all(combined$backward_diagnostics$converged)) &&
  isTRUE(all(combined$backward_diagnostics$all_iterations_monotonic))
any_transient_m_step_failure <- all_backward_fits_succeeded &&
  isTRUE(any(combined$backward_diagnostics$any_m_step_failure))

simulation_complete <- all(c(
  all_scenarios_returned,
  all_forward_fits_succeeded,
  all_forward_optimizers_converged,
  all_backward_fits_returned,
  all_backward_fits_succeeded,
  all_backward_fits_converged
))
failure_reasons <- names(which(!c(
  all_scenarios_returned = all_scenarios_returned,
  all_forward_fits_succeeded = all_forward_fits_succeeded,
  all_forward_optimizers_converged = all_forward_optimizers_converged,
  all_backward_fits_returned = all_backward_fits_returned,
  all_backward_fits_succeeded = all_backward_fits_succeeded,
  all_backward_fits_converged = all_backward_fits_converged
)))
status <- data.frame(
  simulation = simulation,
  simulation_complete = simulation_complete,
  all_scenarios_returned = all_scenarios_returned,
  all_forward_fits_succeeded = all_forward_fits_succeeded,
  all_forward_optimizers_converged = all_forward_optimizers_converged,
  all_backward_fits_returned = all_backward_fits_returned,
  all_backward_fits_succeeded = all_backward_fits_succeeded,
  all_backward_fits_converged = all_backward_fits_converged,
  any_transient_m_step_failure = any_transient_m_step_failure,
  n_scenarios = nrow(combined$scenario_status),
  n_forward_transitions = nrow(combined$forward_metrics),
  n_backward_fits = nrow(combined$backward_fit_status),
  expected_forward_transitions = expected_forward_transitions,
  expected_backward_fits = expected_backward_fits,
  failure_reasons = if (length(failure_reasons)) {
    paste(failure_reasons, collapse = "; ")
  } else {
    NA_character_
  },
  started = format(started, "%Y-%m-%d %H:%M:%S %Z"),
  finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  elapsed_minutes = as.numeric(difftime(Sys.time(), started, units = "mins")),
  slurm_job_id = Sys.getenv("SLURM_JOB_ID", unset = NA_character_),
  slurm_array_task_id = Sys.getenv(
    "SLURM_ARRAY_TASK_ID", unset = NA_character_
  ),
  cores = n_cores
)
output <- list(status = status, results = combined)

destination <- if (simulation_complete) complete_file else failed_file
temporary_file <- paste0(destination, ".tmp_", Sys.getpid())
saveRDS(output, temporary_file)
if (!file.rename(temporary_file, destination)) {
  stop("Could not atomically move the completed result to ", destination)
}

if (simulation_complete) {
  if (file.exists(failed_file)) file.remove(failed_file)
  message(
    "Simulation ", simulation, " completed in ",
    round(status$elapsed_minutes, 2), " minutes"
  )
  quit(save = "no", status = 0L)
}

message("Simulation ", simulation, " failed completeness requirements")
quit(save = "no", status = 2L)
