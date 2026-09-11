# Run the primary source-detection analysis.

library(dplyr)
library(here)
library(parallel)
library(purrr)
library(readr)

source(here("Code", "R", "source_detection", "source_detection_functions.R"))

mod_dat <- readRDS(here("DataProcessed", "experimental", "mod_dat_arrays.rds"))

# The known number of occupied source groups is used at each resolution.
# Individual plants (config 64) are retained for one-source inference;
# multi-source fits use at most the 16-cell resolution to control combinatorial
# growth.
blocks <- dimnames(mod_dat$intensity)[["block"]]
treatments <- dimnames(mod_dat$intensity)[["treat"]]
available_configs <- dimnames(mod_dat$groups)[["config"]]
single_source_configs <- intersect(c("4", "8h", "8v", "16", "64"), available_configs)
multi_source_configs <- intersect(c("4", "8h", "8v", "16"), available_configs)

analysis_grid <- map_dfr(treatments, function(treat) {
  configs <- if (as.integer(treat) == 1L) {
    single_source_configs
  } else {
    multi_source_configs
  }
  expand.grid(
    block = blocks,
    treat = treat,
    config = configs,
    stringsAsFactors = FALSE
  ) |>
    rowwise() |>
    mutate(
      n_sources = length(get_true_source_groups(
        mod_dat, block, treat, config
      ))
    ) |>
    ungroup()
})

kappa_grid <- exp(seq(log(0.5), log(2), length.out = 5))

# Fit each full forward model once and reuse its transition-specific estimates
# across every spatial resolution for that experiment.
forward_grid <- distinct(analysis_grid, block, treat)
forward_checkpoint <- here(
  "DataProcessed", "results", "source_detection_main", "forward_fits.rds"
)
expected_forward_names <- paste(
  forward_grid$block, forward_grid$treat, sep = "_"
)
if (file.exists(forward_checkpoint)) {
  saved_forward_results <- readRDS(forward_checkpoint)
} else {
  saved_forward_results <- NULL
}
if (
  !is.null(saved_forward_results) &&
  identical(names(saved_forward_results), expected_forward_names)
) {
  forward_results <- saved_forward_results
  message("Using ", length(forward_results), " checkpointed forward fits")
} else {
  forward_results <- pmap(forward_grid, function(block, treat) {
    message("Forward fit: block ", block, ", treatment ", treat)
    fit_forward_experiment(block, treat, mod_dat, kappa_grid = kappa_grid)
  })
  names(forward_results) <- expected_forward_names
  dir.create(dirname(forward_checkpoint), recursive = TRUE, showWarnings = FALSE)
  saveRDS(forward_results, forward_checkpoint)
}

# Older checkpoints may contain the optimized forward fits but not uncertainty
# diagnostics. Add numerical Hessians without repeating the optimizations.
has_hessian_diagnostics <- vapply(forward_results, function(result) {
  "hessian_usable_vcov" %in% names(result$summary)
}, logical(1))
if (!all(has_hessian_diagnostics)) {
  for (i in which(!has_hessian_diagnostics)) {
    forward_results[[i]] <- add_forward_hessian_diagnostics(
      forward_fit = forward_results[[i]],
      block = forward_grid$block[i],
      treat = forward_grid$treat[i],
      mod_dat = mod_dat
    )
  }
  saveRDS(forward_results, forward_checkpoint)
  message("Added numerical Hessian diagnostics to the forward-fit checkpoint")
}

# Fit both complementary temporal views. Cumulative fits use all information
# available through an endpoint; transition fits use only the indicated
# interval. The final cumulative endpoint is the whole-study analysis.
transition_visits <- dimnames(mod_dat$intensity)[["visit"]][-1]
fit_grid <- bind_rows(
  merge(
    analysis_grid,
    data.frame(endpoint_visit = transition_visits),
    by = NULL
  ) |>
    mutate(analysis_mode = "cumulative"),
  merge(
    analysis_grid,
    data.frame(endpoint_visit = transition_visits),
    by = NULL
  ) |>
    mutate(analysis_mode = "transition")
) |>
  arrange(block, treat, config, analysis_mode, endpoint_visit)

requested_cores <- suppressWarnings(as.integer(Sys.getenv(
  "SOURCE_DETECTION_CORES",
  unset = NA_character_
)))
detected_cores <- parallel::detectCores(logical = FALSE)
if (is.na(detected_cores)) detected_cores <- 1L
n_cores <- if (is.na(requested_cores)) {
  max(1L, detected_cores - 1L)
} else {
  max(1L, requested_cores)
}
n_cores <- min(n_cores, nrow(fit_grid))

output_directory <- here("DataProcessed", "results", "source_detection_main")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

# Each row is an independent EM fit and reads the shared forward results and
# model data without modifying them. Inner transition optimizations remain
# sequential, avoiding nested parallelism and CPU oversubscription.
checkpoint_file <- file.path(output_directory, "em_fits.rds")
if (file.exists(checkpoint_file)) {
  saved_results <- readRDS(checkpoint_file)
  if (length(saved_results) == nrow(fit_grid)) {
    em_results <- saved_results
  } else {
    em_results <- vector("list", nrow(fit_grid))
  }
} else {
  em_results <- vector("list", nrow(fit_grid))
}

fits_to_run <- which(vapply(em_results, function(result) {
  is.null(result) || inherits(result, "try-error")
}, logical(1)))

message(
  "Running ", length(fits_to_run), " of ", nrow(fit_grid),
  " EM fits on ", n_cores, " core(s)"
)
new_results <- parallel::mclapply(fits_to_run, function(i) {
  fit_spec <- fit_grid[i, ]
  block <- fit_spec$block[[1]]
  treat <- fit_spec$treat[[1]]
  config <- fit_spec$config[[1]]
  n_sources <- fit_spec$n_sources[[1]]
  endpoint_visit <- fit_spec$endpoint_visit[[1]]
  analysis_mode <- fit_spec$analysis_mode[[1]]

  message(
    "EM fit: block ", block, ", treatment ", treat,
    ", config ", config, ", sources ", n_sources,
    ", mode ", analysis_mode, ", endpoint ", endpoint_visit
  )
  forward_fit <- forward_results[[paste(block, treat, sep = "_")]]
  fit_source_detection_experiment(
    block = block,
    treat = treat,
    config = config,
    n_sources = n_sources,
    mod_dat = mod_dat,
    theta_initial = forward_fit$theta,
    through_visit = if (analysis_mode == "cumulative") {
      endpoint_visit
    } else {
      NULL
    },
    transition_visit = if (analysis_mode == "transition") {
      endpoint_visit
    } else {
      NULL
    }
  )
}, mc.cores = n_cores, mc.preschedule = FALSE)
em_results[fits_to_run] <- new_results

# Preserve the expensive fits before constructing any combined summaries.
saveRDS(em_results, checkpoint_file)

fit_status <- fit_grid |>
  mutate(
    success = !vapply(em_results, inherits, logical(1), "try-error"),
    error = vapply(em_results, function(result) {
      if (inherits(result, "try-error")) as.character(result) else NA_character_
    }, character(1))
  )
write_csv(fit_status, file.path(output_directory, "fit_status.csv"))
if (any(!fit_status$success)) {
  stop(
    sum(!fit_status$success),
    " EM fits failed. See fit_status.csv; successful fits remain checkpointed."
  )
}

forward_summary <- bind_rows(map(forward_results, "summary"))
candidate_summary <- bind_rows(map(em_results, "candidate_summary"))
parameter_summary <- bind_rows(map(em_results, "parameter_summary")) |>
  select(-any_of(c("n_parameters", "n_observations", "bic")))
responsibility_history <- bind_rows(map(em_results, "responsibility_history"))
diagnostic_posteriors <- bind_rows(map(em_results, "diagnostic_posteriors"))
em_history <- bind_rows(map(em_results, "history"))
accuracy_summary <- bind_rows(map(em_results, "accuracy_summary")) |>
  select(-any_of(c("n_parameters", "n_observations", "bic")))
all_primary_assignments <- bind_rows(map(em_results, "whole_study_assignment"))

# Keep retrospective transition diagnostics only from the final study fit.
# Earlier endpoints are represented by genuine cumulative refits above.
whole_study_posteriors <- candidate_summary |>
  filter(fit_scope == "whole_study")
cumulative_fit_posteriors <- candidate_summary |>
  filter(fit_scope %in% c("cumulative_refit", "whole_study"))
transition_fit_posteriors <- candidate_summary |>
  filter(fit_scope == "transition_refit")
primary_fit_accuracy <- accuracy_summary |>
  filter(evidence %in% c(
    "transition_refit", "cumulative_refit", "whole_study"
  ))
whole_study_assignments_only <- all_primary_assignments |>
  filter(fit_scope == "whole_study")
full_study_diagnostics <- diagnostic_posteriors |>
  filter(fit_scope == "whole_study")

write_csv(forward_summary, file.path(output_directory, "forward_parameters.csv"))
write_csv(
  cumulative_fit_posteriors,
  file.path(output_directory, "cumulative_fit_posteriors.csv")
)
write_csv(
  transition_fit_posteriors,
  file.path(output_directory, "transition_fit_posteriors.csv")
)
write_csv(
  candidate_summary,
  file.path(output_directory, "all_primary_fit_posteriors.csv")
)
write_csv(
  whole_study_posteriors,
  file.path(output_directory, "whole_study_posteriors.csv")
)
write_csv(
  parameter_summary,
  file.path(output_directory, "em_parameters.csv")
)
write_csv(
  full_study_diagnostics,
  file.path(output_directory, "transition_and_cumulative_posteriors.csv")
)
write_csv(
  responsibility_history,
  file.path(output_directory, "em_responsibility_history.csv")
)
write_csv(
  em_history,
  file.path(output_directory, "em_convergence_history.csv")
)
write_csv(
  primary_fit_accuracy,
  file.path(output_directory, "source_accuracy.csv")
)
write_csv(
  accuracy_summary,
  file.path(output_directory, "all_fit_diagnostic_accuracy.csv")
)
write_csv(
  whole_study_assignments_only,
  file.path(output_directory, "whole_study_source_assignments.csv")
)
write_csv(
  all_primary_assignments,
  file.path(output_directory, "all_primary_fit_source_assignments.csv")
)
saveRDS(forward_results, file.path(output_directory, "forward_fits.rds"))

message("Main source-detection analysis complete: ", output_directory)
