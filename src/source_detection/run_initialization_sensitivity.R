# Whole-study sensitivity analysis for direct EM initialization.
#
# Unlike the primary analysis, this script does not fit the full forward
# model. For each starting kappa, it constructs transition-specific
# method-of-moments estimates and passes them directly to the same EM fitter
# used by the primary analysis. All other model components are unchanged.

library(dplyr)
library(here)
library(parallel)
library(purrr)
library(readr)

source(here("src", "source_detection", "source_detection_functions.R"))

mod_dat <- readRDS(here("data", "processed", "experimental", "mod_dat_arrays.rds"))

blocks <- dimnames(mod_dat$intensity)[["block"]]
treatments <- dimnames(mod_dat$intensity)[["treat"]]
available_configs <- dimnames(mod_dat$groups)[["config"]]
single_source_configs <- intersect(
  c("4", "8h", "8v", "16", "64"), available_configs
)
multi_source_configs <- intersect(
  c("4", "8h", "8v", "16"), available_configs
)

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

# These match the range used to initialize the primary forward fits. Here,
# however, no forward optimization is performed and every start is retained.
initial_kappa_grid <- exp(seq(log(0.5), log(2), length.out = 5))
fit_grid <- merge(
  analysis_grid,
  data.frame(initial_kappa = initial_kappa_grid),
  by = NULL
) |>
  arrange(block, treat, config, initial_kappa)

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

output_directory <- here(
  "output", "source_detection_initialization_sensitivity"
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
checkpoint_file <- file.path(output_directory, "initialization_fits.rds")

message(
  "Running ", nrow(fit_grid),
  " direct-initialization whole-study fits on ", n_cores, " core(s)"
)

fit_results <- parallel::mclapply(seq_len(nrow(fit_grid)), function(i) {
  fit_spec <- fit_grid[i, ]
  block <- fit_spec$block[[1]]
  treat <- fit_spec$treat[[1]]
  config <- fit_spec$config[[1]]
  n_sources <- fit_spec$n_sources[[1]]
  initial_kappa <- fit_spec$initial_kappa[[1]]

  message(
    "Direct EM start: block ", block, ", treatment ", treat,
    ", config ", config, ", sources ", n_sources,
    ", initial kappa ", signif(initial_kappa, 4)
  )

  tryCatch({
    theta_initial <- initialize_experiment_theta(
      block = block,
      treat = treat,
      mod_dat = mod_dat,
      initial_kappa = initial_kappa
    )
    fit <- fit_source_detection_experiment(
      block = block,
      treat = treat,
      config = config,
      n_sources = n_sources,
      mod_dat = mod_dat,
      theta_initial = theta_initial,
      include_naive = FALSE
    )
    list(
      success = TRUE,
      initial_kappa = initial_kappa,
      theta_initial = theta_initial,
      fit = fit,
      error = NA_character_
    )
  }, error = function(e) {
    list(
      success = FALSE,
      initial_kappa = initial_kappa,
      theta_initial = NULL,
      fit = NULL,
      error = conditionMessage(e)
    )
  })
}, mc.cores = n_cores, mc.preschedule = FALSE)

# Checkpoint the expensive model fits before constructing any summaries.
# If downstream reporting fails, this object can be loaded without refitting.
saveRDS(fit_results, checkpoint_file)

fit_status <- bind_cols(
  fit_grid,
  tibble(
    success = vapply(fit_results, `[[`, logical(1), "success"),
    error = vapply(fit_results, `[[`, character(1), "error")
  )
)

successful_results <- fit_results[vapply(
  fit_results, `[[`, logical(1), "success"
)]
if (!length(successful_results)) {
  stop("Every direct-initialization sensitivity fit failed.")
}

add_start <- function(result, component) {
  result$fit[[component]] |>
    mutate(initial_kappa = result$initial_kappa, .before = 1)
}

candidate_summary <- bind_rows(map(
  successful_results, add_start, "candidate_summary"
))
parameter_summary <- bind_rows(map(
  successful_results, add_start, "parameter_summary"
))
responsibility_history <- bind_rows(map(
  successful_results, add_start, "responsibility_history"
))
em_history <- bind_rows(map(successful_results, add_start, "history"))
accuracy_summary <- bind_rows(map(
  successful_results, add_start, "accuracy_summary"
)) |>
  filter(evidence == "whole_study")
assignment_summary <- bind_rows(map(
  successful_results, add_start, "whole_study_assignment"
))
initial_parameter_summary <- bind_rows(map(successful_results, function(x) {
  as.data.frame(x$theta_initial) |>
    mutate(
      initial_kappa = x$initial_kappa,
      visit = rownames(x$theta_initial),
      .before = 1
    )
}))

# Compare the direct starts with one another at the scenario level.
start_comparison <- candidate_summary |>
  filter(predicted) |>
  select(
    initial_kappa, block, treat, config, n_sources,
    predicted_candidate_id = candidate_id,
    winning_posterior = posterior,
    observed_candidate_loglik = loglik
  ) |>
  left_join(
    accuracy_summary |>
      select(
        initial_kappa, block, treat, config, n_sources,
        standard_accuracy, distance_weighted_accuracy
      ),
    by = c("initial_kappa", "block", "treat", "config", "n_sources")
  ) |>
  left_join(
    parameter_summary |>
      distinct(
        initial_kappa, block, treat, config, n_sources,
        observed_loglik, status, iterations
      ),
    by = c("initial_kappa", "block", "treat", "config", "n_sources")
  )

# If the primary analysis is available, record agreement with its whole-study
# prediction and the difference in the maximized observed-data likelihood.
main_directory <- here("output", "source_detection")
main_posterior_file <- file.path(main_directory, "whole_study_posteriors.csv")
main_parameter_file <- file.path(main_directory, "em_parameters.csv")
if (file.exists(main_posterior_file) && file.exists(main_parameter_file)) {
  main_winners <- read_csv(main_posterior_file, show_col_types = FALSE) |>
    filter(predicted) |>
    mutate(
      block = as.character(block),
      treat = as.character(treat),
      config = as.character(config),
      n_sources = as.integer(n_sources)
    ) |>
    select(
      block, treat, config, n_sources,
      main_candidate_id = candidate_id,
      main_winning_posterior = posterior
    )
  main_loglik <- read_csv(main_parameter_file, show_col_types = FALSE) |>
    filter(fit_scope == "whole_study") |>
    mutate(
      block = as.character(block),
      treat = as.character(treat),
      config = as.character(config),
      n_sources = as.integer(n_sources)
    ) |>
    distinct(block, treat, config, n_sources, main_observed_loglik = observed_loglik)

  start_comparison <- start_comparison |>
    mutate(
      block = as.character(block),
      treat = as.character(treat),
      config = as.character(config),
      n_sources = as.integer(n_sources)
    ) |>
    left_join(main_winners, by = c("block", "treat", "config", "n_sources")) |>
    left_join(main_loglik, by = c("block", "treat", "config", "n_sources")) |>
    mutate(
      agrees_with_main = predicted_candidate_id == main_candidate_id,
      observed_loglik_difference = observed_loglik - main_observed_loglik
    )
}

write_csv(fit_status, file.path(output_directory, "fit_status.csv"))
write_csv(
  initial_parameter_summary,
  file.path(output_directory, "initial_parameters.csv")
)
write_csv(
  candidate_summary,
  file.path(output_directory, "whole_study_posteriors.csv")
)
write_csv(parameter_summary, file.path(output_directory, "em_parameters.csv"))
write_csv(
  responsibility_history,
  file.path(output_directory, "em_responsibility_history.csv")
)
write_csv(em_history, file.path(output_directory, "em_convergence_history.csv"))
write_csv(accuracy_summary, file.path(output_directory, "source_accuracy.csv"))
write_csv(
  assignment_summary,
  file.path(output_directory, "source_assignments.csv")
)
write_csv(
  start_comparison,
  file.path(output_directory, "initialization_comparison.csv")
)
# The complete fit objects were checkpointed immediately after fitting.

message("Initialization sensitivity analysis complete: ", output_directory)
