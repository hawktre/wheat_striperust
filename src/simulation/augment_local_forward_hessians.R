# Add forward-model Hessian diagnostics to a completed local simulation run.
# This recreates simulated responses from deterministic seeds but does not
# repeat either the forward optimization or any backward EM fit.

library(dplyr)
library(here)
library(readr)

source(here("src", "simulation", "simulation_study_functions.R"))

result_directory <- here(
  "output", "simulation", "whole_epidemic_study", "local"
)
checkpoint_file <- file.path(result_directory, "local_checkpoint.rds")
forward_file <- file.path(result_directory, "forward_metrics.csv")
archive_directory <- here(
  "output", "simulation", "whole_epidemic_study", "archive"
)
dir.create(archive_directory, recursive = TRUE, showWarnings = FALSE)

mod_dat <- readRDS(here("data", "processed", "experimental", "mod_dat_arrays.rds"))
generating_fits <- readRDS(here(
  "output", "source_detection", "forward_fits.rds"
))
results <- readRDS(checkpoint_file)
forward_metrics <- read_csv(forward_file, show_col_types = FALSE)

scenario_keys <- forward_metrics |>
  distinct(simulation, block, treat) |>
  arrange(simulation, block, treat)

updated_metrics <- vector("list", nrow(scenario_keys))
for (i in seq_len(nrow(scenario_keys))) {
  key <- scenario_keys[i, ]
  simulation <- key$simulation
  block <- key$block
  treat <- as.character(key$treat)
  message(
    "Hessian diagnostics: simulation ", simulation,
    ", block ", block, ", treatment ", treat
  )

  generating_fit <- generating_fits[[paste(block, treat, sep = "_")]]
  scenario_seed <- 100000L * as.integer(simulation) +
    1000L * match(block, dimnames(mod_dat$intensity)[["block"]]) +
    10L * match(treat, dimnames(mod_dat$intensity)[["treat"]])
  simulated <- simulate_whole_epidemic(
    mod_dat, generating_fit, block, treat,
    mechanism = "secondary_spread", seed = scenario_seed
  )
  simulated_mod_dat <- replace_experiment_intensity(
    mod_dat, block, treat, simulated
  )
  experiment <- extract_experiment(simulated_mod_dat, block, treat)

  scenario_metrics <- forward_metrics |>
    filter(
      .data$simulation == .env$simulation,
      .data$block == .env$block,
      as.character(.data$treat) == .env$treat
    ) |>
    arrange(as.numeric(visit))
  theta <- cbind(
    beta = scenario_metrics$beta,
    delta = scenario_metrics$delta,
    gamma = scenario_metrics$gamma,
    log_kappa = log(scenario_metrics$kappa),
    log_phi = log(scenario_metrics$phi)
  )
  rownames(theta) <- as.character(scenario_metrics$visit)

  diagnostics <- lapply(seq_along(experiment$transition_visits), function(t) {
    forward_hessian_diagnostics(
      theta = theta[t, ],
      y_current = experiment$y[[t + 1L]],
      y_previous = experiment$y[[t]],
      wind_matrix = experiment$wind[[t]],
      distance_matrix = experiment$distance
    )
  })
  standard_errors <- do.call(rbind, lapply(
    diagnostics, `[[`, "standard_errors"
  ))

  scenario_metrics$hessian_rank <- vapply(
    diagnostics, `[[`, integer(1), "rank"
  )
  scenario_metrics$hessian_minimum_eigenvalue <- vapply(
    diagnostics, `[[`, numeric(1), "minimum_eigenvalue"
  )
  scenario_metrics$hessian_maximum_eigenvalue <- vapply(
    diagnostics, `[[`, numeric(1), "maximum_eigenvalue"
  )
  scenario_metrics$hessian_condition_number <- vapply(
    diagnostics, `[[`, numeric(1), "condition_number"
  )
  scenario_metrics$hessian_positive_definite <- vapply(
    diagnostics, `[[`, logical(1), "positive_definite"
  )
  scenario_metrics$hessian_inversion_succeeded <- vapply(
    diagnostics, `[[`, logical(1), "inversion_succeeded"
  )
  scenario_metrics$hessian_usable_vcov <- vapply(
    diagnostics, `[[`, logical(1), "usable_vcov"
  )
  scenario_metrics$se_beta <- standard_errors[, "beta"]
  scenario_metrics$se_delta <- standard_errors[, "delta"]
  scenario_metrics$se_gamma <- standard_errors[, "gamma"]
  scenario_metrics$se_log_kappa <- standard_errors[, "log_kappa"]
  scenario_metrics$se_log_phi <- standard_errors[, "log_phi"]
  scenario_metrics$se_kappa <- scenario_metrics$kappa *
    scenario_metrics$se_log_kappa
  scenario_metrics$se_phi <- scenario_metrics$phi *
    scenario_metrics$se_log_phi
  scenario_metrics$gamma_log_kappa_correlation <- vapply(
    diagnostics, `[[`, numeric(1), "gamma_log_kappa_correlation"
  )
  scenario_metrics <- scenario_metrics |>
    mutate(
      beta_wald_covered = abs(beta_bias) <= 1.96 * se_beta,
      delta_wald_covered = abs(delta_bias) <= 1.96 * se_delta,
      gamma_wald_covered = abs(gamma_bias) <= 1.96 * se_gamma,
      kappa_wald_covered = abs(log_kappa_bias) <= 1.96 * se_log_kappa,
      phi_wald_covered = abs(log_phi_bias) <= 1.96 * se_log_phi,
      beta_standardized_bias = beta_bias / se_beta,
      delta_standardized_bias = delta_bias / se_delta,
      gamma_standardized_bias = gamma_bias / se_gamma,
      log_kappa_standardized_bias = log_kappa_bias / se_log_kappa,
      log_phi_standardized_bias = log_phi_bias / se_log_phi
    )
  updated_metrics[[i]] <- scenario_metrics
}

forward_metrics <- bind_rows(updated_metrics) |>
  arrange(simulation, block, treat, as.numeric(visit))

# Preserve the completed pre-Hessian checkpoint before updating its compact
# forward summaries. The expensive backward results are unchanged.
file.copy(
  checkpoint_file,
  file.path(archive_directory, "local_checkpoint_before_hessian_augmentation.rds"),
  overwrite = TRUE
)
for (i in seq_along(results)) {
  if (is.null(results[[i]]) || is.null(results[[i]]$forward_metrics)) next
  key <- results[[i]]$scenario_status[1, ]
  results[[i]]$forward_metrics <- forward_metrics |>
    filter(
      simulation == key$simulation,
      block == key$block,
      as.character(treat) == as.character(key$treat)
    )
}

write_csv(forward_metrics, forward_file)
saveRDS(results, checkpoint_file)
message("Added Hessian diagnostics to: ", forward_file)
