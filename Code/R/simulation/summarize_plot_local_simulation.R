# Tables and diagnostic plots for the local whole-epidemic simulation study.

library(dplyr)
library(ggplot2)
library(here)
library(readr)
library(tidyr)

result_directory <- here(
  "DataProcessed", "results", "simulation", "whole_epidemic_study", "local"
)
table_directory <- file.path(result_directory, "tables")
plot_directory <- file.path(result_directory, "plots")
dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_directory, recursive = TRUE, showWarnings = FALSE)

forward <- read_csv(
  file.path(result_directory, "forward_metrics.csv"), show_col_types = FALSE
)
accuracy <- read_csv(
  file.path(result_directory, "backward_accuracy.csv"), show_col_types = FALSE
)
backward <- read_csv(
  file.path(result_directory, "backward_diagnostics.csv"), show_col_types = FALSE
)
fit_status <- read_csv(
  file.path(result_directory, "backward_fit_status.csv"), show_col_types = FALSE
)
progression <- read_csv(
  file.path(result_directory, "progression_diagnostics.csv"), show_col_types = FALSE
)
epidemic <- read_csv(
  file.path(result_directory, "epidemic_summary.csv"), show_col_types = FALSE
)
mod_dat <- readRDS(here("DataProcessed", "experimental", "mod_dat_arrays.rds"))

fit_scope_labels <- c(
  transition_refit = "Independent transition",
  cumulative_refit = "Cumulative through visit",
  whole_study = "Whole study"
)
config_levels <- c("4", "8h", "8v", "16", "64")
config_labels <- c("4 cells", "8 horizontal", "8 vertical", "16 cells", "64 plants")
base_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "bottom"
  )

mcse_mean <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))
safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

# Forward-model parameter recovery --------------------------------------
parameter_specification <- data.frame(
  parameter = c("beta", "delta", "gamma", "log(kappa)", "log(phi)"),
  estimate = c("beta", "delta", "gamma", "kappa", "phi"),
  truth = c("true_beta", "true_delta", "true_gamma", "true_kappa", "true_phi"),
  bias = c("beta_bias", "delta_bias", "gamma_bias", "log_kappa_bias", "log_phi_bias"),
  se = c("se_beta", "se_delta", "se_gamma", "se_log_kappa", "se_log_phi"),
  covered = c("beta_wald_covered", "delta_wald_covered", "gamma_wald_covered",
              "kappa_wald_covered", "phi_wald_covered"),
  stringsAsFactors = FALSE
)

forward_long <- bind_rows(lapply(seq_len(nrow(parameter_specification)), function(i) {
  spec <- parameter_specification[i, ]
  estimate <- forward[[spec$estimate]]
  truth <- forward[[spec$truth]]
  if (spec$parameter == "log(kappa)") {
    estimate <- log(estimate)
    truth <- log(truth)
  }
  if (spec$parameter == "log(phi)") {
    estimate <- log(estimate)
    truth <- log(truth)
  }
  data.frame(
    simulation = forward$simulation,
    block = forward$block,
    treat = as.character(forward$treat),
    visit = forward$visit,
    parameter = spec$parameter,
    estimate = estimate,
    truth = truth,
    bias = forward[[spec$bias]],
    reported_se = forward[[spec$se]],
    covered = forward[[spec$covered]],
    hessian_usable = forward$hessian_usable_vcov
  )
}))

forward_parameter_summary <- forward_long |>
  group_by(block, treat, visit, parameter) |>
  summarise(
    n_simulations = n(),
    n_usable_hessians = sum(hessian_usable),
    hessian_usable_proportion = mean(hessian_usable),
    mean_estimate = mean(estimate),
    empirical_sd = sd(estimate),
    mean_bias = mean(bias),
    median_bias = median(bias),
    rmse = sqrt(mean(bias^2)),
    bias_mcse = mcse_mean(bias),
    mean_reported_se = safe_mean(reported_se),
    se_to_empirical_sd = mean_reported_se / empirical_sd,
    wald_coverage = safe_mean(covered),
    .groups = "drop"
  )
write_csv(
  forward_parameter_summary,
  file.path(table_directory, "forward_parameter_recovery.csv")
)

forward_hessian_summary <- forward |>
  group_by(block, treat, visit) |>
  summarise(
    n_simulations = n(),
    optimizer_convergence_rate = mean(convergence == 0),
    positive_definite_rate = mean(hessian_positive_definite),
    inversion_success_rate = mean(hessian_inversion_succeeded),
    usable_vcov_rate = mean(hessian_usable_vcov),
    median_hessian_rank = median(hessian_rank),
    median_condition_number = median(hessian_condition_number, na.rm = TRUE),
    maximum_condition_number = max(hessian_condition_number, na.rm = TRUE),
    median_abs_gamma_log_kappa_correlation = median(
      abs(gamma_log_kappa_correlation), na.rm = TRUE
    ),
    mean_function_evaluations = mean(function_evaluations),
    mean_gradient_evaluations = mean(gradient_evaluations),
    .groups = "drop"
  )
write_csv(
  forward_hessian_summary,
  file.path(table_directory, "forward_hessian_diagnostics.csv")
)

# Backward-model prediction and numerical diagnostics -------------------
primary_accuracy <- accuracy |>
  filter(evidence == fit_scope) |>
  mutate(
    fit_type = factor(fit_scope_labels[fit_scope], levels = fit_scope_labels),
    config = factor(config, levels = config_levels, labels = config_labels),
    complete_source_set = standard_accuracy == 1,
    naive_complete_source_set = naive_standard_accuracy == 1
  )

summarize_accuracy <- function(data, include_block) {
  grouping <- c(
    if (include_block) "block",
    "treat", "config", "n_sources", "visit", "fit_scope"
  )
  data |>
    group_by(across(all_of(grouping))) |>
    summarise(
      n_simulations = n(),
      exact_set_recovery = mean(complete_source_set),
      exact_set_recovery_mcse = mcse_mean(as.numeric(complete_source_set)),
      mean_component_accuracy = mean(standard_accuracy),
      component_accuracy_mcse = mcse_mean(standard_accuracy),
      mean_distance_weighted_accuracy = mean(distance_weighted_accuracy),
      distance_weighted_accuracy_mcse = mcse_mean(distance_weighted_accuracy),
      naive_exact_set_recovery = mean(naive_complete_source_set),
      naive_mean_component_accuracy = mean(naive_standard_accuracy),
      naive_mean_distance_weighted_accuracy = mean(
        naive_distance_weighted_accuracy
      ),
      mean_model_minus_naive_distance_accuracy = mean(
        distance_weighted_accuracy - naive_distance_weighted_accuracy
      ),
      .groups = "drop"
    )
}
write_csv(
  summarize_accuracy(primary_accuracy, TRUE),
  file.path(table_directory, "backward_prediction_by_block.csv")
)
write_csv(
  summarize_accuracy(primary_accuracy, FALSE),
  file.path(table_directory, "backward_prediction_overall.csv")
)

backward_diagnostic_summary <- backward |>
  group_by(block, treat, config, n_sources, through_visit, fit_scope) |>
  summarise(
    n_simulations = n(),
    convergence_rate = mean(converged),
    monotonic_rate = mean(all_iterations_monotonic),
    m_step_failure_rate = mean(any_m_step_failure),
    median_iterations = median(iterations),
    maximum_iterations = max(iterations),
    median_true_rank = median(true_rank),
    mean_true_posterior = mean(true_posterior),
    mean_winning_posterior = mean(winning_posterior),
    mean_effective_candidates = mean(effective_candidates),
    .groups = "drop"
  )
write_csv(
  backward_diagnostic_summary,
  file.path(table_directory, "backward_em_diagnostics.csv")
)

group_count <- c(`4` = 4, `8h` = 8, `8v` = 8, `16` = 16, `64` = 64)
runtime_summary <- fit_status |>
  mutate(
    n_candidates = choose(group_count[config], n_sources),
    fit_scope = if_else(
      analysis_mode == "transition",
      "transition_refit",
      if_else(endpoint_visit == max(endpoint_visit), "whole_study", "cumulative_refit")
    )
  ) |>
  group_by(config, n_sources, endpoint_visit, fit_scope, n_candidates) |>
  summarise(
    n_fits = n(),
    success_rate = mean(success),
    median_seconds = median(elapsed_seconds),
    maximum_seconds = max(elapsed_seconds),
    .groups = "drop"
  )
write_csv(runtime_summary, file.path(table_directory, "runtime_summary.csv"))

progression_summary <- progression |>
  group_by(block, treat, visit) |>
  summarise(
    n_simulations = n(),
    mean_proportion_decreased = mean(proportion_decreased),
    maximum_proportion_decreased = max(proportion_decreased),
    mean_decrease_when_present = safe_mean(if_else(
      n_decreased > 0, mean_decrease, NA_real_
    )),
    maximum_decrease = max(maximum_decrease),
    total_positive_to_zero = sum(n_positive_to_zero),
    positive_to_zero_rate = mean(proportion_positive_to_zero),
    .groups = "drop"
  )
write_csv(
  progression_summary,
  file.path(table_directory, "nonmonotone_progression_diagnostics.csv")
)

# Figures ---------------------------------------------------------------
bias_plot <- ggplot(
  forward_parameter_summary,
  aes(as.factor(visit), mean_bias, color = treat, group = treat)
) +
  geom_hline(yintercept = 0, color = "grey50") +
  geom_errorbar(
    aes(ymin = mean_bias - 1.96 * bias_mcse,
        ymax = mean_bias + 1.96 * bias_mcse),
    width = 0.15
  ) +
  geom_point(size = 2) +
  geom_line() +
  facet_grid(parameter ~ block, scales = "free_y") +
  labs(
    title = "Forward-model parameter bias",
    subtitle = "Points are means; intervals show Monte Carlo uncertainty (five-replicate pilot)",
    x = "Transition ending at visit", y = "Estimate minus generating value",
    color = "Sources"
  ) + base_theme
ggsave(file.path(plot_directory, "forward_parameter_bias.png"), bias_plot,
       width = 12, height = 11, dpi = 300)

estimate_plot <- ggplot(
  forward_long,
  aes(truth, estimate, color = as.factor(visit))
) +
  geom_abline(slope = 1, intercept = 0, color = "grey50") +
  geom_point(alpha = 0.65) +
  facet_grid(parameter ~ block, scales = "free") +
  labs(
    title = "Estimated versus generating forward parameters",
    x = "Generating value", y = "Estimate", color = "Visit"
  ) + base_theme
ggsave(file.path(plot_directory, "forward_estimate_vs_truth.png"), estimate_plot,
       width = 12, height = 11, dpi = 300)

hessian_plot <- ggplot(
  forward_hessian_summary,
  aes(as.factor(visit), factor(treat), fill = usable_vcov_rate)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", usable_vcov_rate))) +
  facet_wrap(~block) +
  scale_fill_viridis_c(limits = c(0, 1)) +
  labs(
    title = "Proportion of forward fits with a usable Hessian",
    x = "Transition ending at visit", y = "Planted sources", fill = "Proportion"
  ) + base_theme
ggsave(file.path(plot_directory, "forward_hessian_usable_rate.png"), hessian_plot,
       width = 10, height = 7, dpi = 300)

coverage_data <- forward_parameter_summary |>
  mutate(parameter = factor(
    parameter, levels = c("beta", "delta", "gamma", "log(kappa)", "log(phi)")
  ))
coverage_plot <- ggplot(
  coverage_data,
  aes(as.factor(visit), parameter, fill = wald_coverage)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = if_else(
    is.na(wald_coverage), "NA", sprintf("%.2f", wald_coverage)
  )), size = 3) +
  facet_grid(block ~ treat) +
  scale_fill_viridis_c(limits = c(0, 1), na.value = "grey85") +
  labs(
    title = "Forward-model 95% Wald interval coverage",
    subtitle = "Coverage is conditional on a usable covariance matrix",
    x = "Transition ending at visit", y = NULL, fill = "Coverage"
  ) + base_theme
ggsave(file.path(plot_directory, "forward_wald_coverage.png"), coverage_plot,
       width = 11, height = 10, dpi = 300)

accuracy_heatmap_data <- primary_accuracy |>
  group_by(treat, config, visit, fit_type) |>
  summarise(exact_set_recovery = mean(complete_source_set), .groups = "drop")
accuracy_heatmap <- ggplot(
  accuracy_heatmap_data,
  aes(as.factor(visit), config, fill = exact_set_recovery)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", exact_set_recovery)), size = 3) +
  facet_grid(treat ~ fit_type, scales = "free_x", space = "free_x") +
  scale_fill_viridis_c(limits = c(0, 1)) +
  labs(
    title = "Exact recovery of the complete source set",
    subtitle = "Averaged over blocks and five simulations per block",
    x = "Endpoint visit", y = "Spatial resolution", fill = "Recovery"
  ) + base_theme
ggsave(file.path(plot_directory, "backward_exact_recovery_heatmap.png"),
       accuracy_heatmap, width = 13, height = 8, dpi = 300)

whole_accuracy <- primary_accuracy |>
  filter(fit_scope == "whole_study")
distance_plot <- ggplot(
  whole_accuracy,
  aes(config, distance_weighted_accuracy, fill = block)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.45) +
  geom_point(
    aes(color = block),
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
    size = 1.4, alpha = 0.8
  ) +
  facet_wrap(~treat, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Whole-study distance-weighted source accuracy",
    subtitle = "Each point is one simulated epidemic",
    x = "Spatial resolution", y = "Distance-weighted accuracy",
    fill = "Block", color = "Block"
  ) + base_theme
ggsave(file.path(plot_directory, "whole_study_distance_accuracy.png"),
       distance_plot, width = 12, height = 6, dpi = 300)

posterior_plot <- backward |>
  mutate(
    fit_type = factor(fit_scope_labels[fit_scope], levels = fit_scope_labels),
    config = factor(config, levels = config_levels, labels = config_labels)
  ) |>
  ggplot(aes(as.factor(through_visit), true_posterior, color = block)) +
  geom_boxplot(aes(group = interaction(block, through_visit)), outlier.shape = NA) +
  geom_point(
    position = position_jitter(width = 0.1, height = 0),
    alpha = 0.55, size = 1
  ) +
  facet_grid(treat + fit_type ~ config, scales = "free_x", space = "free_x") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Posterior probability of the true source set",
    x = "Endpoint visit", y = "True-source posterior", color = "Block"
  ) + base_theme
ggsave(file.path(plot_directory, "backward_true_source_posterior.png"),
       posterior_plot, width = 14, height = 12, dpi = 300)

iteration_plot <- backward |>
  mutate(
    fit_type = factor(fit_scope_labels[fit_scope], levels = fit_scope_labels),
    config = factor(config, levels = config_levels, labels = config_labels)
  ) |>
  ggplot(aes(as.factor(through_visit), iterations, color = block)) +
  geom_boxplot(aes(group = interaction(block, through_visit)), outlier.shape = NA) +
  geom_point(
    position = position_jitter(width = 0.1, height = 0),
    alpha = 0.5, size = 1
  ) +
  facet_grid(treat + fit_type ~ config, scales = "free") +
  labs(
    title = "Backward EM iteration counts",
    x = "Endpoint visit", y = "Iterations", color = "Block"
  ) + base_theme
ggsave(file.path(plot_directory, "backward_em_iterations.png"), iteration_plot,
       width = 14, height = 12, dpi = 300)

runtime_plot <- fit_status |>
  mutate(
    n_candidates = choose(group_count[config], n_sources),
    fit_type = if_else(analysis_mode == "transition", "Transition", "Cumulative")
  ) |>
  ggplot(aes(n_candidates, elapsed_seconds, color = fit_type)) +
  geom_point(alpha = 0.35) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~treat) +
  labs(
    title = "Source-detection runtime and candidate-set size",
    x = "Number of candidate source sets (log scale)",
    y = "Elapsed seconds (log scale)", color = "Fit"
  ) + base_theme
ggsave(file.path(plot_directory, "backward_runtime.png"), runtime_plot,
       width = 10, height = 6, dpi = 300)

progression_plot <- ggplot(
  progression,
  aes(as.factor(visit), proportion_decreased, color = block)
) +
  geom_boxplot(aes(group = interaction(block, visit)), outlier.shape = NA) +
  geom_point(
    position = position_jitter(width = 0.1, height = 0), alpha = 0.7
  ) +
  facet_wrap(~treat) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Non-monotone disease progression in simulated epidemics",
    x = "Transition ending at visit", y = "Proportion of plants decreasing",
    color = "Block"
  ) + base_theme
ggsave(file.path(plot_directory, "nonmonotone_progression.png"), progression_plot,
       width = 10, height = 6, dpi = 300)

# Compare simulated mean severity with the observed generating scenarios.
observed_epidemic <- expand.grid(
  block = dimnames(mod_dat$intensity)[["block"]],
  treat = dimnames(mod_dat$intensity)[["treat"]],
  visit = dimnames(mod_dat$intensity)[["visit"]],
  stringsAsFactors = FALSE
) |>
  rowwise() |>
  mutate(mean_severity = mean(mod_dat$intensity[, block, treat, visit])) |>
  ungroup() |>
  mutate(
    block = as.character(block),
    treat = as.character(treat),
    visit = as.character(visit)
  )
simulated_epidemic <- epidemic |>
  mutate(
    block = as.character(block),
    treat = as.character(treat),
    visit = as.character(visit)
  ) |>
  group_by(block, treat, visit) |>
  summarise(
    simulated_mean = mean(mean_severity),
    simulated_minimum = min(mean_severity),
    simulated_maximum = max(mean_severity),
    .groups = "drop"
  ) |>
  left_join(observed_epidemic, by = c("block", "treat", "visit"))
epidemic_plot <- ggplot(
  simulated_epidemic,
  aes(as.numeric(visit), simulated_mean, color = block, group = block)
) +
  geom_ribbon(
    aes(ymin = simulated_minimum, ymax = simulated_maximum, fill = block),
    alpha = 0.12, color = NA
  ) +
  geom_line() + geom_point() +
  geom_point(aes(y = mean_severity), shape = 4, size = 2.5) +
  facet_wrap(~treat) +
  labs(
    title = "Observed and simulated epidemic trajectories",
    subtitle = "Lines are simulation means, ribbons are ranges, and crosses are observed values",
    x = "Visit", y = "Mean disease severity", color = "Block", fill = "Block"
  ) + base_theme
ggsave(file.path(plot_directory, "epidemic_trajectory_check.png"), epidemic_plot,
       width = 11, height = 6, dpi = 300)

message("Simulation tables saved to: ", table_directory)
message("Simulation plots saved to: ", plot_directory)
