# Diagnostics and results plots for the whole-study initialization sensitivity
# analysis. Run after run_initialization_sensitivity.R has completed.

library(dplyr)
library(ggplot2)
library(here)
library(readr)
library(tidyr)

result_directory <- here(
  "output", "source_detection_initialization_sensitivity"
)
main_directory <- here("output", "source_detection")
plot_directory <- file.path(result_directory, "plots")
dir.create(plot_directory, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  status = "fit_status.csv",
  comparison = "initialization_comparison.csv",
  parameters = "em_parameters.csv",
  history = "em_convergence_history.csv",
  posteriors = "whole_study_posteriors.csv"
)
required_paths <- setNames(
  file.path(result_directory, unname(required_files)),
  names(required_files)
)
missing_files <- required_files[!file.exists(required_paths)]
if (length(missing_files)) {
  stop(
    "Run run_initialization_sensitivity.R first. Missing result files: ",
    paste(missing_files, collapse = ", ")
  )
}

fit_status <- read_csv(required_paths[["status"]], show_col_types = FALSE)
comparison <- read_csv(required_paths[["comparison"]], show_col_types = FALSE)
parameters <- read_csv(required_paths[["parameters"]], show_col_types = FALSE)
em_history <- read_csv(required_paths[["history"]], show_col_types = FALSE)
posteriors <- read_csv(required_paths[["posteriors"]], show_col_types = FALSE)

config_levels <- c("4", "8h", "8v", "16", "64")
config_labels <- c(
  "4 cells", "8 horizontal cells", "8 vertical cells",
  "16 cells", "64 plants"
)
treat_levels <- c("1", "2", "4")
treat_labels <- c(
  "1 planted source", "2 planted sources", "4 planted sources"
)
kappa_breaks <- sort(unique(comparison$initial_kappa))

format_results <- function(data) {
  if ("block" %in% names(data)) {
    data <- data |> mutate(block = factor(block))
  }
  data |>
    mutate(
      treat = factor(
        as.character(treat), levels = treat_levels, labels = treat_labels
      ),
      config = factor(
        as.character(config), levels = config_levels, labels = config_labels
      )
    )
}

base_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "bottom"
  )

# Compact numerical summaries -------------------------------------------
fit_status_summary <- fit_status |>
  count(treat, config, success) |>
  group_by(treat, config) |>
  mutate(proportion = n / sum(n)) |>
  ungroup()

convergence_summary <- em_history |>
  group_by(initial_kappa, block, treat, config, n_sources) |>
  arrange(iteration, .by_group = TRUE) |>
  summarise(
    iterations = max(iteration),
    converged = any(converged),
    final_observed_loglik = observed_loglik[which.max(iteration)],
    final_loglik_change = loglik_change[which.max(iteration)],
    final_responsibility_change =
      max_responsibility_change[which.max(iteration)],
    all_iterations_monotonic = all(monotonic[!is.na(monotonic)]),
    any_m_step_failure = any(m_step_convergence > 0, na.rm = TRUE),
    .groups = "drop"
  )

posterior_summary <- posteriors |>
  group_by(initial_kappa, block, treat, config, n_sources) |>
  summarise(
    winner_posterior = max(posterior),
    true_posterior = posterior[true_candidate][1],
    true_rank = rank[true_candidate][1],
    posterior_entropy = -sum(if_else(
      posterior > 0, posterior * log(posterior), 0
    )),
    effective_candidates = exp(posterior_entropy),
    .groups = "drop"
  )

# Stability is evaluated within each block x treatment x resolution across
# the five direct starting values.
scenario_stability <- comparison |>
  group_by(block, treat, config, n_sources) |>
  summarise(
    n_starts = n_distinct(initial_kappa),
    n_distinct_predictions = n_distinct(predicted_candidate_id),
    all_starts_same_prediction = n_distinct_predictions == 1L,
    proportion_agreeing_with_main = mean(agrees_with_main, na.rm = TRUE),
    observed_loglik_range = diff(range(observed_loglik, na.rm = TRUE)),
    standard_accuracy_range = diff(range(standard_accuracy, na.rm = TRUE)),
    distance_accuracy_range =
      diff(range(distance_weighted_accuracy, na.rm = TRUE)),
    .groups = "drop"
  )

scenario_summary <- scenario_stability |>
  group_by(treat, config) |>
  summarise(
    n_scenarios = n(),
    proportion_stable = mean(all_starts_same_prediction),
    proportion_agreeing_with_main = mean(proportion_agreeing_with_main),
    mean_loglik_range = mean(observed_loglik_range),
    maximum_loglik_range = max(observed_loglik_range),
    mean_distance_accuracy_range = mean(distance_accuracy_range),
    .groups = "drop"
  )

write_csv(
  fit_status_summary,
  file.path(result_directory, "fit_status_summary.csv")
)
write_csv(
  convergence_summary,
  file.path(result_directory, "convergence_diagnostics.csv")
)
write_csv(
  posterior_summary,
  file.path(result_directory, "posterior_diagnostics.csv")
)
write_csv(
  scenario_stability,
  file.path(result_directory, "initialization_stability_by_scenario.csv")
)
write_csv(
  scenario_summary,
  file.path(result_directory, "initialization_stability_summary.csv")
)

# Agreement with the primary warm-start analysis ------------------------
agreement_plot_data <- scenario_summary |>
  format_results()

agreement_plot <- ggplot(
  agreement_plot_data,
  aes(x = config, y = treat, fill = proportion_agreeing_with_main)
) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(
    aes(label = scales::percent(proportion_agreeing_with_main, accuracy = 1)),
    size = 3.5
  ) +
  scale_fill_viridis_c(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = scales::label_percent(),
    name = "Agreement"
  ) +
  labs(
    title = "Agreement with the forward warm-start source prediction",
    subtitle = "Proportion across blocks and direct starting-kappa values",
    x = "Spatial resolution",
    y = NULL
  ) +
  base_theme +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

ggsave(
  file.path(plot_directory, "agreement_with_main_analysis.png"),
  agreement_plot,
  width = 9,
  height = 4.5,
  units = "in",
  dpi = 300
)

# Within-scenario prediction stability ----------------------------------
prediction_stability_plot <- scenario_stability |>
  format_results() |>
  ggplot(aes(x = config, y = n_distinct_predictions, color = block)) +
  geom_point(
    position = position_jitter(width = 0.12, height = 0),
    size = 2.2,
    alpha = 0.85
  ) +
  facet_wrap(~treat, nrow = 1, scales = "free_x") +
  scale_y_continuous(breaks = seq_len(length(kappa_breaks))) +
  labs(
    title = "Sensitivity of the selected source set to initialization",
    subtitle = paste0(
      "Number of distinct winners across ", length(kappa_breaks),
      " direct starting-kappa values; 1 indicates complete stability"
    ),
    x = "Spatial resolution",
    y = "Distinct predicted source sets",
    color = "Block"
  ) +
  base_theme +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

ggsave(
  file.path(plot_directory, "prediction_stability_across_starts.png"),
  prediction_stability_plot,
  width = 11,
  height = 4.5,
  units = "in",
  dpi = 300
)

# Objective values ------------------------------------------------------
likelihood_plot <- comparison |>
  format_results() |>
  ggplot(aes(
    x = initial_kappa,
    y = observed_loglik_difference,
    group = block,
    color = block
  )) +
  geom_hline(yintercept = 0, color = "grey35", linetype = 2) +
  geom_line(alpha = 0.65, linewidth = 0.6) +
  geom_point(size = 1.8) +
  facet_grid(treat ~ config, scales = "free_y", drop = TRUE) +
  scale_x_log10(breaks = kappa_breaks, labels = scales::label_number()) +
  labs(
    title = "Final likelihood relative to the forward warm-start fit",
    subtitle = "Positive values indicate a higher observed-data likelihood",
    x = expression("Direct starting " * kappa),
    y = "Observed log-likelihood difference",
    color = "Block"
  ) +
  base_theme

ggsave(
  file.path(plot_directory, "likelihood_difference_from_main.png"),
  likelihood_plot,
  width = 13,
  height = 8,
  units = "in",
  dpi = 300
)

# Accuracy sensitivity --------------------------------------------------
accuracy_plot <- comparison |>
  select(
    initial_kappa, block, treat, config,
    standard_accuracy, distance_weighted_accuracy
  ) |>
  pivot_longer(
    cols = c(standard_accuracy, distance_weighted_accuracy),
    names_to = "metric",
    values_to = "accuracy"
  ) |>
  mutate(
    metric = factor(
      metric,
      levels = c("standard_accuracy", "distance_weighted_accuracy"),
      labels = c("Exact accuracy", "Distance-weighted accuracy")
    )
  ) |>
  format_results() |>
  ggplot(aes(
    x = initial_kappa,
    y = accuracy,
    group = block,
    color = block
  )) +
  geom_line(alpha = 0.55, linewidth = 0.6) +
  geom_point(size = 1.7) +
  facet_grid(treat + metric ~ config, drop = TRUE) +
  scale_x_log10(breaks = kappa_breaks, labels = scales::label_number()) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(
    title = "Whole-study source accuracy across direct initializations",
    x = expression("Direct starting " * kappa),
    y = "Accuracy",
    color = "Block"
  ) +
  base_theme

ggsave(
  file.path(plot_directory, "accuracy_across_starts.png"),
  accuracy_plot,
  width = 13,
  height = 12,
  units = "in",
  dpi = 300
)

# Convergence behavior --------------------------------------------------
iteration_plot <- convergence_summary |>
  format_results() |>
  ggplot(aes(
    x = initial_kappa,
    y = iterations,
    group = block,
    color = block,
    shape = converged
  )) +
  geom_line(alpha = 0.5, linewidth = 0.55) +
  geom_point(size = 2) +
  facet_grid(treat ~ config, scales = "free_y", drop = TRUE) +
  scale_x_log10(breaks = kappa_breaks, labels = scales::label_number()) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 4)) +
  labs(
    title = "EM convergence across direct initializations",
    x = expression("Direct starting " * kappa),
    y = "EM iterations",
    color = "Block",
    shape = "Converged"
  ) +
  base_theme

ggsave(
  file.path(plot_directory, "iterations_across_starts.png"),
  iteration_plot,
  width = 13,
  height = 8,
  units = "in",
  dpi = 300
)

# Does the optimized distance parameter retain memory of its start? ------
kappa_plot <- parameters |>
  format_results() |>
  ggplot(aes(
    x = initial_kappa,
    y = kappa,
    group = interaction(block, visit),
    color = factor(visit)
  )) +
  geom_line(alpha = 0.45, linewidth = 0.5) +
  geom_point(alpha = 0.7, size = 1.3) +
  facet_grid(treat ~ config, scales = "free_y", drop = TRUE) +
  scale_x_log10(breaks = kappa_breaks, labels = scales::label_number()) +
  labs(
    title = expression("Optimized " * kappa * " across direct starts"),
    subtitle = "Each line represents one block-transition",
    x = expression("Direct starting " * kappa),
    y = expression("Final " * kappa),
    color = "Transition\nending visit"
  ) +
  base_theme

ggsave(
  file.path(plot_directory, "optimized_kappa_across_starts.png"),
  kappa_plot,
  width = 13,
  height = 8,
  units = "in",
  dpi = 300
)

message("Initialization sensitivity plots saved to: ", plot_directory)
