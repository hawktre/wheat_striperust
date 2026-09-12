# Diagnostic and results plots for the primary source-detection analysis.
# Run after src/source_detection/run_main_analysis.R has completed.

library(dplyr)
library(ggplot2)
library(here)
library(readr)
library(sf)
library(tidyr)

result_directory <- here("output", "source_detection")
plot_directory <- file.path(result_directory, "plots")
dir.create(plot_directory, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  accuracy = "source_accuracy.csv",
  posteriors = "all_primary_fit_posteriors.csv",
  parameters = "em_parameters.csv",
  history = "em_convergence_history.csv"
)
required_paths <- setNames(
  file.path(result_directory, unname(required_files)),
  names(required_files)
)
missing_files <- required_files[!file.exists(required_paths)]
if (length(missing_files)) {
  stop(
    "Run run_main_analysis.R first. Missing result files: ",
    paste(missing_files, collapse = ", ")
  )
}

accuracy <- read_csv(required_paths[["accuracy"]], show_col_types = FALSE)
posteriors <- read_csv(required_paths[["posteriors"]], show_col_types = FALSE)
parameters <- read_csv(required_paths[["parameters"]], show_col_types = FALSE)
em_history <- read_csv(required_paths[["history"]], show_col_types = FALSE)

config_levels <- c("4", "8h", "8v", "16", "64")
config_labels <- c(
  "4 cells", "8 horizontal cells", "8 vertical cells",
  "16 cells", "64 plants"
)
treat_levels <- c("1", "2", "4")
treat_labels <- c(
  "1 planted source", "2 planted sources", "4 planted sources"
)

format_results <- function(data) {
  data |>
    mutate(
      treat_id = as.character(treat),
      config_id = as.character(config),
      block = factor(block),
      treat = factor(as.character(treat), levels = treat_levels,
                     labels = treat_labels),
      config = factor(config, levels = config_levels, labels = config_labels),
      visit = factor(as.character(visit), levels = c("2", "3", "4", "5"),
                     ordered = TRUE),
      fit_type = case_when(
        fit_scope == "transition_refit" ~ "Independent transition",
        fit_scope %in% c("cumulative_refit", "whole_study") ~
          "Cumulative through visit",
        TRUE ~ fit_scope
      ),
      fit_type = factor(
        fit_type,
        levels = c("Independent transition", "Cumulative through visit")
      )
    )
}

accuracy <- format_results(accuracy)

base_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "bottom"
  )

# Accuracy trajectories --------------------------------------------------
# All four blocks are shown directly; no across-block trend is superimposed.
for (treatment_name in levels(accuracy$treat)) {
  accuracy_subset <- accuracy |>
    filter(treat == treatment_name)

  distance_plot <- ggplot(
    accuracy_subset,
    aes(
      x = visit,
      y = distance_weighted_accuracy,
      group = block,
      color = block
    )
  ) +
    geom_line(alpha = 0.55, linewidth = 0.6) +
    geom_point(size = 2) +
    facet_grid(fit_type ~ config, drop = TRUE) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(
      title = treatment_name,
      subtitle = "Each line represents one experimental block",
      x = "Endpoint visit",
      y = "Distance-weighted accuracy",
      color = "Block"
    ) +
    base_theme

  ggsave(
    file.path(
      plot_directory,
      paste0("distance_weighted_accuracy_treatment_", sub(" .*", "", treatment_name), ".png")
    ),
    distance_plot,
    width = 12,
    height = 5.5,
    units = "in",
    dpi = 300
  )
}

# Exact accuracy stays on its discrete support; point area shows overlapping
# block observations.
exact_accuracy_plot <- ggplot(
  accuracy,
  aes(x = visit, y = standard_accuracy, color = fit_type)
) +
  geom_count(position = position_dodge(width = 0.45)) +
  facet_grid(treat ~ config, drop = TRUE) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_size_area(max_size = 6, breaks = 1:4) +
  labs(
    title = "Exact source-identification accuracy",
    subtitle = "Point area represents the number of overlapping blocks",
    x = "Endpoint visit",
    y = "Proportion of source cells identified exactly",
    color = "Fit",
    size = "Blocks"
  ) +
  base_theme

ggsave(
  file.path(plot_directory, "exact_accuracy_all_scenarios.png"),
  exact_accuracy_plot,
  width = 13,
  height = 8,
  units = "in",
  dpi = 300
)

# Heatmaps preserve the discrete observed values and make differences across
# visits, blocks, resolutions, and temporal fit scopes easy to scan.
exact_accuracy_heatmap <- ggplot(
  accuracy,
  aes(x = visit, y = block, fill = standard_accuracy)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = format(round(standard_accuracy, 2), nsmall = 2)),
            size = 2.5) +
  facet_grid(treat + fit_type ~ config, drop = TRUE) +
  scale_fill_viridis_c(limits = c(0, 1), option = "C", direction = -1) +
  labs(
    title = "Exact source-identification accuracy",
    subtitle = "Each tile is one observed block-level result",
    x = "Endpoint visit",
    y = "Block",
    fill = "Exact\naccuracy"
  ) +
  base_theme

ggsave(
  file.path(plot_directory, "exact_accuracy_heatmap.png"),
  exact_accuracy_heatmap,
  width = 13,
  height = 10,
  units = "in",
  dpi = 300
)

# Boxplots are retained as compact distribution summaries, with every one of
# the four contributing block values overlaid so their discreteness is clear.
distance_boxplot <- ggplot(
  accuracy,
  aes(x = visit, y = distance_weighted_accuracy, fill = fit_type)
) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.6,
    alpha = 0.45,
    outlier.shape = NA
  ) +
  geom_point(
    aes(color = block),
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.75
    ),
    size = 1.7
  ) +
  facet_grid(treat ~ config, drop = TRUE) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(
    title = "Distance-weighted source-identification accuracy",
    subtitle = "Boxplots summarize four blocks; points show every observation",
    x = "Endpoint visit",
    y = "Distance-weighted accuracy",
    fill = "Fit",
    color = "Block"
  ) +
  base_theme

ggsave(
  file.path(plot_directory, "distance_weighted_accuracy_boxplots.png"),
  distance_boxplot,
  width = 13,
  height = 8,
  units = "in",
  dpi = 300
)

# Model versus naive maximum-severity prediction ------------------------
accuracy_comparison <- accuracy |>
  select(
    block, treat, config, visit, fit_type,
    model_exact = standard_accuracy,
    naive_exact = naive_standard_accuracy,
    model_distance = distance_weighted_accuracy,
    naive_distance = naive_distance_weighted_accuracy
  ) |>
  pivot_longer(
    cols = c(model_exact, naive_exact, model_distance, naive_distance),
    names_to = c("method", "metric"),
    names_sep = "_",
    values_to = "accuracy"
  ) |>
  mutate(
    method = factor(
      method,
      levels = c("model", "naive"),
      labels = c("Source-detection model", "Maximum severity")
    ),
    metric = factor(
      metric,
      levels = c("exact", "distance"),
      labels = c("Exact accuracy", "Distance-weighted accuracy")
    )
  )

model_naive_plot <- ggplot(
  accuracy_comparison,
  aes(x = visit, y = accuracy, group = interaction(block, method), color = method)
) +
  geom_line(alpha = 0.4, linewidth = 0.5) +
  geom_point(alpha = 0.7, size = 1.5) +
  facet_grid(treat + fit_type + metric ~ config, drop = TRUE) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(
    title = "Source-detection model versus maximum severity",
    subtitle = "Each line represents one experimental block",
    x = "Endpoint visit",
    y = "Accuracy",
    color = "Prediction method"
  ) +
  base_theme

ggsave(
  file.path(plot_directory, "model_vs_naive_accuracy.png"),
  model_naive_plot,
  width = 14,
  height = 16,
  units = "in",
  dpi = 300
)

accuracy_difference <- accuracy |>
  mutate(
    exact_difference = standard_accuracy - naive_standard_accuracy,
    distance_difference =
      distance_weighted_accuracy - naive_distance_weighted_accuracy
  ) |>
  select(
    block, treat, config, visit, fit_type,
    `Exact accuracy` = exact_difference,
    `Distance-weighted accuracy` = distance_difference
  ) |>
  pivot_longer(
    cols = c(`Exact accuracy`, `Distance-weighted accuracy`),
    names_to = "metric",
    values_to = "difference"
  )

accuracy_difference_heatmap <- ggplot(
  accuracy_difference,
  aes(x = visit, y = block, fill = difference)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  facet_grid(treat + fit_type + metric ~ config, drop = TRUE) +
  scale_fill_gradient2(
    low = "#B2182B",
    mid = "white",
    high = "#2166AC",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  labs(
    title = "Source-detection model minus maximum-severity accuracy",
    subtitle = "Blue favors the source-detection model; red favors maximum severity",
    x = "Endpoint visit",
    y = "Block",
    fill = "Accuracy\ndifference"
  ) +
  base_theme

ggsave(
  file.path(plot_directory, "model_minus_naive_accuracy_heatmap.png"),
  accuracy_difference_heatmap,
  width = 14,
  height = 16,
  units = "in",
  dpi = 300
)

# Whole-study spatial predictions ---------------------------------------
# Spatial panels show the actual discrete source choices underlying the
# accuracy scores. One figure is written for each spatial resolution.
spatial_grids <- readRDS(here(
  "data", "processed", "experimental", "grids_sp.rds"
))
whole_study_accuracy <- accuracy |>
  filter(fit_scope == "whole_study")

spatial_plot_parts <- vector("list", nrow(whole_study_accuracy))
for (row_index in seq_len(nrow(whole_study_accuracy))) {
  result_row <- whole_study_accuracy[row_index, ]
  resolution <- result_row$config_id[[1]]
  predicted_sources <- strsplit(
    as.character(result_row$predicted_candidate_id[[1]]), "_", fixed = TRUE
  )[[1]]
  true_sources <- strsplit(
    as.character(result_row$true_candidate_id[[1]]), "_", fixed = TRUE
  )[[1]]

  spatial_plot_parts[[row_index]] <- spatial_grids[[resolution]]$grid |>
    mutate(
      grid_id = as.character(grid_id),
      block = result_row$block[[1]],
      treat = result_row$treat[[1]],
      config = result_row$config[[1]],
      config_id = resolution,
      source_status = case_when(
        grid_id %in% predicted_sources & grid_id %in% true_sources ~
          "Predicted and true",
        grid_id %in% predicted_sources ~ "Predicted only",
        grid_id %in% true_sources ~ "True only",
        TRUE ~ "Neither"
      )
    )
}
spatial_plot_data <- do.call(rbind, spatial_plot_parts) |>
  mutate(
    source_status = factor(
      source_status,
      levels = c(
        "Neither", "Predicted only", "True only", "Predicted and true"
      )
    )
  )

for (resolution_name in levels(accuracy$config)) {
  resolution_data <- spatial_plot_data |>
    filter(config == resolution_name)
  if (nrow(resolution_data) == 0) next

  spatial_source_plot <- ggplot(resolution_data) +
    geom_sf(aes(fill = source_status), color = "grey35", linewidth = 0.3) +
    geom_sf_text(
      data = resolution_data |> filter(source_status != "Neither"),
      aes(label = grid_id),
      size = 2.7
    ) +
    facet_grid(treat ~ block) +
    scale_fill_manual(values = c(
      "Neither" = "white",
      "Predicted only" = "#2166AC",
      "True only" = "#B2182B",
      "Predicted and true" = "#7B3294"
    )) +
    labs(
      title = paste("Whole-study source predictions:", resolution_name),
      subtitle = "Labels identify cells selected as predicted, true, or both",
      x = "East",
      y = "North",
      fill = "Source status"
    ) +
    coord_sf(expand = FALSE) +
    base_theme

  ggsave(
    file.path(
      plot_directory,
      paste0(
        "whole_study_spatial_predictions_",
        resolution_data$config_id[[1]],
        ".png"
      )
    ),
    spatial_source_plot,
    width = 10,
    height = 7.5,
    units = "in",
    dpi = 300
  )
}

model_naive_summary <- accuracy |>
  group_by(treat, config, fit_type, visit) |>
  summarise(
    n_blocks = n(),
    model_exact_accuracy = mean(standard_accuracy),
    naive_exact_accuracy = mean(naive_standard_accuracy),
    exact_accuracy_difference =
      model_exact_accuracy - naive_exact_accuracy,
    model_distance_accuracy = mean(distance_weighted_accuracy),
    naive_distance_accuracy = mean(naive_distance_weighted_accuracy),
    distance_accuracy_difference =
      model_distance_accuracy - naive_distance_accuracy,
    .groups = "drop"
  )

write_csv(
  model_naive_summary,
  file.path(result_directory, "model_vs_naive_accuracy_summary.csv")
)

# Posterior concentration ------------------------------------------------
scenario_columns <- c(
  "block", "treat", "config", "n_sources", "through_visit", "fit_scope"
)
posterior_diagnostics <- posteriors |>
  group_by(across(all_of(scenario_columns))) |>
  summarise(
    winner_posterior = max(posterior),
    true_posterior = posterior[true_candidate][1],
    true_rank = rank[true_candidate][1],
    loglik_gap_from_winner = loglik[true_candidate][1] - max(loglik),
    posterior_entropy = -sum(if_else(
      posterior > 0,
      posterior * log(posterior),
      0
    )),
    effective_candidates = exp(posterior_entropy),
    .groups = "drop"
  ) |>
  mutate(visit = through_visit) |>
  format_results()

write_csv(
  posterior_diagnostics,
  file.path(result_directory, "posterior_diagnostics.csv")
)

true_posterior_plot <- ggplot(
  posterior_diagnostics,
  aes(x = visit, y = true_posterior, group = block, color = block)
) +
  geom_line(alpha = 0.55) +
  geom_point(size = 1.8) +
  facet_grid(treat + fit_type ~ config, drop = TRUE) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(
    title = "Posterior probability assigned to the true source set",
    x = "Endpoint visit",
    y = "True-candidate posterior",
    color = "Block"
  ) +
  base_theme

ggsave(
  file.path(plot_directory, "true_candidate_posterior.png"),
  true_posterior_plot,
  width = 13,
  height = 11,
  units = "in",
  dpi = 300
)

true_rank_plot <- ggplot(
  posterior_diagnostics,
  aes(x = visit, y = true_rank, group = block, color = block)
) +
  geom_line(alpha = 0.55) +
  geom_point(size = 1.8) +
  facet_grid(treat + fit_type ~ config, scales = "free_y", drop = TRUE) +
  scale_y_reverse() +
  labs(
    title = "Rank of the true source set",
    subtitle = "Rank 1 is best",
    x = "Endpoint visit",
    y = "True-candidate rank",
    color = "Block"
  ) +
  base_theme

ggsave(
  file.path(plot_directory, "true_candidate_rank.png"),
  true_rank_plot,
  width = 13,
  height = 11,
  units = "in",
  dpi = 300
)

# EM convergence ---------------------------------------------------------
convergence_summary <- em_history |>
  group_by(
    block, treat, config, n_sources, through_visit, fit_scope
  ) |>
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
  ) |>
  mutate(visit = through_visit) |>
  format_results()

write_csv(
  convergence_summary,
  file.path(result_directory, "convergence_diagnostics.csv")
)

iteration_plot <- ggplot(
  convergence_summary,
  aes(
    x = visit,
    y = iterations,
    color = treat,
    shape = converged
  )
) +
  geom_point(
    position = position_jitter(width = 0.12, height = 0),
    alpha = 0.8,
    size = 2
  ) +
  facet_grid(fit_type ~ config, scales = "free_y", drop = TRUE) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 4)) +
  labs(
    title = "EM iterations and convergence",
    x = "Endpoint visit",
    y = "EM iterations",
    color = "Treatment",
    shape = "Converged"
  ) +
  base_theme

ggsave(
  file.path(plot_directory, "em_iterations_and_convergence.png"),
  iteration_plot,
  width = 12,
  height = 6,
  units = "in",
  dpi = 300
)

# Whole-study parameter estimates ---------------------------------------
whole_study_parameters <- parameters |>
  filter(fit_scope == "whole_study") |>
  format_results() |>
  pivot_longer(
    cols = c(beta, delta, gamma, kappa, phi),
    names_to = "parameter",
    values_to = "estimate"
  ) |>
  mutate(
    parameter = factor(
      parameter,
      levels = c("beta", "delta", "gamma", "kappa", "phi"),
      labels = c("beta", "delta", "gamma", "kappa", "phi")
    )
  )

parameter_plot <- ggplot(
  whole_study_parameters,
  aes(
    x = visit,
    y = estimate,
    group = interaction(block, config),
    color = config
  )
) +
  geom_line(alpha = 0.55) +
  geom_point(alpha = 0.75, size = 1.5) +
  facet_grid(parameter ~ treat, scales = "free_y") +
  labs(
    title = "Whole-study transition-specific parameter estimates",
    x = "Transition ending at visit",
    y = "Estimate",
    color = "Resolution"
  ) +
  base_theme

ggsave(
  file.path(plot_directory, "whole_study_parameter_estimates.png"),
  parameter_plot,
  width = 11,
  height = 10,
  units = "in",
  dpi = 300
)

message("Plots and diagnostic summaries saved to: ", plot_directory)
