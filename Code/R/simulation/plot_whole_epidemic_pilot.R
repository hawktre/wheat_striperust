# Diagnostic plots for the local whole-epidemic pilot.

library(dplyr)
library(ggplot2)
library(here)
library(readr)

result_directory <- here(
  "DataProcessed", "results", "simulation", "whole_epidemic_pilot"
)
plot_directory <- file.path(result_directory, "plots")
dir.create(plot_directory, recursive = TRUE, showWarnings = FALSE)

epidemic <- read_csv(file.path(result_directory, "epidemic_summary.csv"), show_col_types = FALSE)
plants <- read_csv(file.path(result_directory, "plant_severity.csv"), show_col_types = FALSE)
oracle <- read_csv(file.path(result_directory, "oracle_source_diagnostics.csv"), show_col_types = FALSE)

mechanism_labels <- c(
  observed = "Observed",
  secondary_spread = "Recursive secondary spread",
  persistent_sources = "Persistent original sources only"
)
epidemic <- epidemic |>
  mutate(mechanism_label = factor(mechanism_labels[mechanism], levels = mechanism_labels))

trajectory_plot <- ggplot(
  epidemic,
  aes(x = as.numeric(visit), y = mean_severity, group = interaction(block, replicate))
) +
  geom_line(aes(color = block), alpha = 0.55) +
  geom_point(aes(color = block), size = 1.5) +
  facet_grid(treat ~ mechanism_label) +
  labs(
    title = "Observed and recursively simulated epidemic trajectories",
    x = "Visit", y = "Mean disease severity", color = "Block"
  ) +
  theme_bw() + theme(legend.position = "bottom")
ggsave(file.path(plot_directory, "mean_severity_trajectories.png"), trajectory_plot,
       width = 11, height = 7, dpi = 300)

zero_plot <- ggplot(
  epidemic,
  aes(x = as.numeric(visit), y = zero_fraction, group = interaction(block, replicate))
) +
  geom_line(aes(color = block), alpha = 0.55) +
  geom_point(aes(color = block), size = 1.5) +
  facet_grid(treat ~ mechanism_label) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Zero-inflation through time",
    x = "Visit", y = "Fraction of plants with zero severity", color = "Block"
  ) +
  theme_bw() + theme(legend.position = "bottom")
ggsave(file.path(plot_directory, "zero_fraction_trajectories.png"), zero_plot,
       width = 11, height = 7, dpi = 300)

# One compact spatial example from each treatment and mechanism.
spatial_data <- plants |>
  filter(block == "A", replicate %in% c(0, 1)) |>
  filter((mechanism == "observed" & replicate == 0) |
           (mechanism != "observed" & replicate == 1)) |>
  mutate(mechanism_label = factor(mechanism_labels[mechanism], levels = mechanism_labels))
spatial_plot <- ggplot(spatial_data, aes(east, north)) +
  geom_tile(aes(fill = severity), color = "grey80") +
  geom_point(data = filter(spatial_data, true_source), shape = 4, size = 2.2) +
  facet_grid(treat + mechanism_label ~ visit) +
  coord_equal() +
  scale_fill_viridis_c(limits = c(0, 1)) +
  labs(
    title = "Spatial development of observed and simulated epidemics",
    subtitle = "Crosses mark the experimentally inoculated plants (block A)",
    x = NULL, y = NULL, fill = "Severity"
  ) +
  theme_bw() + theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave(file.path(plot_directory, "spatial_epidemic_examples.png"), spatial_plot,
       width = 12, height = 15, dpi = 300)

oracle <- oracle |>
  mutate(
    mechanism = factor(mechanism_labels[mechanism], levels = mechanism_labels[-1]),
    visit = factor(visit, ordered = TRUE),
    config = factor(config, levels = c("4", "8h", "8v", "16", "64"))
  )
rank_plot <- ggplot(oracle, aes(visit, true_rank, color = block, group = interaction(block, replicate))) +
  geom_line(alpha = 0.45) + geom_point(alpha = 0.7) +
  facet_grid(treat + mechanism ~ config, scales = "free_y") +
  scale_y_reverse() +
  labs(
    title = "Oracle rank of the true source set",
    subtitle = "Uses the generating transition parameters; rank 1 is best",
    x = "Evidence through visit", y = "True-source rank", color = "Block"
  ) +
  theme_bw() + theme(legend.position = "bottom")
ggsave(file.path(plot_directory, "oracle_true_source_rank.png"), rank_plot,
       width = 13, height = 11, dpi = 300)

entropy_plot <- ggplot(
  oracle,
  aes(visit, effective_candidates, color = block, group = interaction(block, replicate))
) +
  geom_line(alpha = 0.45) + geom_point(alpha = 0.7) +
  facet_grid(treat + mechanism ~ config, scales = "free_y") +
  labs(
    title = "Oracle posterior ambiguity",
    subtitle = "Larger values indicate that evidence is spread over more candidate source sets",
    x = "Evidence through visit", y = "Effective number of candidates", color = "Block"
  ) +
  theme_bw() + theme(legend.position = "bottom")
ggsave(file.path(plot_directory, "oracle_effective_candidates.png"), entropy_plot,
       width = 13, height = 11, dpi = 300)

message("Whole-epidemic pilot plots saved to: ", plot_directory)
