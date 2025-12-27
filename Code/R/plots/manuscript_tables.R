## ---------------------------
##
## Script name: manuscript_tables.R
##
## Purpose of script: Create Professional Tables for Manuscript
##
## Author: Trent VanHawkins
##
## Date Created: 2025-12-26
##
##
## ---------------------------

library(here)
library(tidyverse)
library(gt)
library(kableExtra)

backward <- readRDS(here("DataProcessed/results/backward_model/all_backward_fits.rds"))

estimation <- rbind(backward$individual, backward$shared)

initialization <- rbind(backward$shared |> filter(treat == "1"), backward$sensitivity |> select(names(backward$shared)))

# Sensitivity Analysis For Parameter Estimation --------------------------

estimation |> 
  mutate(config = factor(config, levels = c("4", "8h", "8v", "16", "64")),
visit = as.character(visit)) |> 
  group_by(result_type, config, n_src, visit) |> 
  summarise(N = n(),
error = mean(mean_error),
            acc = mean(acc),
            dist_acc = mean(dist_acc)) |> 
  pivot_wider(names_from = result_type, values_from = c("error", "acc", "dist_acc")) |> 
  ungroup() |>
  mutate(group_label = paste0("Config: ", config, " | Sources: ", n_src, " | N = ", N)) |>
  select(group_label, visit, contains("dist_acc"), contains("acc"), contains("error")) |> 
  gt(groupname_col = "group_label") |> 
  cols_label(
    dist_acc_individual = "Individual",
    dist_acc_shared = "Shared",
    acc_individual = "Individual",
    acc_shared = "Shared",
    error_individual = "Individual",
    error_shared = "Shared",
    visit = "Visit"
  ) |>
  tab_spanner(columns = contains("dist_acc"), label = "Distance-Weighted Acc.") |> 
  tab_spanner(columns = starts_with("acc"), label = "Accuracy") |> 
  tab_spanner(columns = contains("error"), label = "Error (m)") |> 
  fmt_number(decimals = 3)

initialization |> 
  mutate(config = factor(config, levels = c("4", "8h", "8v", "16", "64")),
visit = as.character(visit)) |> 
  group_by(result_type, config, n_src, visit) |> 
  summarise(N = n(),
error = mean(mean_error),
            acc = mean(acc),
            dist_acc = mean(dist_acc)) |> 
  pivot_wider(names_from = result_type, values_from = c("error", "acc", "dist_acc")) |> 
  ungroup() |>
  mutate(group_label = paste0("Config: ", config, " | Sources: ", n_src, " | N = ", N)) |>
  select(group_label, visit, contains("dist_acc"), contains("acc"), contains("error")) |> 
  gt(groupname_col = "group_label") |> 
  cols_label(
    dist_acc_sensitivity = "MOM",
    dist_acc_shared = "MLE",
    acc_sensitivity = "MOM",
    acc_shared = "MLE",
    error_sensitivity = "MOM",
    error_shared = "MLE",
    visit = "Visit"
  ) |>
  tab_spanner(columns = contains("dist_acc"), label = "Distance-Weighted Acc.") |> 
  tab_spanner(columns = starts_with("acc"), label = "Accuracy") |> 
  tab_spanner(columns = contains("error"), label = "Error (m)") |> 
  fmt_number(decimals = 3)
