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
library(sf)
library(kableExtra)
library(wesanderson)
# Read in experimental data -----------------------------------------------
## Experiment Data
stripe_2024 <- readRDS(here("DataProcessed/experimental/stripe_clean.rds")) %>%
  st_as_sf(coords = c("east", "north"))

## Inoculations
inocs_2024 <- readRDS(here("DataProcessed/experimental/inoc_sp.rds"))

## Grids
design_2024 <- readRDS(here("DataProcessed/experimental/grids_sp.rds"))


# Read in Results --------------------------------------------------------
forward <- readRDS(here("DataProcessed/results/forward_model/forward_fits.rds"))
backward <- readRDS(here(
  "DataProcessed/results/backward_model/all_backward_fits.rds"
))

# Separate Sensitivity Analyses ------------------------------------------
backward_shared <- backward$shared |>
  select(-c(Q_track, pi, component_dist_acc))
estimation <- rbind(backward$individual, backward$shared)
initialization <- rbind(
  backward$shared,
  backward$sensitivity |> select(names(backward$shared))
)

## Set plotting colors
intensity_cols <- wesanderson::wes_palette("Zissou1", 7, type = c("continuous"))
parameter_cols <- wesanderson::wes_palette("Darjeeling1", 5, "discrete")

# Experimental Design Plot -----------------------------------------------
grids_2024 <- rbind(
  design_2024$`4`$grid,
  design_2024$`8h`$grid,
  design_2024$`8v`$grid,
  design_2024$`16`$grid,
  design_2024$`64`$grid
) %>%
  mutate(
    config = case_when(
      name == "K = 4" ~ "4",
      name == "K = 8 (Horizontal)" ~ "8h",
      name == "K = 8 (Vertical)" ~ "8v",
      name == "K = 16" ~ "16",
      name == "K = 64" ~ "64"
    )
  )

grids_2024 %>%
  mutate(
    name = factor(
      name,
      levels = c(
        "K = 4",
        "K = 8 (Horizontal)",
        "K = 8 (Vertical)",
        "K = 16",
        "K = 64"
      )
    )
  ) %>%
  ggplot() +
  geom_sf(fill = "transparent") +
  geom_sf_text(aes(label = grid_id), size = 2) +
  geom_sf(
    data = inocs_2024 %>% filter(block == "C"),
    aes(shape = "Inoculation"),
    fill = "#CC0000",
    size = 1.5
  ) +
  labs(x = "East (m)", y = "North (m)", fill = "Inoculation") +
  facet_grid(
    treat ~ name,
    labeller = labeller(
      treat = function(x) paste0(x, " inoculation(s)")
    )
  ) +
  scale_shape_manual(values = 23, name = "Legend") +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave(
  filename = "experimental_design_2024_blockC.png",
  path = here("manuscript/Figures/"),
  width = 8,
  height = 6,
  units = "in"
)

# Boxplot Raw Data ---------------------------------------------------------
stripe_2024 %>%
  ggplot(aes(x = visit, y = intensity, fill = treat)) +
  geom_boxplot() +
  facet_wrap(~block) +
  labs(
    x = "Visit",
    y = "Disease Intensity",
    fill = "Inoculation Sites"
  ) +
  scale_fill_manual(values = intensity_cols[c(1, 4, 7)]) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(
  "disease_boxplot.png",
  path = here("manuscript/Figures/"),
  width = 8,
  height = 7,
  units = "in"
)

# Raw Data summary (tabular) -----------------------------------------

stripe_summary <- stripe_2024 %>%
  st_drop_geometry() |>
  group_by(block, treat, visit) %>%
  summarise(
    mean_disease = paste0(
      round(mean(intensity), 3),
      " (",
      round(sd(intensity), 2),
      ")"
    ),
    zeros = paste0(round(sum(intensity == 0) / n(), 3)),
    .groups = "drop"
  ) %>%
  ungroup()


# pivot to wide format
summary_wide <- stripe_summary %>%
  pivot_wider(
    names_from = visit,
    values_from = c(mean_disease, zeros)
  ) %>%
  dplyr::select(-c(zeros_3, zeros_4, zeros_5)) |>
  rename("Block" = block, "Inoculations" = treat)

# make gt table
summary_wide %>%
  select(-Block) |>
  kbl(
    digits = 3,
    booktabs = TRUE,
    format = "latex",
    col.names = c("Inoculations", paste0("Visit ", 1:5), paste0("Visit ", 1:2))
  ) |>
  add_header_above(c(
    " " = 1,
    "Disease Intensity (Prop.)" = 5,
    "Zeros (Prop.)" = 2
  )) |>
  pack_rows("Block A", 1, 3, hline_after = T) |>
  pack_rows("Block B", 4, 6, hline_before = T, hline_after = T) |>
  pack_rows("Block C", 7, 9, hline_before = T, hline_after = T) |>
  pack_rows("Block D", 10, 12, hline_before = T, hline_after = T)


# Forward Model Parameters ------------------------------------------------
params_2024 <- forward %>% 
  select(block, treat, visit, theta) %>% 
  unnest_longer(theta, indices_to = "param", values_to = "value") %>% 
  mutate(value = case_when(param == "phi" ~ exp(value),
                           T ~ value),
         param = case_when(param == "beta" ~ "Intercept",
                           param == "delta" ~ "Auto-Infection",
                           param == "gamma" ~ "Cross-Infection",
                           param == "kappa" ~ "Distance",
                           param == "phi" ~ "Precision"))

params_2024 %>% 
  ggplot(aes(x = as.numeric(visit), y = value))+
  geom_point(aes(color = param))+
  geom_line(aes(color = param))+
  facet_grid(block ~ treat, scales = "free_y", label = labeller(
    treat = function(x) paste0("Inocs: ", x)
  ))+
  scale_color_manual(values = parameter_cols) +
  labs(x = "Visit", y = "Estimate", color = "Parameter")+
  theme(legend.position = "bottom")

ggsave(filename = "parameter_estimates.png", path = here("manuscript/Figures/"), width = 8, height = 5, units = "in")



# Sensitivity Analysis (Plots) --------------------------------------------
acc_diff_estimation <- estimation |>
  filter(!(n_src == 4 & config == "64")) %>% 
  mutate(source_type = factor(if_else(n_src == 1, "Single-source", "Multi-source"), 
                              levels = c("Single-source", "Multi-source")),
    config = factor(
      config,
      levels = c("4", "8h", "8v", "16", "64"),
      labels = c("4", "8 Horizontal", "8 Vertical", "16", "64")
    ),
    visit = factor(visit)) |>
  group_by(config, source_type, visit, result_type) |>
  summarize(
    dist_acc = median(dist_acc),
    .groups = "drop"
  ) |>
  pivot_wider(
    names_from = result_type,
    values_from = dist_acc
  ) |>
  mutate(
    diff = shared - individual
  )


estimation_heat <- ggplot(acc_diff_estimation,
       aes(x = visit, y = config, fill = diff)) +
  geom_tile(color = "black") +
  facet_grid(~ source_type, scales = "free_y", space = "free_y",
             labeller = labeller(
               n_src = function(x) paste("True Sources:", x)
             )) +
  scale_fill_viridis_c(
    na.value = "gray",
    limits = c(-0.5, 0.5),
    breaks = c(-0.5, 0, 0.5),
    labels = c("Specific", "0", "Shared"),
    option = "viridis",
    name = expression(Delta*" Distance-weighted Accuracy")
  ) +
  labs(
    x = "Visit",
    y = "K",
    title = "Estimation: Shared vs. Configuration-specific"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank()
  )

acc_diff_initialization <- initialization |>
  filter(!(n_src == 4 & config == "64")) %>% 
  mutate(source_type = factor(if_else(n_src == 1, "Single-source", "Multi-source"), 
                              levels = c("Single-source", "Multi-source")),
    config = factor(
      config,
      levels = c("4", "8h", "8v", "16", "64"),
      labels = c("4", "8 Horizontal", "8 Vertical", "16", "64")
    ),
    visit = factor(visit)) |>
  group_by(config, source_type, visit, result_type) |>
  summarize(
    dist_acc = median(dist_acc),
    .groups = "drop"
  ) |>
  pivot_wider(
    names_from = result_type,
    values_from = dist_acc
  ) |>
  mutate(
    diff = shared - sensitivity
  )


init_heat <- ggplot(acc_diff_initialization,
       aes(x = visit, y = config, fill = diff)) +
  geom_tile(color = "black") +
  facet_grid(~ source_type, scales = "free_y", space = "free_y",
             labeller = labeller(
               n_src = function(x) paste("True Sources:", x)
             )) +
  scale_fill_viridis_c(
    na.value = "gray",
    limits = c(-0.5, 0.5),
    breaks = c(-0.5, 0, 0.5),
    labels = c("Naive", "0", "Informed"),
    option = "magma",
    name = expression(Delta*" Distance-weighted Accuracy")
  ) +
  labs(
    title = "Initialization: Informed vs. Naive",
    x = "Visit",
    y = ""
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank()
  )

library(ggpubr)

ggarrange(estimation_heat, init_heat, nrow = 1)
ggsave("Sensitivity_heatmap.png", path = here("manuscript/Figures/"), width = 10, height = 4)
# Sensitivity Analysis (Tables) -------------------------------------------
## Estimation
acc_diff_estimation %>% 
  pivot_wider(names_from = source_type, values_from = c(shared, individual, diff)) %>% 
  select(visit, contains("Single"), contains("Multi")) %>% 
  kbl(format = "latex", digits = 3, booktabs = TRUE,
      col.names = c("Visit", rep(c("Shared", "Component-Specific", "Difference"), 2))) %>% 
  add_header_above(c(" " = 1, "Single-source" = 3, "Multi-source" = 3)) %>% 
  pack_rows("K = 4", 1, 4, hline_after = T) %>% 
  pack_rows("K = 8 Horizontal", 5, 8, indent = T, hline_before = T, hline_after = T) %>% 
  pack_rows("K = 8 Vertical", 9, 12, indent = T, hline_before = T, hline_after = T) %>% 
  pack_rows("K = 16", 13, 16, indent = T, hline_before = T, hline_after = T) %>% 
  pack_rows("K = 64", 17, 20, indent = T, hline_before = T, hline_after = T) 

## Initialization
acc_diff_initialization %>% 
  pivot_wider(names_from = source_type, values_from = c(shared, sensitivity, diff)) %>% 
  select(visit, contains("Single"), contains("Multi")) %>% 
  kbl(format = "latex", digits = 3, booktabs = TRUE,
      col.names = c("Visit", rep(c("Informed", "Naive", "Difference"), 2))) %>% 
  add_header_above(c(" " = 1, "Single-source" = 3, "Multi-source" = 3)) %>% 
  pack_rows("K = 4", 1, 4, hline_after = T) %>% 
  pack_rows("K = 8 Horizontal", 5, 8, indent = T, hline_before = T, hline_after = T) %>% 
  pack_rows("K = 8 Vertical", 9, 12, indent = T, hline_before = T, hline_after = T) %>% 
  pack_rows("K = 16", 13, 16, indent = T, hline_before = T, hline_after = T) %>% 
  pack_rows("K = 64", 17, 20, indent = T, hline_before = T, hline_after = T) 
  

# N-sources  -------------------------------------------------------------
backward_shared |>
  filter(treat != 1) %>% 
  mutate(
    config = factor(
      config,
      levels = c("4", "8h", "8v", "16", "64"),
      labels = c("4", "8 Horizontal", "8 Vertical", "16", "64")
    )
  ) |>
  select(config, block, treat, n_src) |>
  distinct() %>% 
  pivot_wider(values_from = n_src, names_from = c(block, treat)) %>% 
  arrange(config) %>% 
  kbl(digits = 3, format = "latex", booktabs = TRUE,
      col.names = c("K", rep(LETTERS[1:4], 2))) %>% 
  add_header_above(c(" " = 1, "2 Inoculations" = 4, "4 Inoculations" = 4)) 


# Overall Results (Figure) ------------------------------------------------
backward_shared |>
  mutate(source_type = if_else(n_src == 1, "Single-source", "Multi-source")) %>% 
  mutate(visit = factor(visit),
         config = factor(
           config,
           levels = c("4", "8h", "8v", "16", "64"),
           labels = c("4", "8 Horizontal", "8 Vertical", "16", "64")
         )
  ) %>% 
  ggplot(aes(x = visit, y = dist_acc))+
  geom_boxplot()+
  facet_grid(source_type~config, 
             labeller = labeller(
               n_src = function(x) paste0(x, " True Source(s)")
             )
  )+
  labs(y = "Distance-weighted Accuracy",
       x = "Visit") +
  theme_bw()
  
ggsave(filename = "experiment_boxplot_main.png", path = here("manuscript/Figures/"), width = 14, height = 6, units = "in")

backward_shared |>
  filter(treat %in% c(2,4)) %>% 
  mutate(visit = factor(visit),
         config = factor(
           config,
           levels = c("4", "8h", "8v", "16", "64"),
           labels = c("4", "8 Horizontal", "8 Vertical", "16", "64")
         )
  ) %>% 
  ggplot(aes(x = visit, y = dist_acc))+
  geom_boxplot()+
  facet_wrap(~config, nrow = 1,
             labeller = labeller(
               n_src = function(x) paste0(x, " True Source(s)")
             )
  )+
  labs(y = "Distance-weighted Accuracy",
       x = "Visit",
       title = "Multi-source") +
  theme_bw()

ggsave(filename = "multi_source_boxplot_main.png", path = here("manuscript/Figures/"), width = 14, height = 4, units = "in")

# Overall Results (Tabular) --------------------------------------------------------
## Averaged within configurations
backward_shared |>
  mutate(source_type = if_else(n_src == 1, "Single-source", "Multi-source")) %>%
  mutate(
    config = factor(
      config,
      levels = c("4", "8h", "8v", "16", "64"),
      labels = c("4", "8 Horizontal", "8 Vertical", "16", "64")
    )
  ) |>
  group_by(config, visit, source_type) |>
  summarise(N = n(),
            dist_acc = median(dist_acc),
    acc = median(acc),
    error = median(mean_error)
  ) |>
  ungroup() |>
  pivot_wider(names_from = source_type, values_from = c(N, dist_acc, acc, error)) %>% 
  select(visit, contains("Single"), contains("Multi")) %>% 
  kbl(format = "latex",
    digits = 3,
    booktabs = TRUE,
    col.names = c("Visit", rep(c("N","Distance-weighted Acc.", "Accuracy", "Error (m)"), 2))
  ) %>% 
  add_header_above(c(" " = 1, "Single-source" = 4, "Multi-source" = 4)) %>% 
  pack_rows("K = 4", 1, 4, indent = T, hline_after = T) |>
  pack_rows("K = 8 (Horizontal)", 5, 8, indent = T, hline_after = T, hline_before = T) |>
  pack_rows("K = 8 (Vertical)", 9, 12, indent = T, hline_after = T, hline_before = T) |>
  pack_rows("K = 16", 13, 16, indent = T, hline_after = T, hline_before = T) |>
  pack_rows("K = 64", 17, 20, indent = T, hline_after = T, hline_before = T) 

## Avereged within True Sources
backward_shared |>
  mutate(
    config = factor(
      config,
      levels = c("4", "8h", "8v", "16", "64"),
      labels = c("4", "8 Horizontal", "8 Vertical", "16", "64")
    )
  ) |>
  group_by(n_src, visit) |>
  summarise(N = n(),
            dist_acc = mean(dist_acc),
            acc = mean(acc),
            error = mean(mean_error)
  ) |>
  ungroup() |>
  select(-c(n_src, N)) %>% 
  kbl(
    digits = 3,
    booktabs = TRUE,
    col.names = c("Visit", "Distance-weighted Acc.", "Accuracy", "Error (m)")
  ) |>
  pack_rows("True Sources: 1 | 22 Predictions" , 1, 4, indent = T, hline_after = T) |>
  pack_rows("True Sources: 2 | 22 Predictions", 5, 8, indent = T, hline_after = T, hline_before = T) |>
  pack_rows("True Sources: 3 | 10 Predictions", 9, 12, indent = T, hline_after = T, hline_before = T) |>
  pack_rows("True Sources: 4 | 6 Predictions", 13, 16, indent = T, hline_after = T, hline_before = T) 


# Simulation Results (forward parameter estimates) ------------------------------------------------------
sims <- readRDS(here("DataProcessed/results/simulation/sims.rds"))

# Summarise convergence results
forward_sims <- sims %>% select(sim, block, treat, visit, converged.forward, neg_loglik, iters, grad_norm, init_kappa, theta.forward) %>% distinct()

greek_cols <- RColorBrewer::brewer.pal(5, "Set1")

forward_sims_long <- forward_sims %>%
  filter(converged.forward == T) %>% 
  # Convert each named theta vector into a tibble (with name = param)
  mutate(theta = map(theta.forward, ~ enframe(.x, name = "param", value = "value"))) %>%
  # Unnest into long format
  unnest(theta) %>% 
  mutate(value = case_when(param == "phi" ~ exp(value),
                           T ~ value))

true_params <- forward %>% 
  mutate(theta = map(theta, ~ enframe(.x, name = "param", value = "value"))) %>%
  # Unnest into long format
  unnest(theta) %>% 
  select(block, treat, visit, param, value)%>% 
  mutate(value = case_when(param == "phi" ~ exp(value),
                           T ~ value))

forward_sims_long <- left_join(forward_sims_long, true_params %>% mutate(visit = as.numeric(visit)), by = c("block", "treat", "visit", "param"), suffix = c(".sim", ".true")) %>% 
  mutate(param = case_when(param == "beta" ~ "Intercept",
                          param == "delta" ~ "Auto-Infection",
                          param == "gamma" ~ "Cross-Infection",
                          param == "kappa" ~ "Distance",
                          param == "phi" ~ "Precision"))

plot_forward_sims <- function(blk, drop_gamma = F){
  if (drop_gamma) {
    forward_sims_long %>%
      mutate(treat_lab = paste0(treat," Inoculation(s)"),
             visit_lab = paste0("Visit ", visit)) %>% 
      filter(block == blk & param != "Cross-Infection") %>%
      ggplot(aes(x = param, y = value.sim)) +
      geom_boxplot(outlier.alpha = 0.3) +
      geom_point(aes(y = value.true, color = param), 
                 shape = 18, size = 3, position = position_dodge(width = 0.75)) +
      facet_grid(visit_lab ~ treat_lab, scales = "free_y") +
      labs(y = "Estimated value", x = "Parameter", color = "True Value",
           title = paste0("Forward Model Simulated Parameter Estimates (Block ", blk, ")")) +
      theme_bw()+
      scale_color_manual(values = parameter_cols[c(1,2,4,5)])+
      theme(legend.position = "none")
  }else{
    forward_sims_long %>%
      mutate(treat_lab = paste0(treat," Inoculation(s)"),
             visit_lab = paste0("Visit ", visit)) %>% 
      filter(block == blk) %>%
      ggplot(aes(x = param, y = value.sim)) +
      geom_boxplot(outlier.alpha = 0.3) +
      geom_point(aes(y = value.true, color = param), 
                 shape = 18, size = 3, position = position_dodge(width = 0.75)) +
      facet_grid(visit_lab ~ treat_lab, scales = "free_y") +
      labs(y = "Estimated value", x = "Parameter", color = "True Value",
           title = paste0("Forward Model Simulated Parameter Estimates (Block ", blk, ")")) +
      theme_bw()+
      scale_color_manual(values = parameter_cols)+
      theme(legend.position = "none")
  }
}

for (blk in LETTERS[1:4]) {
  tmp <- plot_forward_sims(blk)
  ggsave(filename = paste0("forward_sims_bias_blk", blk,".png"), path = here("manuscript/Figures/"), width = 14, height = 6, units = "in")
  
}

for (blk in LETTERS[1:4]) {
  tmp <- plot_forward_sims(blk, drop_gamma = TRUE)
  ggsave(filename = paste0("forward_sims_bias_nogamma_blk", blk,".png"), path = here("manuscript/Figures/"), width = 10, height = 6, units = "in")
}

blkA_simbias <- forward_sims_long %>%
  mutate(treat_lab = paste0(treat," Inoculation(s)"),
         visit_lab = paste0("Visit ", visit)) %>% 
  filter(block == "A") %>%
  ggplot(aes(x = param, y = value.sim)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_point(aes(y = value.true, color = param), 
             shape = 18, size = 3, position = position_dodge(width = 0.75)) +
  facet_grid(visit_lab ~ treat_lab, scales = "free_y") +
  labs(y = "Estimated Value", x = "Parameter", color = "True Value") +
  theme_bw()+
  scale_color_manual(values = parameter_cols)+
  theme(legend.position = "none")

blkA_simbias_nogamma <- forward_sims_long %>%
  mutate(treat_lab = paste0(treat," Inoculation(s)"),
         visit_lab = paste0("Visit ", visit)) %>% 
  filter(block == "A" & param != "Cross-Infection") %>%
  ggplot(aes(x = param, y = value.sim)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_point(aes(y = value.true, color = param), 
             shape = 18, size = 3, position = position_dodge(width = 0.75)) +
  facet_grid(visit_lab ~ treat_lab, scales = "free_y") +
  labs(y = "", x = "Parameter", color = "True Value") +
  theme_bw()+
  scale_color_manual(values = parameter_cols[c(1,3,4,5)])+
  theme(legend.position = "none")

ggsave(blkA_simbias, filename = "parameter_bias_blkA.png", path = here("manuscript/Figures/"), width = 14, height = 6, units = "in")
ggsave(blkA_simbias_nogamma, filename = "parameter_bias_nogamma_blkA.png", path = here("manuscript/Figures/"), width = 14, height = 6, units = "in")

# Bias Table --------------------------------------------------------------
tbl_bias <- forward_sims_long %>% 
  mutate(bias = value.sim - value.true) %>% 
  group_by(block, treat, visit, param) %>% 
  summarise(mean_bias = mean(bias),
            med_bias = median(bias),
            sd_bias = sd(bias), .groups = "drop")

tbl_bias %>% 
  filter(block == "A", param %in% c("Cross-Infection", "Precision")) %>% 
  arrange(desc(abs(med_bias))) %>% 
  select(-block) %>% 
  kbl(format = "latex", 
      digits = 3, 
      booktabs = T,
      col.names = c("Inoculations", "Visit", "Parameter", "Average Bias", "Median Bias", "Std. Dev."))

tbl_bias %>% 
  filter(block == "A") %>% 
  arrange(desc(abs(med_bias))) %>% 
  select(-block) %>% 
  kbl(
      digits = 3, 
      booktabs = T,
      col.names = c("Inoculations", "Visit", "Parameter", "Average Bias", "Median Bias", "Std. Dev."))

# Simulation Prediction ---------------------------------------------------

# Repeat convergence summary for backwards fit
backward_sims <- sims %>% 
  select(sim, config, block, treat, visit, n_src.x, converged.backward, Q_final, em_iters, dist_acc, mean_error, acc, component_dist_acc)

## Summarise not-convered
not_converged <- backward_sims %>% 
  group_by(config, block, treat, visit) %>% 
  summarise(prop_failed = mean(!converged.backward),
            n_failed = sum(!converged.backward)) %>% 
  ungroup() %>% 
  filter(prop_failed > 0) %>% 
  arrange(block, treat, visit)

## Boxplot (n_src)
backward_sims %>% 
  filter(converged.backward == TRUE) %>% 
  mutate(
    visit = factor(visit),
    config = factor(
      config,
      levels = c("4", "8h", "8v", "16", "64"),
      labels = c("4", "8 Horizontal", "8 Vertical", "16", "64")
    )
  ) %>% 
  ggplot(aes(x = visit, dist_acc)) +
  geom_boxplot() +
  facet_grid(config~n_src.x, labeller = labeller(
    n_src.x = function(x) paste0(x, " True Source(s)")
  ))+
  labs(
    x = "Visit",
    y = "Distance-weighted Accuracy"
  ) +
  theme_bw()

ggsave(filename = "sims_distacc.png", path = here("manuscript/Figures/"), width = 10, height = 6)
backward_sims %>% 
  filter(converged.backward == T) %>% 
  mutate(
    config = factor(
      config,
      levels = c("4", "8h", "8v", "16", "64"),
      labels = c("4", "8 Horizontal", "8 Vertical", "16", "64")
    )
  ) |>
  group_by(config, visit) |>
  summarise(dist_acc = mean(dist_acc),
            acc = mean(acc),
            error = mean(mean_error)
  ) |>
  ungroup() |>
  select(-config) %>% 
  kbl(
    digits = 3,
    booktabs = TRUE,
    col.names = c("Visit", "Distance-weighted Acc.", "Accuracy", "Error (m)")
  ) |>
  pack_rows("K = 4", 1, 4, indent = T, hline_after = T) |>
  pack_rows("K = 8 (Horizontal)", 5, 8, indent = T, hline_after = T, hline_before = T) |>
  pack_rows("K = 8 (Vertical)", 9, 12, indent = T, hline_after = T, hline_before = T) |>
  pack_rows("K = 16", 13, 16, indent = T, hline_after = T, hline_before = T) |>
  pack_rows("K = 64", 17, 20, indent = T, hline_after = T, hline_before = T) %>% 
  kable_styling()

