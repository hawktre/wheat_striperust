## ---------------------------
##
## Script name: presentation_plots.R
##
## Purpose of script: Create plots for presentation
##
## Author: Trent VanHawkins
##
## Date Created: 2025-11-02
##
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## ---------------------------

## view outputs in non-scientific notation

options(scipen = 6, digits = 4) 

## ---------------------------

## load up the packages we will need:  (uncomment as required)

library(tidyverse)
library(here)
library(sf)
library(kableExtra)

# Read in experimental data -----------------------------------------------
## Experiment Data
stripe_2024 <- readRDS(here("DataProcessed/experimental/stripe_clean.rds")) %>% st_as_sf(coords = c("east", "north"))
stripe_2025 <- readRDS(here("DataProcessed/experimental/2025/stripe_clean_2025.rds")) %>% st_as_sf(coords = c("east", "north"))

## Inoculations
inocs_2024 <- readRDS(here("DataProcessed/experimental/inoc_sp.rds"))
inocs_2025 <- readRDS(here("DataProcessed/experimental/2025/inocs_2025.rds")) %>% st_as_sf(coords = c("east", "north"))

## Grids
design_2024 <- readRDS(here("DataProcessed/experimental/grids_sp.rds"))
design_2025 <- readRDS(here("DataProcessed/experimental/2025/grids_sp_2025.rds"))

# Read in results ---------------------------------------------------------
## Forward Model
forward_2024 <- readRDS(here("DataProcessed/results/forward_model/forward_fits.rds"))
forward_2025 <- readRDS(here("DataProcessed/results/forward_model/forward_fits_2025.rds"))

## Backward Model
backward_2024 <- readRDS(here("DataProcessed/results/backward_model/backward_fits.rds"))
backward_2025 <- readRDS(here("DataProcessed/results/backward_model/backward_fits_2025.rds"))

## Set plotting colors
intensity_cols <- wesanderson::wes_palette("Zissou1", 7, type = c("continuous"))

# Threshold plot -------------------------------------------------------
stripe_2024 %>% 
  mutate(diseased = if_else(intensity > 0, "Yes", "No")) %>% 
  filter(block == "A") %>% 
  ggplot()+
  geom_sf(aes(color = diseased), alpha = 0.75) +
  geom_sf(data = inocs_2024 %>% filter(block == "A"), 
          aes(fill= "Inoculation Point"), 
          shape = 23, size = 2)+
  facet_grid(treat ~ visit, 
             labeller = labeller(
    treat = function(x) paste0("Inocs: ", x),
    visit = function(x) paste0("Visit ", x)
  ))+
  theme_bw() +
  scale_color_manual(values = intensity_cols[c(1,5)]) +
  labs(x = "East (m)", y = "North (m)", color = "Diseased (> 0%)", fill = "") +
  theme(legend.position = "bottom")

ggsave("disease_status.png", path = here("Reports/EEID_Presentation/figures/"), width = 6, height = 8, units = "in")


# True vs. Predicted Plots ------------------------------------------------
## 2024
forward_2024_plt <- forward_2024 %>% 
  select(block, treat, visit, fitted) %>% 
  unnest_longer(fitted, values_to = "fitted", indices_to = "plant_id") %>%
  mutate(plant_id = as.numeric(plant_id)) %>% 
  left_join(stripe_2024 %>% select(block, treat, visit, plant_id, intensity), by = c("block", "treat", "visit", "plant_id")) %>% 
  select(block, treat, visit, plant_id, fitted, intensity, geometry) %>% 
  rename("Predicted" = "fitted",
         "True" = "intensity") %>% 
  pivot_longer(c("Predicted", "True"), names_to = "type", values_to = "intensity")

for (blk in unique(forward_2024_plt$block)) {
  for (trt in unique(forward_2024_plt$treat)) {
    forward_2024_plt %>% 
      filter(block == blk, treat == trt) %>% 
      st_as_sf() %>% 
      ggplot()+
      geom_sf(aes(color = intensity))+
      geom_sf(data = inocs_2024 %>% filter(block == blk, treat == trt), 
              aes(fill= "Inoculation Point"), 
              shape = 23, size = 2, alpha = 0.75)+
      facet_grid(type ~ visit,
                 labeller = labeller(
                   visit = function(x) paste0("Visit ", x)
                 ))+
      theme_bw()+
      labs(x = "East (m)", y = "North (m)", color = "Intensity", fill = "", title = paste0("2024 Fitted Values (Block: ", blk, " Treat: ", trt,")"))+
      scale_colour_gradientn(colours = intensity_cols, limits = c(0,1)) +
      theme(legend.position = "bottom")
    
    ggsave(paste0("forward_preds_2024_",blk,"_",trt,".png"), path = here("Reports/EEID_Presentation/figures/forward/preds"), width = 8, height = 5, units = "in")
  }
}

## 2025
forward_2025_plt <- forward_2025 %>% 
  select(block, treat, visit, fitted) %>% 
  unnest_longer(fitted, values_to = "fitted", indices_to = "plant_id") %>%
  mutate(plant_id = as.numeric(plant_id)) %>% 
  left_join(stripe_2025 %>% select(block, treat, visit, plant_id, intensity), by = c("block", "treat", "visit", "plant_id")) %>% 
  select(block, treat, visit, plant_id, fitted, intensity, geometry) %>% 
  rename("Predicted" = "fitted",
         "True" = "intensity") %>% 
  pivot_longer(c("Predicted", "True"), names_to = "type", values_to = "intensity")

for (blk in unique(forward_2025_plt$block)) {
  for (trt in unique(forward_2025_plt$treat)) {
    forward_2025_plt %>% 
      filter(block == blk, treat == trt) %>% 
      st_as_sf() %>% 
      ggplot()+
      geom_sf(aes(color = intensity))+
      geom_sf(data = inocs_2025 %>% filter(block == blk, treat == trt), 
              aes(fill= "Inoculation Point"), 
              shape = 23, size = 2, alpha = 0.75)+
      facet_grid(type ~ visit,
                 labeller = labeller(
                   visit = function(x) paste0("Visit ", x)
                 ))+
      theme_bw()+
      labs(x = "East (m)", y = "North (m)", color = "Intensity", fill = "", title = paste0("2025 Fitted Values (Block: ", blk, " Treat: ", trt,")"))+
      scale_colour_gradientn(colours = intensity_cols, limits = c(0,1)) +
      theme(legend.position = "bottom")
    
    ggsave(paste0("forward_preds_2025_",blk,"_",trt,".png"), path = here("Reports/EEID_Presentation/figures/forward/preds"), width = 8, height = 5, units = "in")
  }
}



# Residuals ---------------------------------------------------------------
## 2024
forward_2024_resid <- forward_2024 %>%
  mutate(data = map2(fitted, resid, ~tibble(
    plant_id = seq_along(.x),
    fitted = .x,
    resid = .y
  ))) %>%
  select(block, treat, visit, data) %>%
  unnest(data) %>% 
  left_join(stripe_2024, by = c("block", "treat", "visit", "plant_id")) %>% 
  select(block, treat, visit, plant_id, intensity, fitted, resid, geometry) %>% 
  st_as_sf()

for (blk in unique(forward_2024_resid$block)) {
  for (trt in unique(forward_2024_plt$treat)) {
    forward_2024_resid %>% 
      filter(block == blk, treat == trt) %>% 
      ggplot(aes(x = fitted, y = resid)) +
      geom_point()+
      geom_hline(yintercept = 0, linetype = "dashed")+
      facet_wrap(~visit, nrow = 1, 
                 label = labeller(
        visit = function(x) paste0("Visit ", x)
      )) +
      labs(x = "Fitted", y = "Dev. Residual", title = paste0("2024 Residuals (Block: ", blk, " Treat: ", trt,")"))+ 
      lims(x = c(0,1))+
      theme_bw()

    ggsave(paste0("forward_resid_2024_",blk, "_", trt, ".png"), path = here("Reports/EEID_Presentation/figures/forward/residual/"), width = 6, height = 3, units = "in")
  }
}
## 2025
forward_2025_resid <- forward_2025 %>%
  mutate(data = map2(fitted, resid, ~tibble(
    plant_id = seq_along(.x),
    fitted = .x,
    resid = .y
  ))) %>%
  select(block, treat, visit, data) %>%
  unnest(data) %>% 
  left_join(stripe_2025, by = c("block", "treat", "visit", "plant_id")) %>% 
  select(block, treat, visit, plant_id, intensity, fitted, resid, geometry) %>% 
  st_as_sf()

for (blk in unique(forward_2025_resid$block)) {
  for (trt in unique(forward_2025_plt$treat)) {
    forward_2025_resid %>% 
      filter(block == blk, treat == trt) %>% 
      ggplot(aes(x = fitted, y = resid)) +
      geom_point()+
      geom_hline(yintercept = 0, linetype = "dashed")+
      facet_wrap(~visit, nrow = 1, 
                 label = labeller(
                   visit = function(x) paste0("Visit ", x)
                 )) +
      labs(x = "Fitted", y = "Dev. Residual", title = paste0("2025 Residuals (Block: ", blk, " Treat: ", trt,")"))+ 
      lims(x = c(0,1))+
      theme_bw()
    
    ggsave(paste0("forward_resid_2025_",blk, "_", trt, ".png"), path = here("Reports/EEID_Presentation/figures/forward/residual/"), width = 6, height = 3, units = "in")
  }
}


# Spatial Residuals -------------------------------------------------------
## 2024
for (blk in unique(forward_2024_plt$block)) {
  for (trt in unique(forward_2024_plt$treat)) {
    forward_2024_resid %>% 
      filter(block == blk, treat == trt) %>% 
      st_as_sf() %>% 
      ggplot()+
      geom_sf(aes(color = resid))+
      geom_sf(data = inocs_2024 %>% filter(block == blk, treat == trt), 
              aes(fill= "Inoculation Point"), 
              shape = 23, size = 2, alpha = 0.75)+
      geom_sf_label(data = forward_2024_resid %>% 
                      filter(block == blk, treat == trt, abs(resid) > 3), aes(label = plant_id), size = 2, nudge_x = 1.5, nudge_y = -1.5)+
      facet_wrap(~ visit, nrow = 1,
                 labeller = labeller(
                   visit = function(x) paste0("Visit ", x)
                 ))+
      theme_bw()+
      labs(x = "East (m)", y = "North (m)", color = "Intensity", fill = "", title = paste0("2024 Spatial Residuals (Block: ", blk, " Treat: ", trt,")"))+
      scale_colour_gradient(low = intensity_cols[1], high = intensity_cols[length(intensity_cols)], limits = c(-3, 3)) +
      theme(legend.position = "bottom")
    
    ggsave(paste0("forward_spat_resid_2024_",blk,"_",trt,".png"), path = here("Reports/EEID_Presentation/figures/forward/residual/"), width = 8, height = 5, units = "in")
  }
}

## 2025
for (blk in unique(forward_2025_plt$block)) {
  for (trt in unique(forward_2025_plt$treat)) {
    forward_2025_resid %>% 
      filter(block == blk, treat == trt) %>% 
      st_as_sf() %>% 
      ggplot()+
      geom_sf(aes(color = resid))+
      geom_sf(data = inocs_2024 %>% filter(block == blk, treat == trt), 
              aes(fill= "Inoculation Point"), 
              shape = 23, size = 2, alpha = 0.75)+
      geom_sf_label(data = forward_2024_resid %>% 
                      filter(block == blk, treat == trt, abs(resid) > 3), aes(label = plant_id), size = 2, nudge_x = 1, nudge_y = -1)+
      facet_wrap(~ visit, nrow = 1,
                 labeller = labeller(
                   visit = function(x) paste0("Visit ", x)
                 ))+
      theme_bw()+
      labs(x = "East (m)", y = "North (m)", color = "Intensity", fill = "", title = paste0("2025 Spatial Residuals (Block: ", blk, " Treat: ", trt,")"))+
      scale_colour_gradient(low = intensity_cols[1], high = intensity_cols[length(intensity_cols)], limits = c(-3, 3)) +
      theme(legend.position = "bottom")
    
    ggsave(paste0("forward_spat_resid_2025_",blk,"_",trt,".png"), path = here("Reports/EEID_Presentation/figures/forward/residual/"), width = 8, height = 5, units = "in")
  }
}


# Parameter Estimates -----------------------------------------------------
parameter_cols <- wesanderson::wes_palette("Darjeeling1", 5, "discrete")

## 2024
params_2024 <- forward_2024 %>% 
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
  labs(x = "Visit", y = "Estimate", color = "Parameter", title = "2024 Parameter Estimates")+
  theme(legend.position = "bottom")

ggsave(filename = "parameter_estimates_2024.png", path = here("Reports/EEID_Presentation/figures/forward/"), width = 8, height = 5, units = "in")

## 2025
params_2025 <- forward_2025 %>% 
  select(block, treat, visit, theta) %>% 
  unnest_longer(theta, indices_to = "param", values_to = "value") %>% 
  mutate(value = case_when(param == "phi" ~ exp(value),
                           T ~ value),
         param = case_when(param == "beta" ~ "Intercept",
                           param == "delta" ~ "Auto-Infection",
                           param == "gamma" ~ "Cross-Infection",
                           param == "kappa" ~ "Distance",
                           param == "phi" ~ "Precision"))

params_2025 %>% 
  ggplot(aes(x = as.numeric(visit), y = value))+
  geom_point(aes(color = param))+
  geom_line(aes(color = param))+
  facet_grid(block ~ treat, scales = "free_y", label = labeller(
    treat = function(x) paste0("Inocs: ", x)
  ))+
  scale_color_manual(values = parameter_cols) +
  labs(x = "Visit", y = "Estimate", color = "Parameter", title = "2025 Parameter Estimates")+
  theme(legend.position = "bottom")

ggsave(filename = "parameter_estimates_2025.png", path = here("Reports/EEID_Presentation/figures/forward/"), width = 8, height = 5, units = "in")


# Experimental Design Plots -----------------------------------------------
grids_2024 <- rbind(design_2024$`4`$grid,
                    design_2024$`8h`$grid,
                    design_2024$`8v`$grid,
                    design_2024$`16`$grid,
                    design_2024$`64`$grid) %>% mutate(config = case_when(name == "K = 4" ~ "4",
                                                                         name == "K = 8 (Horizontal)" ~ "8h",
                                                                         name == "K = 8 (Vertical)" ~ "8v",
                                                                         name == "K = 16" ~ "16",
                                                                         name == "K = 64" ~ "64"))

grids_2024 %>% 
  mutate(name = factor(name, levels = c("K = 4", "K = 8 (Horizontal)", "K = 8 (Vertical)", "K = 16", "K = 64"))) %>% 
  ggplot()+
  geom_sf(fill = "transparent")+
  geom_sf_text(aes(label = grid_id), size = 1.5)+
  geom_sf(data = inocs_2024 %>% filter(treat == 1), shape = 23, aes(fill = block))+
  labs(x = "East (m)", y = "North (m)", title = "2024 Experimental Design (Single Inoculation)", fill = "Inoculation Point (Block)") +
  facet_grid(block~name)+
  theme_classic()+
  theme(legend.position = "bottom")

ggsave(filename = "experimental_design_2024.png", path = here("Reports/EEID_Presentation/figures/experimental/"), width = 8, height = 6, units = "in")

## 2025
grids_2025 <- rbind(design_2025$`6`$grid,
                    design_2025$`12`$grid,
                    design_2025$`72`$grid)

grids_2025 %>% 
  mutate(name = factor(name, levels = c("K = 6", "K = 12", "K = 72"))) %>% 
  ggplot()+
  geom_sf(fill = "transparent")+
  geom_sf_text(aes(label = grid_id), size = 1.5)+
  geom_sf(data = inocs_2025 %>% filter(treat == "single"), shape = 23, aes(fill = block))+
  labs(x = "East (m)", y = "North (m)", title = "2025 Experimental Design (Single Inoculation)", fill = "Inoculation Point (Block)") +
  facet_grid(block~name)+
  theme_classic()+
  theme(legend.position = "bottom")

ggsave(filename = "experimental_design_2025.png", path = here("Reports/EEID_Presentation/figures/experimental/"), width = 8, height = 6, units = "in")


# Spatial Predictions -----------------------------------------------------
## 2024
backward_long <- backward_2024 %>%
  filter(treat == 1) %>% 
  mutate(
    p_long = map(p_bar, ~{
      # Coerce column names to character to avoid type conflict
      colnames(.x) <- as.character(colnames(.x))
      rownames(.x) <- as.character(rownames(.x))
      
      as.data.frame(.x) %>%
        tibble::rownames_to_column("infected_component") %>%
        pivot_longer(
          cols = -infected_component,
          names_to = "source_component",
          values_to = "prob"
        )
    })
  ) %>%
  select(config, block, treat, visit, p_long) %>%
  unnest(p_long)

spatial_probs_2024 <- left_join(backward_long, grids_2024 %>% mutate(,
                                               grid_id = as.character(grid_id)), by = c("config","infected_component" = "grid_id"), relationship = "many-to-many")

nodes_2024 <- spatial_probs_2024 %>% filter(infected_component == source_component) %>% st_as_sf()
edges_2024 <- spatial_probs_2024 %>% filter(infected_component != source_component) %>% st_as_sf()

 
  ggplot() +
  geom_sf_label(data = nodes_2024 %>% filter(block == "A", config == "4"), aes(fill = prob, label = infected_component))+
  # geom_sf(data = edges %>% filter(plot_id == plt, configuration == config, is_max == 1, member_group != pred_group),
  #         aes(alpha = preds), 
  #         arrow = arrow(type = "closed", length = unit(0.075, "inches")))+
  geom_sf(data = grids_2024 %>% filter(config == "4"), fill = "transparent")+
  geom_sf(data = inocs %>% filter(plot_inoc_id == plt), size = 2, shape = 23, fill = "#F8766D") +
  facet_wrap(~visit, nrow = 1)+
  theme(legend.position = "bottom")+
  scale_alpha_continuous(limits = c(0, 1))+
  scale_fill_gradientn(colors = self_infection, limits = c(0,1)) + 
  labs(title = paste0("Spatial Network Predictions (", config_type, ", ", plt, ")"),
       fill = "Self Infection",
       alpha = "Dispersed Infection",
       x = "East",
       y = "North")+
  theme_classic() +
  theme(legend.position = "bottom")
