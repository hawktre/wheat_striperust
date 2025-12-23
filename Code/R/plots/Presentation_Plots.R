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
library(data.table)
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
backward_2024 <- readRDS(here("DataProcessed/results/backward_model/backward_fits_eval.rds"))
#backward_2025 <- readRDS(here("DataProcessed/results/backward_model/backward_fits_2025.rds"))

## Sensitivity Analysis
#backward_2024_sensitivity <- readRDS(here("DataProcessed/results/backward_model/backward_fits_sensitivity.rds"))

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
  geom_sf_text(aes(label = grid_id), size = 2)+
  geom_sf(data = inocs_2024 %>% filter(block == "B"), aes(shape = "Inoculation"), fill = "#CC0000", size = 1.5)+
  labs(x = "East (m)", y = "North (m)", fill = "Inoculation") +
  facet_grid(treat~name, labeller = labeller(
    treat = function(x) paste0(x, " inoculation(s)")
  ))+
  scale_shape_manual(values = 23, name = "Legend")+
  theme_classic()+
  theme(legend.position = "bottom")

ggsave(filename = "experimental_design_2024_blockB.png", path = here("Reports/EEID_Presentation/figures/experimental/"), width = 8, height = 6, units = "in")

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
backward_long_2024 <- backward_2024 %>%
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

spatial_probs_2024 <- left_join(backward_long_2024, grids_2024 %>% mutate(,
                                               grid_id = as.character(grid_id)), 
                                by = c("config","infected_component" = "grid_id"), relationship = "many-to-many") 

nodes_2024 <- spatial_probs_2024 %>% filter(infected_component == source_component) %>% st_as_sf()

centroids_2024 <- map(design_2024, .f = ~st_centroid(.x[["grid"]])) %>% 
  rbindlist() %>% 
  mutate(config = case_when(name == "K = 4" ~ "4",
                           name == "K = 8 (Horizontal)" ~ "8h",
                           name == "K = 8 (Vertical)" ~ "8v",
                           name == "K = 16" ~ "16",
                           name == "K = 64" ~ "64")) %>% 
  mutate(grid_id = as.character(grid_id)) %>% 
  st_as_sf()


edges_2024 <- left_join(spatial_probs_2024 %>% select(-geometry) %>% filter(infected_component != source_component), centroids_2024, by = c("config", "source_component" = "grid_id")) %>% 
  rename("from_geometry" = "geometry") %>% 
  left_join(centroids_2024, by = c("config", "infected_component" = "grid_id")) %>% 
  rename("to_geometry" = "geometry") %>% 
  select(config, block, treat, visit, name, infected_component, source_component, prob, from_geometry, to_geometry) %>% 
  rowwise() %>%
  mutate(geometry = st_sfc(st_linestring(rbind(st_coordinates(from_geometry),
                                               st_coordinates(to_geometry))),
                           crs = st_crs(from_geometry))) %>%
  ungroup() %>% 
  select(-c(to_geometry, from_geometry)) %>% 
  st_as_sf() %>% 
  group_by(config, block, treat, visit, infected_component) %>% 
  mutate(is_max = if_else(prob == max(prob), 1, 0)) %>% 
  slice_max(prob) %>% 
  ungroup()

spatplot <- function(cfg, blk, trt, nodes, edges, grid, inocs){
  nodes_plt <- nodes %>% filter(block == blk, config == cfg)
  edges_plt <- edges %>% filter(block == blk, config == cfg)
  grid_plt <- grid %>% filter(config == cfg)
  inocs_plt <- inocs %>% filter(block == blk, treat == trt)
  config_type <- unique(grid_plt$name)
  ggplot()+
    geom_sf(data = grid_plt, fill = "transparent")+
    geom_sf_label(data = nodes_plt, aes(fill = prob, label = infected_component))+
    geom_sf(data = edges_plt, aes(alpha = prob), arrow = arrow(type = "closed", length = unit(0.075, "inches")))+
    geom_sf(data = inocs_plt, aes(color = "Inoculation Point"), fill = "#F8766D", shape = 23, size = 3)+
    facet_wrap(~visit, nrow = 1)+
    scale_alpha_continuous(limits = c(0, 1))+
    scale_color_manual(values = "black")+
    scale_fill_gradientn(colors = intensity_cols, limits = c(0,1)) + 
    labs(title = paste0("Spatial Network Predictions (", config_type, ", Block ", blk, ")"),
         fill = "Self Infection",
         alpha = "Cross-Infection",
         x = "East",
         y = "North",
         color = "")+
    theme_classic() +
    theme(legend.position = "bottom")
  
}


for (letter in unique(spatial_probs_2024$block)) {
  for(k in unique(spatial_probs_2024$config)){
    plt <- spatplot(k, letter, 1, nodes_2024, edges_2024, grids_2024, inocs_2024)
    ggsave(plt, filename = paste0("spatplot_2024_", letter, k,".png"), path = here("Reports/EEID_Presentation/figures/backward/SpatialPlots/2024/"), width = 10, height = 4, units = "in")
  }
}

## 2025
backward_long_2025 <- backward_2025 %>%
  filter(treat == "single") %>% 
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

grids_2025 <- grids_2025 %>% mutate(config = case_when(name == "K = 6" ~ "6",
                                                       name == "K = 12" ~ "12",
                                                       name == "K = 72" ~ "72"))

spatial_probs_2025 <- left_join(backward_long_2025, grids_2025 %>% mutate(grid_id = as.character(grid_id)), 
                                by = c("config","infected_component" = "grid_id"), relationship = "many-to-many") 

nodes_2025 <- spatial_probs_2025 %>% filter(infected_component == source_component) %>% st_as_sf()

centroids_2025 <- map(design_2025, .f = ~st_centroid(.x[["grid"]])) %>% 
  rbindlist() %>% 
  mutate(config = case_when(name == "K = 6" ~ "6",
                            name == "K = 12" ~ "12",
                            name == "K = 72" ~ "72")) %>% 
  mutate(grid_id = as.character(grid_id)) %>% 
  st_as_sf()


edges_2025 <- left_join(spatial_probs_2025 %>% select(-geometry) %>% filter(infected_component != source_component), centroids_2025, by = c("config", "source_component" = "grid_id")) %>% 
  rename("from_geometry" = "geometry") %>% 
  left_join(centroids_2025, by = c("config", "infected_component" = "grid_id")) %>% 
  rename("to_geometry" = "geometry") %>% 
  select(config, block, treat, visit, name, infected_component, source_component, prob, from_geometry, to_geometry) %>% 
  rowwise() %>%
  mutate(geometry = st_sfc(st_linestring(rbind(st_coordinates(from_geometry),
                                               st_coordinates(to_geometry))),
                           crs = st_crs(from_geometry))) %>%
  ungroup() %>% 
  select(-c(to_geometry, from_geometry)) %>% 
  st_as_sf() %>% 
  group_by(config, block, treat, visit, infected_component) %>% 
  mutate(is_max = if_else(prob == max(prob), 1, 0)) %>% 
  slice_max(prob) %>% 
  ungroup()

for (letter in unique(spatial_probs_2025$block)) {
  for(k in unique(spatial_probs_2025$config)){
    plt <- spatplot(k, letter, "single", nodes_2025, edges_2025, grids_2025, inocs_2025)
    ggsave(plt, filename = paste0("spatplot_2025_", letter, k,".png"), path = here("Reports/EEID_Presentation/figures/backward/SpatialPlots/2025/"), width = 10, height = 4, units = "in")
  }
}


# General Structure Plots -------------------------------------------------

ggplot()+
  geom_sf(data = stripe_2024 %>% filter(block == "A"), aes(color = intensity))+
  geom_sf(data = inocs_2024 %>% filter(block == "A"), aes(fill = "Inoculation Point"), shape = 23, size = 3)+
  facet_grid(treat ~ date, labeller = labeller(
    treat = function(x) paste0("Inocs: ", x)
  ))+
  labs(x = "East (m)",
       y = "North (m)",
       fill = "", color = "Intensity")+
  scale_color_gradientn(colors = intensity_cols, limits = c(0,1))+
  theme_bw()

ggsave(filename = "rawdata_blockA_2024.png", path = here("Reports/EEID_Presentation/figures/experimental/"), width = 10, height = 5, units = "in")

ggplot()+
  geom_sf(data = stripe_2025 %>% filter(block == "B"), aes(color = intensity))+
  geom_sf(data = inocs_2025 %>% filter(block == "B"), aes(fill = "Inoculation Point"), shape = 23, size = 3)+
  facet_grid(treat ~ date)+
  labs(x = "East (m)",
       y = "North (m)",
       fill = "", color = "Intensity")+
  scale_color_gradientn(colors = intensity_cols, limits = c(0,1))+
  theme_bw()

ggsave(filename = "rawdata_blockA_2025.png", path = here("Reports/EEID_Presentation/figures/experimental/"), width = 10, height = 5, units = "in")


# Table for Backwards Model Results ---------------------------------------
backward_2024 %>% 
  mutate(accuracy = if_else(predicted_source == true_source, 1, 0),
         config = factor(config, levels = c("4", "8h", "8v", "16", "64"), labels = c("2 x 2", "4 x 2", "2 x 4", "4 x 4", "8 x 8"))) %>% 
  filter(treat == 1) %>% 
  group_by(config, visit) %>% 
  summarise(acc = mean(accuracy),
            dist_acc = mean(dist_acc),
            dist_error = mean(dist_error)) %>% 
  ungroup() %>% 
  kable(format = "latex", digits = 3, col.names = c("Configuration", "Visit", "Accuracy", "Accuracy (Distance-Weighted)", "Error (m)")) %>% 
  kableExtra::kable_styling()

backward_2025 %>% 
  mutate(accuracy = if_else(predicted_source == true_source, 1, 0),
         config = factor(config, levels = c("6", "12", "72"), labels = c("2 x 3", "4 x 3", "7 x 9"))) %>% 
  filter(treat == "single") %>% 
  group_by(config, visit) %>% 
  summarise(acc = mean(accuracy),
            dist_acc = mean(dist_acc),
            dist_error = mean(dist_error)) %>% 
  ungroup() %>% 
  kable(format = "latex", 
        digits = 3, 
        booktabs = T,
        col.names = c("Configuration", "Visit", "Accuracy", "Accuracy (Distance-Weighted)", "Error (m)")) %>% 
  kableExtra::kable_styling()

backward_2024_t1 <- backward_2024 %>% 
  mutate(accuracy = if_else(predicted_source == true_source, 1, 0),
         config = factor(config, levels = c("4", "8h", "8v", "16", "64"), 
                         labels = c("2 x 2", "4 x 2", "2 x 4", "4 x 4", "8 x 8")),
         init = "MLE") %>% 
  filter(treat == 1) %>% 
  select(config, block, treat, visit, em_iters, predicted_source, true_source, dist_acc, dist_error, accuracy, init)

backward_2024_sensitivity_min <- backward_2024_sensitivity %>% 
  group_by(config, block, treat, visit) %>% 
  slice_min(Q_final) %>% 
  ungroup() %>% 
  mutate(accuracy = if_else(predicted_source == true_source, 1, 0),
         config = factor(config, levels = c("4", "8h", "8v", "16", "64"), 
                         labels = c("2 x 2", "4 x 2", "2 x 4", "4 x 4", "8 x 8")),
         init = "Sensitivity") %>% 
  select(names(backward_2024_t1))

backward_2024_full <- rbind(backward_2024_t1, backward_2024_sensitivity_min)  

backward_2024_full %>% 
  group_by(init, visit) %>% 
  summarise(acc = mean(accuracy),
            dist_acc = mean(dist_acc),
            dist_error = mean(dist_error)) %>% 
  ungroup() %>% 
  kable() |> 
  kable_styling()
  pivot_wider(names_from = init, values_from = c(dist_acc, acc, dist_error)) %>% 
  select(visit, starts_with("acc"), starts_with("dist_acc"), starts_with("dist_error")) %>% 
  kable(
    #format = "latex",
    booktabs = TRUE,
    digits = 3,
    col.names = c(
      "Visit",
      "MLE", "Sensitivity",
      "MLE", "Sensitivity",
      "MLE", "Sensitivity"
    )
  ) %>%
  add_header_above(c(
    " " = 1,                     # empty for 'Visit'
    "Accuracy" = 2,
    "Distance-Weighted Accuracy" = 2,
    "Error (m)" = 2
  ))

discrete_cols <- wesanderson::wes_palette("Darjeeling1", n = 5, type = "discrete")

backward_2024_full %>% 
  ggplot(aes(x = em_iters))+
  geom_histogram(aes(fill = init), position = "identity", alpha = 0.7, color = "black")+
  labs(x = "EM Iterations", y = "Count", fill = "Initialization")+
  scale_fill_manual(values = discrete_cols)+
  theme_classic()


# Simulation Results ------------------------------------------------------
sims <- readRDS(here("DataProcessed/results/simulation/sims_appended.rds"))

sims_clean <- sims %>% 
  mutate(config = factor(config, levels = c("4", "8h", "8v", "16"),
                         labels = c("2 x 2", "4 x 2", "2 x 4", "4 x 4"))) %>% 
  filter(treat == 1)

sims_clean %>% 
  filter(block == "A", config %in% c("2 x 2", "4 x 4")) %>% 
  select(-block) %>% 
  group_by(config, visit) %>%
  summarise("acc" = mean(predicted_source == true_source),
            "dist_acc" = mean(dist_acc),
            "error" = mean(dist_error),
            "miss_error" = mean(dist_error[predicted_source != true_source])) %>% 
  kable(format = "latex",
        col.names = c("Configuration", "Visit", "Acc.", "Dist. Acc.", "Dist. Error", "Miss Error"),
        booktabs = T)


# Comparing shared vs individual component parameters --------------------
backward_2024 |> 
  mutate(config = factor(config, levels = c("4", "8h", "8v", "16", "64"))) |> 
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
  tab_spanner(columns = contains("error"), label = "Error (m)")

