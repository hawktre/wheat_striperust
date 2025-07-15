## ---------------------------
##
## Script name: 05e_EMsummary.R
##
## Purpose of script: Collect and format EM results
##
## Author: Trent VanHawkins
##
## Date Created: 2025-06-16
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

em_data <- readRDS(here("DataProcessed/experimental/em_dat.rds"))
sensitivity_fits <- readRDS(here("DataProcessed/results/backward_model/backwards_sensitivity_fits.rds"))
clusters <- readRDS(here("DataProcessed/experimental/clusters.rds"))
inocs <- readRDS(here("DataProcessed/experimental/inoc_sp.rds"))

## Create a large list object to store results
sensitivity_results <- list()
# Model selection ---------------------------------------------------------
## Sensitivity Analysis 
sensitivity_fit_summary <- imap_dfr(sensitivity_fits, function(plots, config_name) {
  imap_dfr(plots, function(visits, plot_id) {
    imap_dfr(visits, function(inits, visit) {
      imap_dfr(inits, function(visit_result, init_id) {
        tibble(
          configuration = config_name,
          plot_id = plot_id,
          visit = parse_number(visit),
          inits = init_id,
          em_iters = visit_result$em_iters,
          converged = visit_result$converged,
          final_neg_loglik = visit_result$final_neg_loglik, 
          !!!visit_result$final_theta
        )
      })
    })
  })
})

## Sort into dataframe to find change from MLE inits
params_sensitivity_long <- sensitivity_fit_summary %>% 
  filter(final_neg_loglik > -1e5) %>% 
  group_by(configuration, plot_id, visit) %>% 
  slice_min(final_neg_loglik) %>% 
  ungroup() %>% 
  select(configuration, plot_id, visit, inits, em_iters, final_neg_loglik, beta, delta, gamma, kappa, phi) %>% 
  pivot_longer(beta:phi, names_to = "param", values_to = "value") 

sensitivity_results[["params"]] <- params_sensitivity_long

best_mods <- paste0(params_sensitivity_long$configuration, "_", 
                    params_sensitivity_long$plot_id, "_visit", 
                    params_sensitivity_long$visit, "_", 
                    params_sensitivity_long$inits) %>% unique()

# Posterior predictions ---------------------------------------------------
## Sensitivity 
posterior_sensitivity <- imap_dfr(sensitivity_fits, function(plots, config_name) {
  imap_dfr(plots, function(visits, plot_id) {
    imap_dfr(visits, function(inits, visit) {
      imap_dfr(inits, function(init_info, init_name) {
        
        # Get p_mat
        p_mat <- init_info$p_mat
        if (is.null(p_mat)) return(NULL)
        
        # Get group_id from corresponding em_data object
        group_id <- em_data[[config_name]][[plot_id]][[visit]]$group_id
        if (is.null(group_id)) return(NULL)
        
        # Make data frame
        df <- as.data.frame(p_mat)
        colnames(df) <- paste0("p_group", seq_len(ncol(df)))
        df$group_id <- group_id
        
        # Summarize posterior probabilities by group
        df %>%
          group_by(group_id) %>%
          summarise(across(starts_with("p_group"), mean, .names = "mean_{.col}"),
                    .groups = "drop") %>%
          mutate(configuration = config_name,
                 plot_id = plot_id,
                 visit = visit,
                 init = init_name)
      })
    })
  })
})

# Arrange for easy viewing
posterior_summaries_long <- posterior_sensitivity %>%
  select(configuration, plot_id, visit, everything()) %>% 
  pivot_longer(cols = contains("mean"), names_to = "pred_group", values_to = "preds") %>% 
  drop_na(preds) %>% 
  rename("member_group" = "group_id") %>% 
  mutate(pred_group = parse_number(pred_group), 
         unique_id =  paste0(configuration, "_", 
                             plot_id, "_", 
                             visit, "_", 
                             init)) %>% 
  filter(unique_id %in% best_mods)


# Organize spatial information --------------------------------------------
grid <- bind_rows(clusters$stripe_4$grid,
                  clusters$stripe_8h$grid,
                  clusters$stripe_8v$grid,
                  clusters$stripe_16$grid) %>% 
  mutate(configuration = case_when(grid_type == "2 x 2" ~ "stripe_4",
                                   grid_type == "4 x 2" ~ "stripe_8h",
                                   grid_type == "2 x 4" ~ "stripe_8v",
                                   grid_type == "4 x 4" ~ "stripe_16")) %>% 
  rename("geometry" = ".")


centroid <- bind_rows(clusters$stripe_4$centroid,
                      clusters$stripe_8h$centroid,
                      clusters$stripe_8v$centroid,
                      clusters$stripe_16$centroid) %>% 
  mutate(configuration = case_when(grid_type == "2 x 2" ~ "stripe_4",
                                   grid_type == "4 x 2" ~ "stripe_8h",
                                   grid_type == "2 x 4" ~ "stripe_8v",
                                   grid_type == "4 x 4" ~ "stripe_16")) %>% 
  rename("geometry" = ".")

# Spatial Predictions -----------------------------------------------------

# Join centroids for spatial coordinates of each group
nodes <- posterior_summaries_long %>% 
  filter(member_group == pred_group) %>%
  left_join(centroid, by = c("member_group" = "grid_id", "configuration")) %>% 
  select(plot_id, visit, configuration, member_group, preds, geometry) %>% 
  st_as_sf()


# cross-infection probabilities 
edges <- posterior_summaries_long %>% 
  left_join(centroid, by = c("pred_group" = "grid_id", "configuration")) %>% 
  rename(from_geometry = "geometry") %>% 
  left_join(centroid, by = c("member_group" = "grid_id", "configuration")) %>% 
  rename(to_geometry = "geometry") %>% 
  select(plot_id, visit, configuration, preds, member_group, pred_group, to_geometry, from_geometry) %>% 
  rowwise() %>%
  mutate(geometry = st_sfc(st_linestring(rbind(st_coordinates(from_geometry),
                                               st_coordinates(to_geometry))),
                           crs = st_crs(from_geometry))) %>%
  ungroup() %>% 
  select(-c(to_geometry, from_geometry)) %>% 
  st_as_sf() %>% 
  group_by(configuration, plot_id, visit, member_group) %>% 
  mutate(is_max = if_else(preds == max(preds), 1, 0)) %>% 
  slice_max(preds) %>% 
  ungroup()

# Append to results
sensitivity_results[["preds"]] <- list("nodes" = nodes, "edges" = edges)
# Make Predictions----------------------------------------
single_inocs <- inocs %>%
  filter(grepl("1$", plot_inoc_id))

grid_truth <- grid %>% select(configuration, grid_id) 

for (plt in single_inocs$plot_inoc_id) {
  tmp <- single_inocs %>% filter(plot_inoc_id == plt)
  
  grid_truth[plt] <- as.numeric(lengths(st_contains(grid_truth, tmp)) > 0)
}

ground_truth <- grid_truth %>%
  pivot_longer(cols = A1:D1, names_to = "plot_id") %>% 
  filter(value == 1) %>% 
  select(configuration, plot_id, grid_id, value, geometry) %>% 
  arrange(plot_id)



prediction_df <- posterior_summaries_long %>% 
  filter(str_detect(plot_id, "1")) %>% 
  group_by(plot_id, configuration, visit) %>% 
  filter(preds == max(preds)) %>% 
  ungroup() %>% 
  select(configuration, plot_id, visit, pred_group) %>% 
  left_join(ground_truth, by = c("configuration", "plot_id")) %>% 
  rename("true_group" = grid_id) %>% 
  st_drop_geometry() %>% 
  select(-c(value, geometry))


# Compute error metric for each case --------------------------------------
library(tidymodels)

all_preds <- list()

for (config in unique(prediction_df$configuration)){
  dist_tmp <- centroid %>% filter(configuration == config) %>% st_distance()
  
  preds_tmp <- prediction_df %>% 
    filter(configuration == config) %>% 
    mutate(across(c(true_group, pred_group), ~ factor(.x, levels = 1:parse_number(config))))%>%
    mutate(dist_acc = map2_dbl(as.integer(pred_group), as.integer(true_group),
                                   ~ 1 - (dist_tmp[.x, .y]/max(dist_tmp))))
  
  all_preds[[config]] <- preds_tmp
}

# Append to results list
sensitivity_results[['singlesource_acc']] <- all_preds 

saveRDS(sensitivity_results,here("DataProcessed/results/backward_model/sensitivity_results.rds"))



