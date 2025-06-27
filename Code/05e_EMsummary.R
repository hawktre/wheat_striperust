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
em_summary <- readRDS(here("DataProcessed/results/backwards_all_fits.rds"))
clusters <- readRDS(here("DataProcessed/experimental/clusters.rds"))
inocs <- readRDS(here("DataProcessed/experimental/inoc_sp.rds"))


# Model selection ---------------------------------------------------------
## Extract Model fit metrics
fit_df <- imap_dfr(em_summary, function(plots, config_name) {
  imap_dfr(plots, function(visits, plot_id) {
    imap_dfr(visits, function(visit_result, visit) {
      tibble(
        configuration = config_name,
        plot_id = plot_id,
        visit = parse_number(visit),
        em_iters = visit_result$em_iters,
        converged = visit_result$converged,
        final_neg_loglik = visit_result$final_neg_loglik
      )
    })
  })
})

## Compute BIC
best_mod <- fit_df %>% 
  mutate(n_mix = parse_number(configuration),
         bic_k = (n_mix - 1 + 5),
         bic = bic_k * log(64) - 2 * log(abs(final_neg_loglik))) %>%
  group_by(plot_id, visit) %>%
  slice_min(order_by = bic, n = 1, with_ties = FALSE) %>%
  ungroup()



posterior_summaries <- imap_dfr(em_summary, function(plots, config_name) {
  imap_dfr(plots, function(visits, plot_id) {
    imap_dfr(visits, function(visit_info, visit_name) {
      # Get p_mat
      p_mat <- visit_info$p_mat
      if (is.null(p_mat)) return(NULL)
      
      # Get group_id from corresponding em_data object
      group_id <- em_data[[config_name]][[plot_id]][[visit_name]]$group_id
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
               visit = visit_name)
    })
  })
})

# Arrange for easy viewing (optional)
posterior_summaries_long <- posterior_summaries %>%
  select(configuration, plot_id, visit, everything()) %>% 
  pivot_longer(cols = contains("mean"), names_to = "pred_group", values_to = "preds") %>% 
  drop_na(preds) %>% 
  rename("member_group" = "group_id") %>% 
  mutate(pred_group = parse_number(pred_group))

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


posterior_summaries_self_sp <- left_join(posterior_summaries_long, 
                                         grid, 
                                         by = c("configuration" = "configuration", 
                                                "pred_group" = "grid_id")) %>% 
  filter(member_group == pred_group) %>%
  st_as_sf()

posterior_summaries_others_sp <- left_join(posterior_summaries_long, 
                                           centroid, 
                                           by = c("configuration" = "configuration", "pred_group" = "grid_id")) %>% 
  filter(member_group != pred_group) %>% 
  st_as_sf()

# library(RColorBrewer)
# 
# for (config in unique(grid$configuration)) {
#   for(plt in unique(posterior_summaries_self_sp$plot_id)){
#     ggplot() +
#       geom_sf(data = grid %>% filter(configuration == config), fill = "transparent") +
#       geom_sf(data = posterior_summaries_others_sp %>% filter(plot_id == plt, configuration == config),
#               aes(color = preds, size = preds)) +
#       geom_sf_label(data = posterior_summaries_self_sp %>% filter(plot_id == plt, configuration == config),
#                     aes(label = member_group, fill = preds)) +
#       geom_sf(data = inocs %>% filter(plot_inoc_id == plt), size = 2, shape = 23, fill = "#F8766D") +
#       facet_grid(member_group ~ visit) +
#       scale_fill_gradientn(colours = brewer.pal(9, "YlOrRd")[3:9], limits = c(0,1)) +  # green to blue
#       scale_color_gradientn(colours = brewer.pal(9, "YlOrRd")[3:9], limits = c(0,1)) +  # orange to red
#       scale_size_continuous(limits = c(0,1))+
#       labs(fill = "Infection Prob.",
#            color = "Infection Prob.",
#            size = "Infection Prob.",
#            x = "East",
#            y = "North",
#            title = paste0("Infection Probabilities (", config, ", ", plt, ")")) +
#       theme(legend.position = "bottom") +
#       theme_bw()
#     
#     ggsave(filename = paste0(config, "_", plt, "_preds.png"), path = here("DataProcessed/results/figures/spatial_preds/"), width = 8, height = 8)
#   }
# }

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

# Make ground truth for prediction ----------------------------------------
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
  

pred_4   <- prediction_df %>% filter(configuration == "stripe_4") %>% mutate(across(c(true_group, pred_group), ~ factor(.x, levels = 1:4)))
pred_8h   <- prediction_df %>% filter(configuration == "stripe_8h")%>% mutate(across(c(true_group, pred_group), ~ factor(.x, levels = 1:8)))
pred_8v  <- prediction_df %>% filter(configuration == "stripe_8v")%>% mutate(across(c(true_group, pred_group), ~ factor(.x, levels = 1:8)))
pred_16  <- prediction_df %>% filter(configuration == "stripe_16")%>% mutate(across(c(true_group, pred_group), ~ factor(.x, levels = 1:16)))

library(tidymodels)
multi_metrics <- metric_set(accuracy, kap, f_meas)

metrics_config <- list(
  stripe_4  = pred_4,
  stripe_8h = pred_8h,
  stripe_8v = pred_8v,
  stripe_16 = pred_16
) |> 
  purrr::imap_dfr(~ multi_metrics(.x, truth = true_group, estimate = pred_group) |> 
                    mutate(configuration = .y), .id = NULL)

metrics_config_plot <- list(
  stripe_4  = pred_4,
  stripe_8h = pred_8h,
  stripe_8v = pred_8v,
  stripe_16 = pred_16
) |> 
  purrr::imap_dfr(~ .x %>%
                    group_by(plot_id) %>%
                    multi_metrics(truth = true_group, estimate = pred_group) %>%
                    mutate(configuration = .y),
                  .id = NULL)

metrics_config_visit <- list(
  stripe_4  = pred_4,
  stripe_8h = pred_8h,
  stripe_8v = pred_8v,
  stripe_16 = pred_16
) |> 
  purrr::imap_dfr(~ .x %>%
                    group_by(visit) %>%
                    multi_metrics(truth = true_group, estimate = pred_group) %>%
                    mutate(configuration = .y),
                  .id = NULL)


library(wesanderson)

self_infection <- wes_palette("Zissou1", 7, type = "continuous")

config_match <- data.frame(name = c("stripe_4", "stripe_8h", "stripe_8v", "stripe_16"),
                           config = c("2 x 2", "2 x 4", "4 x 2", "4 x 4"))


for (config in unique(edges$configuration)) {
  for (plt in unique(edges$plot_id)) {
    config_type <- config_match %>% filter(name == !!config) %>% .$config
    
     ggplot() +
      geom_sf_label(data = nodes %>% filter(plot_id == plt, configuration == config), aes(fill = preds, label = member_group))+
      geom_sf(data = edges %>% filter(plot_id == plt, configuration == config, is_max == 1, member_group != pred_group),
              aes(alpha = preds), 
              arrow = arrow(type = "closed", length = unit(0.075, "inches")))+
      geom_sf(data = grid %>% filter(configuration == config), fill = "transparent")+
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
    
    ggsave(filename = paste0(config, "_", plt, "_preds.png"), path = here("DataProcessed/results/figures/spatial_network_preds/"), width = 8, height = 8)
  }
}



