## ---------------------------
##
## Script name: 05_BackwardModel.R
##
## Purpose of script: Implement backward model in R
##
## Author: Trent VanHawkins
##
## Date Created: 2025-06-03
##
##
## ---------------------------

## load up the packages we will need:  (uncomment as required)

library(tidyverse)
library(here)
library(sf)

stripe <- readRDS(here("DataProcessed/experimental/stripe_clean.rds"))

stripe_sp <- stripe %>% 
  st_as_sf(coords = c("east", "north"))

assign_grid <- function(dat, n_rows, n_cols){

  # Step 3: Create a grid (polygon) over the bounding box
  grid <- st_make_grid(dat, n = c(n_cols, n_rows), what = "polygons") %>%
    st_sf(grid_id = seq_len(length(.)))
  
  # Step 4: Spatial join to assign each plant to a grid cell
  plants_partitioned <- st_join(stripe_sp, grid)
  
  # Step 5: Assign centroid of grid
  grid_centroid <- st_centroid(grid)
  
  print(ggplot()+
    geom_sf(data = grid, aes(color = as.factor(grid_id)), fill = "transparent")+
    geom_sf(data = grid_centroid) +
    geom_sf(data = plants_partitioned, aes(color = as.factor(grid_id)))+
    theme_classic()+
    theme(legend.position = "none"))
  
  return(list("stripe_assigned" = plants_partitioned, 
              "grid" = grid,
              "centroid" = grid_centroid))
}

#Create boxes 
stripe_4 <- assign_grid(stripe_sp, n_rows = 2, n_cols = 2)
stripe_8h <- assign_grid(stripe_sp, n_rows = 4, n_cols = 2)
stripe_8v <- assign_grid(stripe_sp, n_rows = 2, n_cols = 4)
stripe_16 <- assign_grid(stripe_sp, n_rows = 4, n_cols = 4)

#Save them in a list
clusters <- list(
  stripe_4 = stripe_4,
  stripe_8h = stripe_8h,
  stripe_8v = stripe_8v,
  stripe_16 = stripe_16
)

saveRDS(clusters, here("DataProcessed/experimental/clusters.rds"))