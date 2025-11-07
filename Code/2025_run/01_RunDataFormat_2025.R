## ---------------------------
##
## Script name: 01_RunDataFormat_2025.R
##
## Purpose of script: Sort data into clean structure for modeling 
##
## Author: Trent VanHawkins
##
## Date Created: 2025-11-01
##
##
## ---------------------------


options(scipen = 6, digits = 4) 

## ---------------------------

## load up the packages we will need:  (uncomment as required)

library(tidyverse)
library(here)
library(sf)
source(here("Code/01a_DataFormat_Fun.R"))
source(here("Code/02a_ForwardGradFun.R"))

# Read in the data --------------------------------------------------------
stripe <- readRDS(here("DataProcessed/experimental/2025/stripe_clean_2025.rds"))
wind <- readRDS(here("DataProcessed/wind/wind_clean_2025.rds"))
inocs <- readRDS(here("DataProcessed/experimental/2025/inocs_2025.rds"))
# Create intensity array --------------------------------------------------
n_plants <- length(unique(stripe$plant_id))
n_blocks <- length(unique(stripe$block))
n_trt <- length(unique(stripe$treat))
n_visits <- length(unique(stripe$visit))

intensity <- array(NA_real_, dim = c(n_plants, n_blocks, n_trt, n_visits),
                   dimnames = list(
                     plant = paste0(1:n_plants),  # no plant names
                     block = LETTERS[1:n_blocks],
                     treat = sort(unique(stripe$treat)),
                     visit = paste0(1:n_visits)))

for (blk in dimnames(intensity)$block) {
  for(trt in dimnames(intensity)$treat){
    for(vst in dimnames(intensity)$visit){
      intensity[,blk,trt,vst] <- stripe %>% 
        filter(block == blk, treat == trt, visit == as.numeric(vst)) %>% 
        arrange(plant_id) |> pull(intensity)
    }
  }
}

# Create Distance and Wind Matrices -------------------

## Create a df of survey periods
survey_periods <- stripe %>% 
  select(block, treat, visit, date) %>% 
  distinct() %>% 
  group_by(block, treat) %>% 
  mutate(date_prev = lag(date)) %>% 
  ungroup() %>% 
  filter(visit != 1)

## Distance Matrix
dist_mat <- stripe %>% 
  select(plant_id, north, east) %>% 
  distinct() %>% 
  arrange(plant_id) |> 
  st_as_sf(coords = c("east", "north")) %>% 
  st_distance()

## Directional Matrix (for wind)
coords <- stripe %>% 
  select(plant_id, east, north) %>% 
  distinct() %>% 
  arrange(plant_id) |> 
  select(-plant_id) %>% 
  as.matrix()

dir_mat <- get_dir(coords)

## Wind array
wind_array <- array(NA_real_, dim = c(n_plants, n_plants, n_blocks, n_trt, n_visits - 1),
                    dimnames = list(
                      NULL,
                      NULL,
                      block = LETTERS[1:n_blocks],
                      treat = sort(unique(stripe$treat)),
                      visit = paste0(2:n_visits)))


## Add wind matrices
for (blk in dimnames(wind_array)$block) {
  for(trt in dimnames(wind_array)$treat){
    for(vst in dimnames(wind_array)$visit){
      first <- survey_periods %>% filter(block == blk, treat == trt, visit == as.numeric(vst)) %>% pull(date_prev)
      last <- survey_periods %>% filter(block == blk, treat == trt, visit == as.numeric(vst)) %>% pull(date)
      wind_array[,,blk,trt,vst] <- get_wind_mat(first_day = first, last_day = last, wind = wind, dir.mat = dir_mat)
    }
  }
}


# Assign Groups (For Backward Model) --------------------------------------
stripe_sp <- stripe %>% select(plant_id, east, north) %>% distinct() %>% st_as_sf(coords = c("east", "north")) %>% arrange(plant_id)

##Create grids (see DataFromat_Funs)
stripe_6 <- get_grid(stripe_sp, 2, 3, "K = 6")
stripe_12 <- get_grid(stripe_sp, 4, 3, "K = 12")
stripe_72 <- get_grid(stripe_sp, 8, 9, "K = 72")

## Merge to list
grids <- list(
  "6" = stripe_6,
  "12" = stripe_12,
  "72" = stripe_72
)

saveRDS(grids, here("DataProcessed/experimental/2025/grids_sp_2025.rds"))

grid_dist <- map(grids, .f = ~st_distance(st_centroid(st_geometry(.x[["grid"]]))))

plant_group <- array(NA_real_, dim = c(n_plants, length(grids)),
                     dimnames = list(
                       plant = paste0(1:n_plants),  # no plant names
                       config = names(grids)))

for (conf in dimnames(plant_group)$config) {
  plant_group[,conf] <- grids[[conf]][["points"]] |> arrange(plant_id) |> pull(grid_id)
}

single_inocs <- inocs %>% st_as_sf(coords = c("east", "north")) %>% select(block, treat, inoc_id, geometry) %>% 
  filter(treat == "single")

true_infect <- array(NA_real_, dim = c(n_blocks, length(grids)),
                     dimnames = list(block = LETTERS[1:n_blocks],
                                     config = names(grids)))
for (blk in dimnames(true_infect)$block) {
  for (conf in dimnames(true_infect)$config) {
    true_infect[blk,conf] <- single_inocs %>% filter(block == blk) %>% 
      st_join(grids[[conf]][["grid"]]) %>% pull(grid_id)
  }
}

# Create Final List -------------------------------------------------------

mod_dat <- list(intensity = intensity,
                dist = dist_mat,
                wind = wind_array,
                groups = plant_group,
                truth = true_infect,
                grid_dist = grid_dist)


saveRDS(mod_dat, here("DataProcessed/experimental/2025/mod_dat_arrays_2025.rds"))
