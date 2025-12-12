## ---------------------------
##
## Script name: 01b_DataModPrep.R
##
## Purpose of script: Sort data into clean structure for modeling 
##
## Author: Trent VanHawkins
##
## Date Created: 2025-07-17
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
stripe <- readRDS(here("DataProcessed/experimental/stripe_clean.rds"))
wind <- readRDS(here("DataProcessed/wind/wind_clean.rds"))
clusters <- readRDS(here("DataProcessed/experimental/clusters.rds"))
inocs <- readRDS(here("DataProcessed/experimental/inoc_sp.rds"))

# Create intensity array --------------------------------------------------
n_plants <- length(unique(stripe$plant_id))
n_blocks <- length(unique(stripe$block))
n_trt <- length(unique(stripe$treat))
n_visits <- length(unique(stripe$visit))
n_inocs <- length(unique(inocs$inoc_id))
intensity <- array(NA_real_, dim = c(n_plants, n_blocks, n_trt, n_visits),
                   dimnames = list(
                     plant = paste0(1:n_plants),  # no plant names
                     block = LETTERS[1:n_blocks],
                     treat = paste0(sort(unique(stripe$treat))),
                     visit = paste0(1:n_visits)))

for (blk in dimnames(intensity)$block) {
  for(trt in dimnames(intensity)$treat){
    for(vst in dimnames(intensity)$visit){
        intensity[,blk,trt,vst] <- stripe %>% 
          filter(block == blk, treat == as.numeric(trt), visit == as.numeric(vst)) %>% 
          arrange(plant_id) |> pull(intensity)
    }
  }
}

intensity[,'A', '1', '2']

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
                      treat = paste0(sort(unique(stripe$treat))),
                      visit = paste0(2:n_visits)))


## Add wind matrices
for (blk in dimnames(wind_array)$block) {
  for(trt in dimnames(wind_array)$treat){
    for(vst in dimnames(wind_array)$visit){
      first <- survey_periods %>% filter(block == blk, treat == as.numeric(trt), visit == as.numeric(vst)) %>% pull(date_prev)
      last <- survey_periods %>% filter(block == blk, treat == as.numeric(trt), visit == as.numeric(vst)) %>% pull(date)
      wind_array[,,blk,trt,vst] <- get_wind_mat(first_day = first, last_day = last, wind = wind, dir.mat = dir_mat)
    }
  }
}


# Get Inits ---------------------------------------------------------------
#Define Parameters
param_names <- c("beta", "delta", "gamma", "kappa", "phi")
n_params <- length(param_names)

#Define values of kappa to try
kappa_try <- seq(0.25,2.5,0.25)

# Assign Groups (For Backward Model) --------------------------------------
stripe_sp <- stripe %>% select(plant_id, east, north) %>% distinct() %>% arrange(plant_id) %>% st_as_sf(coords = c("east", "north"))

##Create grids (see DataFromat_Funs)
stripe_4 <- get_grid(stripe_sp, 2, 2, "K = 4")
stripe_8h <- get_grid(stripe_sp, 4, 2, "K = 8 (Horizontal)")
stripe_8v <- get_grid(stripe_sp, 2, 4, "K = 8 (Vertical)")
stripe_16 <- get_grid(stripe_sp, 4, 4, "K = 16")
stripe_64 <- get_grid(stripe_sp, 8, 8, "K = 64")

## Merge to list
grids <- list(
  "4" = stripe_4,
  "8h" = stripe_8h,
  "8v" = stripe_8v,
  "16" = stripe_16,
  "64" = stripe_64
)

saveRDS(grids, here("DataProcessed/experimental/grids_sp.rds"))

grid_dist <- map(grids, .f = ~st_distance(st_centroid(.x[["grid"]])))

grids_test <- map(grids, .f = ~.x[["grid"]] %>% mutate(centroid = st_centroid(geometry)))

plant_group <- array(NA_real_, dim = c(n_plants, length(grids)),
                   dimnames = list(
                     plant = paste0(1:n_plants),  # no plant names
                     config = names(grids)))

for (conf in dimnames(plant_group)$config) {
  plant_group[,conf] <- grids[[conf]][["points"]] |> arrange(plant_id) |> pull(grid_id)
}

inocs_clean <- inocs %>% select(block, treat, inoc_id, geometry)

true_infect <- array(NA_real_, dim = c(n_blocks, n_trt, n_inocs, length(grids)),
                     dimnames = list(block = LETTERS[1:n_blocks],
                                     treat = paste0(sort(unique(stripe$treat))),
                                     inoc = paste0(sort(unique(inocs$inoc_id))),
                                     config = names(grids)))


for (blk in dimnames(true_infect)$block) {
  for(trt in dimnames(true_infect)$treat){
    for(inoc in dimnames(true_infect)$inoc){
      for (conf in dimnames(true_infect)$config){
        print(paste0("Running Block: ", blk, " Treat: ", trt, " Inoc: ", inoc, "Config: ", conf))
        tmp <- inocs_clean %>% filter(block == blk, treat == trt, inoc_id == inoc) %>% pull(geometry)
        if(length(tmp) == 0){next}
        true_infect[blk,trt,inoc,conf] <- which.min(st_distance(x = tmp, y = st_centroid(grids[[conf]]$grid)))
        
    }
    }
    }
}


# Create Final List -------------------------------------------------------

mod_dat <- list(intensity = intensity,
                dist = dist_mat,
                wind = wind_array,
                groups = plant_group,
                truth = true_infect,
                grid_dist = grid_dist)


saveRDS(mod_dat, here("DataProcessed/experimental/mod_dat_arrays.rds"))
