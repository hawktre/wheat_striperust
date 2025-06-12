## ---------------------------
##
## Script name: 05c_EMFormat.R
##
## Purpose of script: Format data for running backwards model
##
## Author: Trent VanHawkins
##
## Date Created: 2025-06-11
##
##
## ---------------------------

## view outputs in non-scientific notation

options(scipen = 6, digits = 4) 

## ---------------------------

## load up the packages we will need:  (uncomment as required)

library(tidyverse)
library(here)
library(sf)


# Read in the data and functions--------------------------------------------------------
clusters <- readRDS(here("DataProcessed/experimental/clusters.rds"))
forward <- readRDS(here("DataProcessed/results/all_fits_hurdle_logistic.rds"))
wind <- readRDS(here("DataProcessed/wind/wind_clean.rds"))


source(here("Code/01a_DataFormat_Fun.R"))


# Write loop to create data packets ---------------------------------------
em_dat <- list()

for (grid_type in names(clusters)) {
  cat("Creating data list object for ", grid_type, sep = "\n")
  # Create lagged data frame ------------------------------------------------
  stripe.lag <- clusters[[grid_type]][["stripe_assigned"]] %>%
    mutate(east = st_coordinates(.)[, 1],
           north = st_coordinates(.)[, 2]) %>%
    st_drop_geometry() %>%  
    arrange(plotID_new, plant_num, visit) %>% 
    group_by(plant_id) %>%
    mutate(intensity = intensity/100,
           intensity_prev = lag(intensity),
           date_prev = lag(date)) %>% 
    ungroup() %>% 
    filter(visit != "visit1")
  
  
  # Create data object for modeling -----------------------------------------
  ## Create a df of survey periods
  survey_periods <- stripe.lag %>% 
    select(plotID_new, visit, date, date_prev) %>% 
    distinct()
  
  mod_dat <- list()
  
  for (i in 1:nrow(survey_periods)){
    
    #subset the data to the current plot-visit combo
    cur.plot <- survey_periods$plotID_new[i]
    cur.visit <- survey_periods$visit[i]
    
    plt_visit <- paste0(cur.plot, "_", cur.visit, sep = "")
    
    plot_dat <- stripe.lag %>% 
      filter(plotID_new == cur.plot,
             visit == cur.visit)
    
    #Subset y_prev
    y_prev <-  plot_dat$intensity_prev
    
    #Subset y_cur
    y_cur <- plot_dat$intensity
    
    #y_cur > 0?
    z <- if_else(y_cur == 0, 0, 1)
    
    #Get Plant ID
    plant_id <- plot_dat$plant_id
    
    #Get Group_id
    group_id <- plot_dat$grid_id
    
    #distance and direction matrices
    dist_dir <- get_dist_dir(cur.plot, plot_dat)
    
    #wind matrix
    wind_mat <- get_wind_mat(first_day = survey_periods$date_prev[i],
                             last_day = survey_periods$date[i],
                             dir.mat = dist_dir$dir,
                             wind = wind)
    
    #Forward model estimates
    inits <- forward[["free"]][[plt_visit]][["theta"]][[1]]
    
    mod_dat[[cur.plot]][[cur.visit]] <- list("y_prev" = y_prev,
                                             "y_cur" = y_cur,
                                             "group_id" = group_id,
                                             "diseased" = z,
                                             "plant_id" = plant_id,
                                             "dist_mat" = dist_dir$dist,
                                             "dir_mat" = dist_dir$dir,
                                             "wind_mat" = wind_mat,
                                             "inits" = inits)
  }
  
  em_dat[[grid_type]] <- mod_dat
}

saveRDS(em_dat, here("DataProcessed/experimental/em_dat.rds"))
