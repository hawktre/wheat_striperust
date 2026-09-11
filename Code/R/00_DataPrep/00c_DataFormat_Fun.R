## ---------------------------
##
## Script name: 01a_FormatData.R
##
## Purpose of script: Format the data for analysis
##
## Author: Trent VanHawkins
##
## Date Created: 2025-04-07
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

# Get distance and direction summaries ------------------------------------
## function for distance and direction matrices
get_dir <- function(coords) {
  n <- nrow(coords)
  
  # Initialize empty matrix for directions
  dir.mat <- matrix(NA_real_, nrow = n, ncol = n)
  
  # Compute angle from source j to target i
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        dx <- coords[i, 1] - coords[j, 1]  # target x - source x
        dy <- coords[i, 2] - coords[j, 2]  # target y - source y
        dir.mat[i, j] <- atan2(dx, dy)  %% (2*pi) # angle in radians
      }
    }
  }
  
  # Return results as a list
  return(dir.mat)
}

# Summarize Wind Run ------------------------------------------------------
summarize_wind_run <- function(first_day, last_day, wind) {
  cardinal_levels <- c(
    "N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
    "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"
  )

  # Survey dates do not contain times. Comparing calendar dates makes each
  # period start-inclusive and end-exclusive without a Date/POSIXct mismatch.
  wind_interval <- wind %>%
    filter(
      as.Date(datetime) >= as.Date(first_day),
      as.Date(datetime) < as.Date(last_day),
      !is.na(cardinal),
      !is.na(speed)
    )

  if (nrow(wind_interval) == 0) {
    stop("No wind observations found from ", first_day, " to ", last_day, ".")
  }

  # Define all 16 transport directions explicitly so directions not observed
  # during an interval remain in the summary with zero wind run.
  direction_key <- tibble(
    cardinal = factor(cardinal_levels, levels = cardinal_levels),
    cardinal.dir = seq(0, 337.5, by = 22.5) * pi / 180
  )

  observed_summary <- wind_interval %>%
    mutate(cardinal = factor(as.character(cardinal), levels = cardinal_levels)) %>%
    group_by(cardinal) %>%
    summarise(
      n_observations = n(),
      mean_speed = mean(speed),
      .groups = "drop"
    )

  direction_key %>%
    left_join(observed_summary, by = "cardinal") %>%
    mutate(
      n_observations = replace_na(n_observations, 0L),
      mean_speed = replace_na(mean_speed, 0),
      proportion_time = n_observations / nrow(wind_interval),
      wind_run = mean_speed * proportion_time
    )
}

# Construct Wind Matrix --------------------------------------------------
get_wind_mat <- function(first_day, last_day, wind, dir.mat){
  wind.tmp <- summarize_wind_run(first_day, last_day, wind)
  
  wind_angles <- wind.tmp[['cardinal.dir']]
  wind_runs <- wind.tmp[['wind_run']]
  
  #initialize an empty matrix
  wind_projection_matrix <- matrix(0, nrow = nrow(dir.mat), ncol = ncol(dir.mat))
  
  for (i in 1:nrow(dir.mat)) {
    for (j in 1:ncol(dir.mat)) {
      if (i == j) next  # skip self-pairs if needed
      
      angle_ij <- dir.mat[i, j]
      
      # Angular difference between transport direction and source j -> target i.
      angle_diff <- abs(atan2(sin(wind_angles - angle_ij), cos(wind_angles - angle_ij)))
      
      # Select wind vectors within π/2 of the direction from j to i
      in_cone <- angle_diff < (pi / 2)
      
      if (!any(in_cone)) next  # skip if no matching wind bins
      
      # Orthogonally project each eligible directional wind run.
      projections <- wind_runs[in_cone] * cos(angle_diff[in_cone])
      
      # Step 4: take average projected wind speed
      wind_projection_matrix[i, j] <- mean(projections)
    }
  }
  
  return(wind_projection_matrix)
}


# Create spatial Grid (For backward Model) --------------------------------
get_grid <- function(pts_sf, nrow, ncol, name) {

  # build grid covering bbox
  true_bbox <- st_bbox(pts_sf)
  dx <- true_bbox[["xmin"]]
  dy <- true_bbox[["ymin"]]
  bbox_sfc <- st_as_sfc(st_bbox(c(xmin = 0,
                                  ymin = 0,
                                  xmax = true_bbox[["xmax"]] + dx,
                                  ymax = true_bbox[["ymax"]] + dy)))
  grid <- st_make_grid(bbox_sfc,
                       n = c(ncol, nrow),
                       what = "polygons") %>%
    st_sf(grid_id = seq_along(.), name = name, geometry = .)

  # assign points to grid cells
  hits <- st_intersects(pts_sf, grid)
  gid <- map_int(hits, ~ .x[1] )
  pts_out <- pts_sf %>% mutate(grid_id = gid)

  list(points = pts_out, grid = grid)
}
