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

# Construct Wind Matrix ----------------------------------------------------------
get_wind_mat <- function(first_day, last_day, wind, dir.mat){
  #Subset the wind data to be in the appropriate frame
  wind.tmp <- wind %>% 
    filter(datetime >= first_day, datetime < last_day) %>% 
    group_by(cardinal, cardinal.dir) %>% 
    summarise(speed = mean(speed), .groups = "drop") 
  
  wind_angles <- wind.tmp[['cardinal.dir']]
  wind_speeds <- wind.tmp[['speed']]
  
  #initialize an empty matrix
  wind_projection_matrix <- matrix(0, nrow = nrow(dir.mat), ncol = ncol(dir.mat))
  
  for (i in 1:nrow(dir.mat)) {
    for (j in 1:ncol(dir.mat)) {
      if (i == j) next  # skip self-pairs if needed
      
      angle_ij <- dir.mat[i, j]
      
      # Difference between wind direction and direction from i to j
      angle_diff <- abs(atan2(sin(wind_angles - angle_ij), cos(wind_angles - angle_ij)))
      
      # Select wind vectors within π/2 of the direction from j to i
      in_cone <- angle_diff < (pi / 2)
      
      if (!any(in_cone)) next  # skip if no matching wind bins
      
      # Step 3: orthogonal projection: speed × cos(angle difference)
      projections <- wind_speeds[in_cone] * cos(angle_diff[in_cone])
      
      # Step 4: take average projected wind speed
      wind_projection_matrix[i, j] <- mean(projections)
    }
  }
  
  return(wind_projection_matrix)
}


# Create spatial Grid (For backward Model) --------------------------------
get_grid <- function(pts_sf, nrow_pts, ncol_pts, name) {
  
  stopifnot(inherits(pts_sf, "sf") && sf::st_geometry_type(pts_sf)[1] == "POINT")
  coords <- st_coordinates(pts_sf)
  
  # infer spacing (median is robust to tiny jitter)
  dx <- median(diff(sort(unique(coords[,1]))))
  dy <- median(diff(sort(unique(coords[,2]))))
  if (is.na(dx) || dx <= 0) stop("can't infer dx. Are x coordinates identical or only one column?")
  if (is.na(dy) || dy <= 0) stop("can't infer dy. Are y coordinates identical or only one row?")
  
  # define offset so that polygon centers == point centers
  offset_x <- min(coords[,1]) - dx/2
  offset_y <- min(coords[,2]) - dy/2

  # build grid covering bbox
  bbox_sfc <- st_as_sfc(st_bbox(pts_sf))
  grid <- st_make_grid(bbox_sfc,
                       cellsize = c(dx * nrow_pts, dy * ncol_pts),
                       offset = c(offset_x, offset_y),
                       what = "polygons") %>%
    st_sf(grid_id = seq_along(.), name = name, geometry = .)

  # assign points to grid cells
  hits <- st_intersects(pts_sf, grid)
  gid <- map_int(hits, ~ if(length(.x)) .x[1] else NA_integer_)
  pts_out <- pts_sf %>% mutate(grid_id = gid)

  list(points = pts_out, grid = grid)
}