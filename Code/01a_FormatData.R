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


# Read in the data --------------------------------------------------------
stripe <- readRDS(here("DataProcessed/experimental/stripe_clean.rds"))
wind <- readRDS(here("DataProcessed/wind/wind_clean.rds"))


# Create lagged data frame ------------------------------------------------
stripe.lag <- stripe %>% 
  arrange(plant_id, visit) %>% 
  group_by(plant_id) %>%
  mutate(intensity_prev = lag(intensity),
         date_prev = lag(date)) %>% 
  ungroup()

# Get distance and direction summaries ------------------------------------
## function for distance and direction matrices
get_distance <- function(cur.plot, dat) {
  # Filter to current plot and ensure unique plant locations
  pts <- dat %>%
    filter(plotID_new == cur.plot) %>%
    distinct(plot, plant_id, east, north) %>%
    st_as_sf(coords = c("east", "north"))
  
  # Compute pairwise distance matrix
  dist.mat <- st_distance(pts)
  
  # Extract (x, y) coordinates matrix
  coords <- st_coordinates(pts)
  n <- nrow(coords)
  
  # Initialize empty matrix for directions
  dir.mat <- matrix(NA, nrow = n, ncol = n)
  
  # Compute angle from source j to target i
  for (i in 1:n) {
    for (j in 1:n) {
      dx <- coords[j, 1] - coords[i, 1]  # target x - source x
      dy <- coords[j, 2] - coords[i, 2]  # target y - source y
      dir.mat[i, j] <- atan2(dy, dx)     # angle in radians
    }
  }
  
  # Return results as a list
  list(dist = dist.mat, dir = dir.mat)
}

## Compute distance matrix for each plot
plotlist <- sort(unique(stripe.lag$plotID_new))
dist.mats <- lapply(plotlist, function(x) get_distance(x, stripe.lag))
names(dist.mats) <- plotlist


# Survey Periods (for wind) ----------------------------------------------------------
## Create a df of survey periods
survey_periods <- stripe.lag %>% 
  select(plotID_new, visit, date, date_prev) %>% 
  distinct() %>% 
  filter(visit != "visit1") 


## For each period, describe get average wind speed in the range of the provided cardinal direction
wind_mat <- function(first_day, last_day, wind, dir.mat){
  browser()
  #Subset the wind data to be in the appropriate frame
  wind.tmp <- wind %>% 
    filter(datetime >= first_day, datetime < last_day) %>% 
    group_by(cardinal, cardinal.dir) %>% 
    summarise(speed = mean(speed)) %>% 
    ungroup()
  
  wind_angles <- wind.tmp[['cardinal.dir']]
  wind_speeds <- wind.tmp[['speed']]
  
  #initialize an empty matrix
  wind_projection_matrix <- matrix(NA, nrow = nrow(dir.mat), ncol = ncol(dir.mat))
  
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


# Step 1: add dir.mat and wind columns to survey_periods
survey_inputs <- survey_periods %>%
  mutate(
    dir.mat = map(plotID_new, ~ dist.mats[[.x]]$dir),
    wind = list(wind)  # same wind object repeated for all
  )

# Step 2: group by plot
survey_nested <- survey_inputs %>%
  group_by(plotID_new) %>%
  group_split()  # gives list of tibbles, one per plot

# Step 3: map over each plot group
wind.mats <- map(survey_nested, function(df) {
  visits <- df$visit  # for naming
  
  inner <- pmap(
    list(first_day = df$date_prev,
         last_day = df$date,
         wind = df$wind,
         dir.mat = df$dir.mat),
    wind_mat
  )
  
  names(inner) <- visits
  inner
})

# Step 4: name outer list by plot ID
names(wind.mats) <- map_chr(survey_nested, ~ .x$plotID_new[1])

saveRDS(wind.mats, here("DataProcessed/wind/wind_mats.rds"))