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
get_distance <- function(cur.plot, dat){
  pts <- dat %>% 
    st_as_sf(coords = c("east", "north")) %>% 
    filter(plotID_new == cur.plot) %>%
    select(plot, plant_id) %>% 
    distinct() 
  
  dist.mat <- st_distance(pts)
  
  dir.mat <- outer(1:nrow(pts), 1:nrow(pts), Vectorize(function(i, j) {
    dx <- st_coordinates(pts)[j, 1] - st_coordinates(pts)[i, 1]
    dy <- st_coordinates(pts)[j, 2] - st_coordinates(pts)[i, 2]
    atan2(dy, dx)
  }))
  
  return(list(dist = dist.mat,
              dir = dir.mat))
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



## For each period, describe get average wind speed in each cardinal direction
wind_summary <- function(first, last){
  wind %>% 
    filter(datetime >= first & datetime < last) %>% 
    group_by(cardinal) %>% 
    summarise(windrun = mean(speed)) %>% 
    ungroup() %>% 
    pivot_wider(names_from = cardinal, values_from = windrun)
  
}

survey_wind <- bind_cols(survey_periods, 
             bind_rows(apply(survey_periods, 1, function(x) wind_summary(first = x[['date_prev']],
                                                                         last = x[['date']]))))

# Join everything together ------------------------------------------------
## long-format windrun 
stripe.wind <- stripe.lag %>% 
  filter(visit != "visit1") %>% 
  left_join(survey_wind, by = c("plotID_new","visit")) %>% 
  pivot_longer(cols = N:NNW, names_to = "cardinal", values_to = "windrun")

