## ---------------------------
##
## Script name: 00b_DataCleaning_Wind.R
##
## Purpose of script: Clean and organize wind data
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
library(readxl)
library(circular)
library(climaemet)


# Format raw wind data ----------------------------------------------------
#Read in the data
stripe <- readRDS(here("DataProcessed/experimental/stripe_clean.rds"))
wind <- read_xlsx(here("DataRaw/local_wind.xlsx"), 
                  sheet = "Config 1",
                  skip = 2)

#Extract what we need 
wind.clean <- wind %>% 
  select(1, 6, 7) %>% 
  rename("datetime" = 1,
         "direction" = 2,
         "speed" = 3) %>% 
  filter(datetime >= ymd_hms('2024-04-09 00:00:00') & 
           datetime < ymd_hms('2024-06-16 00:00:00')) %>% 
  mutate(period = if_else(datetime < as_datetime(min(stripe$date)), "Before Surveys", "During Surveys"),
         month = month(datetime, label = T, abbr = F),
         direction.circ = circular(direction, units = "degrees", template = "geographics"))


# Compute Wind cardinal_dir ----------------------------------------------------

wind_run <- function(direction, n_directions){
  dir_labs <- c("N", "NNE", "NE", "ENE", "E", "ESE", 
                "SE", "SSE", "S", "SSW", "SW", "WSW", 
                "W", "WNW", "NW", "NNW")
  
  n_directions <- n_directions
  dir_bin_width <- 360 / n_directions
  dir_bin_cuts <- seq(dir_bin_width / 2, 360 - dir_bin_width / 2, dir_bin_width)
  
  dir_intervals <- findInterval(c(direction, dir_bin_cuts), dir_bin_cuts)
  dir_intervals[dir_intervals == n_directions] <- 0
  
  cardinal_dir <- head(factor(dir_intervals, labels = dir_labs), -n_directions)
  
  return(cardinal_dir)
}

wind.clean$cardinal <- wind_run(wind.clean$direction, n_directions = 16)

#Save the raw wind data
saveRDS(wind.clean, here("DataProcessed/wind/wind_clean.rds"))
