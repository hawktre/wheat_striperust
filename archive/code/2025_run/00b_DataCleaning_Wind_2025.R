## ---------------------------
##
## Script name: 00b_DataCleaning_Wind_2025.R
##
## Purpose of script: Clean and organize wind data
##
## Author: Trent VanHawkins
##
## Date Created: 2025-11-01
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
stripe <- readRDS(here("DataProcessed/experimental/2025/stripe_clean_2025.rds"))
wind <- read_xlsx(here("DataRaw/local_wind_2025.xlsx"), 
                  sheet = "Config 1",
                  skip = 2)

#Extract what we need 
wind.clean <- wind %>% 
  select(1, 6, 7) %>% 
  rename("datetime" = 1,
         "direction" = 2,
         "speed" = 3) %>% 
  filter(datetime >= ymd(min(stripe$date)) & 
           datetime < ymd(max(stripe$date)))


# Compute Wind cardinal_dir ----------------------------------------------------

wind_cardinal <- function(direction, n_directions){
  
  dir_labs <- c("N", "NNE", "NE", "ENE", "E", "ESE", 
                "SE", "SSE", "S", "SSW", "SW", "WSW", 
                "W", "WNW", "NW", "NNW")
  
  n_directions <- n_directions
  dir_bin_width <- 360 / n_directions
  dir_bin_midpoint <- seq(0, 360 - dir_bin_width, dir_bin_width)*pi/180
  
  labs_key <- bind_cols(cardinal = dir_labs, rads = as.numeric(dir_bin_midpoint))
  dir_bin_cuts <- seq(dir_bin_width / 2, 360 - dir_bin_width / 2, dir_bin_width)
  
  dir_intervals <- findInterval(c(direction, dir_bin_cuts), dir_bin_cuts)
  dir_intervals[dir_intervals == n_directions] <- 0
  
  cardinal_dir <- head(factor(dir_intervals, labels = dir_labs), -n_directions)
  out <- data.frame("cardinal" = cardinal_dir) %>% 
    left_join(labs_key, by = "cardinal")
  
  return(out)
}

wind.clean$cardinal <- wind_cardinal(wind.clean$direction, n_directions = 16)[,1]
wind.clean$cardinal.dir <- wind_cardinal(wind.clean$direction, n_directions = 16)[,2]

#Save the raw wind data
saveRDS(wind.clean, here("DataProcessed/wind/wind_clean_2025.rds"))
