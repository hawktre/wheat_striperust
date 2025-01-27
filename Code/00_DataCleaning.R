## ---------------------------
##
## Script name: 00_DataCleaning.R
##
## Purpose of script: Read-in and familiarize data
##
## Author: Trent VanHawkins
##
## Date Created: 2024-11-21
##
##
## ---------------------------

## ---------------------------

## load up the packages we will need:  (uncomment as required)

library(tidyverse)
library(here)
library(readxl)



# Read in the data --------------------------------------------------------
## Survey Data
stripe <- list()

blocks <- paste0(rep(LETTERS[1:4], each = 3), 1:3)

for (sheet in blocks) {
  stripe[[sheet]] <- read_xlsx("DataRaw/StripeRust.xlsx", sheet = sheet)
}

stripe.full <- map(stripe, ~ .x %>%
      pivot_longer(cols = 4:8, 
                   names_to = "date", 
                   values_to = "intensity") %>%
      mutate(date = str_extract(date, "^[^(]+") %>% mdy(),
             visit = paste0("visit", rep(1:5, length.out = n())))) %>% 
  bind_rows()

## Lowercase names
names(stripe.full) <- tolower(names(stripe.full))

## Create Unique Plant ID
stripe.full <- stripe.full %>% 
  group_by(plot, north, east) %>% 
  mutate(plant_id = paste0(plot, "_plant", cur_group_id())) %>% 
  ungroup

## Innoculum Locations
innoculum <- read_xlsx(here("DataRaw/StripeRust.xlsx"), sheet = "innoculum")

## Weather Data
dailywind_raw <- read.csv(here("DataRaw/daily_weather.csv")) 
hourlywind_raw <- read.csv(here("DataRaw/hourly_weather.csv"))

names(dailywind_raw) <- c("date", "speed", "direction", "windrun")
names(hourlywind_raw) <- c("datetime", "direction", "speed")

##Clean up weather names 
wind <- hourlywind_raw %>% 
  mutate(datetime = mdy_hm(datetime)) %>% 
  as_tibble()


# Join stripe rust data with innoculum ------------------------------------
stripe.innoc <- full_join(stripe.full, innoculum, by = "plot", 
                          relationship = "many-to-many",
                          suffix = c(".survey", ".innoc"))

## Compute difference in cartesian from each innoculum
stripe.innoc <- stripe.innoc %>% 
  mutate(dx = east.survey - east.innoc,
         dy = north.survey - north.innoc,
         r = sqrt(dx**2 + dy**2),
         phi = atan2(dy, dx))

# ## Pivot wider again
# stripe.innoc <- stripe.innoc %>% pivot_wider(
#   id_cols = c(plant_id, visit, date, intensity),
#   names_from = innoculum_num,     
#   values_from = c(r, phi),        
#   names_sep = "_innoculum"        
# )


# Format wind data --------------------------------------------------------

cardinal_directions <- tibble(
  direction = c("N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
                "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"),
  lower_bound = seq(348.75, by = 22.5, length.out = 16) %% 360,
  upper_bound = seq(11.25, by = 22.5, length.out = 16) %% 360
)

wind <- wind %>% 
  mutate(cardinal = case_when(direction >= cardinal_directions$lower_bound[2] & direction < cardinal_directions$upper_bound[2] ~ cardinal_directions$direction[2],
                                         direction >= cardinal_directions$lower_bound[3] & direction < cardinal_directions$upper_bound[3] ~ cardinal_directions$direction[3],
                                         direction >= cardinal_directions$lower_bound[4] & direction < cardinal_directions$upper_bound[4] ~ cardinal_directions$direction[4],
                                         direction >= cardinal_directions$lower_bound[5] & direction < cardinal_directions$upper_bound[5] ~ cardinal_directions$direction[5],
                                         direction >= cardinal_directions$lower_bound[6] & direction < cardinal_directions$upper_bound[6] ~ cardinal_directions$direction[6],
                                         direction >= cardinal_directions$lower_bound[7] & direction < cardinal_directions$upper_bound[7] ~ cardinal_directions$direction[7],
                                         direction >= cardinal_directions$lower_bound[8] & direction < cardinal_directions$upper_bound[8] ~ cardinal_directions$direction[8],
                                         direction >= cardinal_directions$lower_bound[9] & direction < cardinal_directions$upper_bound[9] ~ cardinal_directions$direction[9],
                                         direction >= cardinal_directions$lower_bound[10] & direction < cardinal_directions$upper_bound[10] ~ cardinal_directions$direction[10],
                                         direction >= cardinal_directions$lower_bound[11] & direction < cardinal_directions$upper_bound[11] ~ cardinal_directions$direction[11],
                                         direction >= cardinal_directions$lower_bound[12] & direction < cardinal_directions$upper_bound[12] ~ cardinal_directions$direction[12],
                                         direction >= cardinal_directions$lower_bound[13] & direction < cardinal_directions$upper_bound[13] ~ cardinal_directions$direction[13],
                                         direction >= cardinal_directions$lower_bound[14] & direction < cardinal_directions$upper_bound[14] ~ cardinal_directions$direction[14],
                                         direction >= cardinal_directions$lower_bound[15] & direction < cardinal_directions$upper_bound[15] ~ cardinal_directions$direction[15],
                                         direction >= cardinal_directions$lower_bound[16] & direction < cardinal_directions$upper_bound[16] ~ cardinal_directions$direction[16],
                                         TRUE ~ "N"))
avg_speed <- mean(wind$speed)

wind %>% 
  group_by(cardinal) %>% 
  summarise(wind_run = avg_speed * (n()/nrow(wind)))
