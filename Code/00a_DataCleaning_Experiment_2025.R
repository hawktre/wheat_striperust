## ---------------------------
##
## Script name: 00a_DataCleaning_Experiment_2025.R
##
## Purpose of script: Clean 2025 Experimental data 
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
## load up the packages we will need:  (uncomment as required)
library(tidyverse)
library(here)
library(readxl)
library(sf)
 
## Survey Data
stripe <- list()

blocks <- paste0(rep(LETTERS[1:3], each = 4), 1:4)

# Read in the data --------------------------------------------------------
## Experimental Data
for (sheet in blocks) {
  stripe[[sheet]] <- read_xlsx(here("DataRaw/StripeRust_2025_norepeat.xlsx"), sheet = sheet)
}

stripe_full <- map(stripe, ~ .x |> 
                     pivot_longer(cols = 4:7, 
                                  names_to = "date", 
                                  values_to = "intensity") %>%
                     mutate(date = mdy(str_remove(date, "\\(\\%\\)")),
                            visit = as.factor(rep(1:4, length.out = n())),
                            intensity = intensity/100))|> 
  bind_rows() |> 
  select(Plot, North, East, date, visit, intensity)

## Clean up names
names(stripe_full) <- tolower(names(stripe_full))

## inoculum Locations
inoculum <- read_xlsx(here("DataRaw/StripeRust_2025_norepeat.xlsx"), sheet = "inoculum")

## Create Unique Plant ID
id_key <- stripe_full %>% 
  select(north, east) %>% 
  distinct() %>% 
  arrange(north, east) |> 
  mutate(plant_id = row_number())
