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
 
## Survey Data
stripe <- list()

blocks <- paste0(rep(LETTERS[1:4], each = 3), 1:3)

# Read in the data --------------------------------------------------------
## Experimental Data
for (sheet in blocks) {
  stripe[[sheet]] <- read_xlsx(here("DataRaw/StripeRust.xlsx"), sheet = sheet)
}

stripe.full <- map(stripe, ~ .x %>%
                     pivot_longer(cols = 4:8, 
                                  names_to = "date", 
                                  values_to = "intensity") %>%
                     mutate(date = str_extract(date, "^[^(]+") %>% mdy(),
                            visit = paste0("visit", rep(1:5, length.out = n())))) %>% 
  bind_rows()

## inoculum Locations
inoculum <- read_xlsx(here("DataRaw/StripeRust.xlsx"), sheet = "inoculum")

## Clean up names
names(stripe.full) <- tolower(names(stripe.full))

## Create Unique Plant ID
id_key <- stripe.full %>% 
  select(plot, north, east) %>% 
  distinct() %>% 
  group_by(plot) %>% 
  reframe(north = north,
          east = east,
          plant_id = paste0(plot, "_plant",row_number())) %>%
  ungroup()

## Get the total number of innocluations for each plot
inoc_key <- inoculum %>% 
  group_by(plot) %>% 
  summarise(inoculum_total = max(inoculum_num)) %>% 
  ungroup()

## Split plots into blocks and tx number
stripe.inoc <- stripe.full %>% 
  left_join(inoc_key, by = "plot") %>% 
  left_join(id_key, by = c("plot", "north", "east")) %>% 
  separate_wider_position(plot, 
                          widths = c("block" = 1, "rep" = 1),
                          cols_remove = F)

saveRDS(stripe.inoc, here("DataProcessed/experimental/stripe_clean.rds"))
