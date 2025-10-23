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
library(sf)
 
## Survey Data
stripe <- list()

blocks <- paste0(rep(LETTERS[1:4], each = 3), 1:3)

# Read in the data --------------------------------------------------------
## Experimental Data
for (sheet in blocks) {
  stripe[[sheet]] <- read_xlsx(here("DataRaw/StripeRust.xlsx"), sheet = sheet)
}

stripe_full <- map(stripe, ~ .x |> 
                     pivot_longer(cols = 4:8, 
                                  names_to = "date", 
                                  values_to = "intensity") %>%
                     mutate(date = mdy(str_remove(date, "\\(\\%\\)")),
                            visit = as.factor(rep(1:5, length.out = n())),
                            intensity = intensity/100))|> 
  bind_rows() |> 
  select(Plot, North, East, date, visit, intensity)

## Clean up names
names(stripe_full) <- tolower(names(stripe_full))

## inoculum Locations
inoculum <- read_xlsx(here("DataRaw/StripeRust.xlsx"), sheet = "inoculum")

## Create Unique Plant ID
id_key <- stripe_full %>% 
  select(north, east) %>% 
  distinct() %>% 
  arrange(north, east) |> 
  mutate(plant_id = row_number())

## Get the total number of innocluations for each plot
inoc_key <- inoculum %>% 
  group_by(plot) %>% 
  summarise(n_inocs = max(inoculum_num)) %>% 
  ungroup() |> 
  mutate(block = substr(plot, 1, 1), 
plot_new = paste0(block, n_inocs))

inoc_sp <- inoculum |> 
  select(plot, inoculum_num, north, east) |> 
  rename("inoc_id" = inoculum_num) |> 
  left_join(inoc_key, by = "plot") |> 
  select(plot_new, inoc_id, north, east) |> 
  mutate(block = substr(plot_new, 1, 1),
treat = as.factor(as.numeric(substr(plot_new, 2, 2))),
inoc_id = as.factor(inoc_id)) |> 
  rename("plot" = plot_new) |> 
  select(plot, block, treat, inoc_id, north, east) |> 
  st_as_sf(coords = c("east", "north"))

inoc_sp |> 
  ggplot()+
  geom_sf()+
  facet_grid(block ~ treat)

saveRDS(inoc_sp, here("DataProcessed/experimental/inoc_sp.rds"))
## Split plots into blocks and tx number
stripe_clean <- stripe_full %>% 
  left_join(inoc_key, by = "plot") %>% 
  left_join(id_key, by = c("north", "east")) |> 
  select(plot_new, north, east, date, visit, plant_id, intensity) |> 
  rename("plot" = plot_new) |> 
  mutate(block = as.factor(substr(plot, 1, 1)),
treat = as.factor(as.numeric(substr(plot, 2, 2)))) |> 
  select(plot, block, treat, everything())

saveRDS(stripe_clean, here("DataProcessed/experimental/stripe_clean.rds"))
