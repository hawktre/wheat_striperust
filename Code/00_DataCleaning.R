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

## view outputs in non-scientific notation

options(scipen = 6, digits = 4) 

## ---------------------------

## load up the packages we will need:  (uncomment as required)

library(tidyverse)
library(here)
library(readxl)



# Read in the data --------------------------------------------------------
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
             visit = paste0("visit", rep(1:5, length.out = n()))))


