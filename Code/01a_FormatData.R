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


# Read in the data --------------------------------------------------------
stripe <- readRDS(here("DataProcessed/experimental/stripe_clean.rds"))
wind <- readRDS(here("DataProcessed/wind/wind_clean.rds"))


# Create lagged data frame ------------------------------------------------
stripe %>% 
  

