## ---------------------------
##
## Script name: 05c_EMFormat.R
##
## Purpose of script: Add grid metadata for backward model
##
## Author: Trent VanHawkins
##
## Date Created: 2025-06-11
##
##
## ---------------------------

## view outputs in non-scientific notation

options(scipen = 6, digits = 4) 

## ---------------------------

## load up the packages we will need:  (uncomment as required)

library(tidyverse)
library(here)

# Read in the data and functions--------------------------------------------------------
clusters <- readRDS(here("DataProcessed/experimental/clusters.rds"))
forward_dat <- readRDS(here("DataProcessed/experimental/mod_dat_arrays.rds"))

# Write loop to create data packets ---------------------------------------
em_dat <- list()

forward_dat$

saveRDS(em_dat, here("DataProcessed/experimental/em_dat.rds"))
