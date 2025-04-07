## ---------------------------
##
## Script name: 00a_TrevData_Exploaration.R
##
## Purpose of script: Explore archived data and code from Trevor
##
## Author: Trent VanHawkins
##
## Date Created: 2025-02-24
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

#Load data from the archive
cdm_data <- load(here("Code/Archive/cdm_project/cdm_data.RData"))
