## ---------------------------
##
## Script name: 04_ModelFit_Hurdle_Logistic.R
##
## Purpose of script: Fit model with logistic autoinfection term
##
## Author: Trent VanHawkins
##
## Date Created: 2025-05-16
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
library(here)
library(tidyverse)
library(furrr)
library(data.table)
library(MASS)

# Read in the data --------------------------------------------------------
mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat_arrays.rds"))
source(here("Code/R/forward_model/02b_ForwardModelFun.R"))


# Set up indices ----------------------------------------------------------
blocks <- dimnames(mod_dat$intensity)[["block"]]
treats <- dimnames(mod_dat$intensity)[["treat"]]
visits <- dimnames(mod_dat$intensity)[["visit"]]
kappa_try <- exp(seq(log(0.5), log(4.0), length.out = 8))

combos <- expand.grid(
  block = blocks,
  treat = treats,
  visit = visits[2:length(visits)],
  stringsAsFactors = FALSE
)

# Fit the model -----------------------------------------------------------

start <- Sys.time()
free_fits <- pmap(
  combos,
  ~ forward_fit(..1, ..2, ..3, mod_dat, kappa_try)
) %>%
  rbindlist()
end <- Sys.time()
runtime <- difftime(end, start, units = "mins") # could be "mins", "hours", etc.
message("Runtime = ", round(runtime, 2), " minutes")


# Get Fitted values & residuals-----------------------------------------------------------
free_fits$fitted <- pmap(
  free_fits,
  ~ get_fitted(
    blk = ..1,
    trt = ..2,
    vst = ..3,
    par = ..7,
    alpha = ..8,
    mod_dat
  )
)
free_fits$resid <- pmap(
  free_fits,
  ~ compute_deviance_resid(
    blk = ..1,
    trt = ..2,
    vst = ..3,
    par = ..7,
    alpha = ..8,
    fitted = ..9,
    data = mod_dat
  )
)

saveRDS(free_fits, here("DataProcessed/results/forward_model/forward_fits.rds"))
