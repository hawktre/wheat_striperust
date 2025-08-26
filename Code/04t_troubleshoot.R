library(here)
library(data.table)
library(dplyr)
library(purrr)
library(parallel)  # Use base parallel package
source(here("Code/02b_ForwardModelFun.R"))
source(here("Code/03b_BackwardModelFun.R"))
source(here("Code/04a_SimFunc.R"))

# Read in data
forward_fits <- readRDS(here("DataProcessed/results/forward_model/forward_fits.rds"))
mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat_arrays.rds"))



test <- single_sim(84, mod_dat, forward_fits, kappa_try = seq(0.25, 2.5, 0.25))


