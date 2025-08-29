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
sims <- readRDS(here("DataProcessed/results/simulation/sims.rds"))


sim_9921 <- single_sim(9921, mod_dat, forward_fits, kappa_try = seq(0.25, 2.5, 0.25))
sim_4994 <- single_sim(4994, mod_dat, forward_fits, kappa_try = seq(0.25, 2.5, 0.25))

sims_final <- data.table::rbindlist(list(sims, sim_9921, sim_4994))

saveRDS(sims_final, here("DataProcessed/results/simulation/sims_appended.rds"))
