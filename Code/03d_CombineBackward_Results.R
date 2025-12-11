## ---------------------------
##
## Script name: 03d_CombineBackwardResults.R
##
## Purpose of script: Combine backward model array results
##
## ---------------------------

library(dplyr)
library(purrr)
library(here)
library(data.table)

source(here("Code/03b_BackwardModelFun.R"))

forward <- readRDS(here("DataProcessed/results/forward_model/forward_fits.rds"))
mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat_arrays.rds"))

# Read all individual results
array_dir <- here("DataProcessed/results/backward_model/array_results")
result_files <- list.files(array_dir, pattern = "^backward_.*\\.rds$", full.names = TRUE)

message("Found ", length(result_files), " result files")

backward <- map(result_files, readRDS) %>% rbindlist()

# Source prediction for treat == 1
sources_predicted <- backward %>%
  select(config, block, treat, visit, p_mat, n_src) %>%
  pmap(~source_pred(config = ..1,
                    blk = ..2,
                    trt = ..3,
                    vst = ..4,
                    p_mat = ..5,
                    n_src = ..6,
                    mod_dat = mod_dat)) %>%
  rbindlist()

results_merge <- left_join(backward |> mutate(treat = as.factor(treat),
                                               visit = as.factor(visit)), 
                           forward |> mutate(treat = as.factor(treat),
                                            visit = as.factor(visit)), 
                           by = c("block", "treat", "visit"), 
                           suffix = c(".backward", ".forward")) %>%
  left_join(sources_predicted |> mutate(treat = as.factor(treat),
                                        visit = as.factor(visit)), 
            by = c("config", "block", "treat", "visit"))

saveRDS(results_merge, here("DataProcessed/results/backward_model/backward_fits.rds"))

message("Successfully combined all results")