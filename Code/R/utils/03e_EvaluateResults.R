library(here)
library(tidyverse)
library(data.table)

source(here("Code/R/backward_model_individual/03b_BackwardModelFun.R"))

mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat_arrays.rds"))
forward <- readRDS(here("DataProcessed/results/forward_model/forward_fits.rds"))
backward_individual <- readRDS(here("DataProcessed/results/backward_model/backward_fits.rds"))
backward_shared <- readRDS(here("DataProcessed/results/backward_model/backward_fits_shared.rds"))

test <- backward_shared |> mutate(config = factor(config, levels = c("4", "8h", "8v", "16", "64"))) |> arrange(visit, treat, block, config)
# Add labels to results and merge them --------------------------------------------------
backward_individual <- backward_individual |> mutate(result_type = "individual")
backward_shared <- backward_shared |> mutate(result_type = "shared")


individual_predictions <- backward_individual %>% 
  select(config, block, treat, visit, n_src, p_mat) %>% 
  pmap(~source_pred(config = ..1,
                    blk = ..2,
                    trt = ..3,
                    vst = ..4,
                    n_src = ..5,
                    p_mat = ..6,
                    mod_dat = mod_dat)) %>% 
  rbindlist()

shared_predictions <- backward_shared %>% 
  select(config, block, treat, visit, n_src, p_mat) %>% 
  pmap(~source_pred(config = ..1,
                    blk = ..2,
                    trt = ..3,
                    vst = ..4,
                    n_src = ..5,
                    p_mat = ..6,
                    mod_dat = mod_dat)) %>% 
  rbindlist()

backward_individual <- left_join(backward_individual, individual_predictions, by = c("config", "block", "treat", "visit", "n_src"))
backward_shared <- left_join(backward_shared, shared_predictions, by = c("config", "block", "treat", "visit", "n_src"))

backward <- rbind(backward_individual, backward_shared)
saveRDS(backward, here("DataProcessed/results/backward_model/backward_fits_eval.rds"))
