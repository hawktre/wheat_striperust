library(here)
library(tidyverse)
library(data.table)

source(here("Code/R/backward_model_shared/03b_BackwardModelFunShared.R"))

mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat_arrays.rds"))
forward <- readRDS(here("DataProcessed/results/forward_model/forward_fits.rds"))
backward_individual <- readRDS(here("DataProcessed/results/backward_model/backward_fits_individual.rds"))
backward_shared <- readRDS(here("DataProcessed/results/backward_model/backward_fits_shared.rds"))
backward_sensitivity <- readRDS(here("DataProcessed/results/backward_model/backward_fits_sensitivity_sequential.rds"))
# Add labels to results and merge them -------------------------------------------------
backward_individual <- backward_individual |> mutate(result_type = "individual")
backward_shared <- backward_shared |> mutate(result_type = "shared")
backward_sensitivity <- backward_sensitivity |> 
  group_by(config, block, treat, visit) |> 
  slice_min(Q_final) |> 
  ungroup() |> 
  mutate(result_type = "sensitivity")

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

shared_no64 <- backward_shared %>% 
  filter(treat %in% c("1","2") & config %in% c("4", "16", "64")) %>% 
  select(config, block, treat, visit, n_src, p_mat) %>% 
  pmap(~source_pred(config = ..1,
                    blk = ..2,
                    trt = ..3,
                    vst = ..4,
                    n_src = ..5,
                    p_mat = ..6,
                    mod_dat = mod_dat)) %>% 
  rbindlist()

sensitivity_predictions <- backward_sensitivity |> 
  select(config, block, treat, visit, n_src, p_mat) %>% 
  pmap(~source_pred(config = ..1,
                    blk = ..2,
                    trt = ..3,
                    vst = ..4,
                    n_src = ..5,
                    p_mat = ..6,
                    mod_dat = mod_dat)) %>% 
  rbindlist()

backward_individual <- left_join(backward_individual%>% select(-p_mat), individual_predictions %>% select(-p_bar), by = c("config", "block", "treat", "visit", "n_src"))
backward_shared <- left_join(backward_shared %>% select(-p_mat), shared_predictions %>% select(-p_bar), by = c("config", "block", "treat", "visit", "n_src"))
backward_sensitivity <- left_join(backward_sensitivity%>% select(-p_mat), sensitivity_predictions %>% select(-p_bar), by = c("config", "block", "treat", "visit", "n_src"))

backward <- list("individual" = backward_individual,
"shared" = backward_shared,
"shared_short" = shared_no64,
"sensitivity" = backward_sensitivity)

saveRDS(backward, here("DataProcessed/results/backward_model/all_backward_fits.rds"))
