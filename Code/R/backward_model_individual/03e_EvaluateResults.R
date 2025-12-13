library(here)
library(tidyverse)
library(data.table)

source(here("Code/R/backward_model/03b_BackwardModelFun.R"))

mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat_arrays.rds"))
forward <- readRDS(here("DataProcessed/results/forward_model/forward_fits.rds"))
backward <- readRDS(here("DataProcessed/results/backward_model/backward_fits.rds"))

sources_predicted <- backward %>% 
  select(config, block, treat, visit, n_src, p_mat) %>% 
  pmap(~source_pred(config = ..1,
                    blk = ..2,
                    trt = ..3,
                    vst = ..4,
                    n_src = ..5,
                    p_mat = ..6,
                    mod_dat = mod_dat)) %>% 
  rbindlist()

results_merge <- left_join(backward |> mutate(treat = as.factor(treat), 
                                              visit = as.factor(visit)), forward |> mutate(treat = as.factor(treat), 
                                                                                           visit = as.factor(visit)), by = c("block", "treat", "visit"), suffix = c(".backward", ".forward")) %>% 
  left_join(sources_predicted |> mutate(treat = as.factor(treat), 
                                        visit = as.factor(visit)), by = c("config", "block", "treat", "visit")) 

saveRDS(results_merge, here("DataProcessed/results/backward_model/backward_fits_full.rds"))
