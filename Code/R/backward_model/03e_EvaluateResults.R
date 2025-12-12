library(here)
library(tidyverse)

source(here("Code/R/backward_model/03b_BackwardModelFun.R"))

backward <- readRDS(here("DataProcessed/results/backward_model/backward_fits.rds"))

sources_predicted <- backward %>% 
  select(config, block, treat, visit, p_mat) %>% 
  pmap(~source_pred(config = ..1,
                    blk = ..2,
                    trt = ..3,
                    vst = ..4,
                    p_mat = ..5,
                    mod_dat = mod_dat)) %>% 
  rbindlist()

results_merge <- left_join(backward |> mutate(treat = as.factor(treat), 
                                              visit = as.factor(visit)), forward |> mutate(treat = as.factor(treat), 
                                                                                           visit = as.factor(visit)), by = c("block", "treat", "visit"), suffix = c(".backward", ".forward")) %>% 
  left_join(sources_predicted |> mutate(treat = as.factor(treat), 
                                        visit = as.factor(visit)), by = c("config", "block", "treat", "visit")) 

saveRDS(results_merge, here("DataProcessed/results/backward_model/backward_fits_sensitivity.rds"))