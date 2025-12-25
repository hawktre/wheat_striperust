library(here)
library(tidyverse)

## Read in forward fits
mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat_arrays.rds"))

# Set up Indices ----------------------------------------------------------
blocks <- dimnames(mod_dat$intensity)[["block"]]
treats <- dimnames(mod_dat$intensity)[["treat"]]
visits <- dimnames(mod_dat$intensity)[["visit"]]
configs <- dimnames(mod_dat$groups)[["config"]]
kappa_try <- c(0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0)


# Get the inits for each value of kappa -----------------------------------

# Fit the model -----------------------------------------------------------
# 2. Fit backward model
combos_backward <- expand.grid(config = configs, blk = blocks, trt = treats, vst = visits[-1], kappa = kappa_try, stringsAsFactors = FALSE)
combos_backward_t1 <- combos_backward |> 
  filter(trt == 1) |> 
  arrange(config, blk, trt, vst, kappa) |> 
  group_by(config, blk, trt, vst) |> 
  mutate(batch = cur_group_id()) |> 
  ungroup()

cat(max(combos_backward_t1$batch))
