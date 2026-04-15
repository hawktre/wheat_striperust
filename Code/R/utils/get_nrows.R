library(dplyr)
library(here)

## Read in forward fits
forward <- readRDS(here("DataProcessed/results/forward_model/forward_fits.rds"))
mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat_arrays.rds"))

# Set up Indices ----------------------------------------------------------
blocks <- dimnames(mod_dat$intensity)[["block"]]
treats <- dimnames(mod_dat$intensity)[["treat"]]
visits <- dimnames(mod_dat$intensity)[["visit"]]
configs <- dimnames(mod_dat$groups)[["config"]]
n_src <- c(1, 2, 3, 4)

# Create full combinations
combos_backward <- expand.grid(config = configs, blk = blocks, trt = treats, n_src = n_src, vst = visits[-1], stringsAsFactors = FALSE)
combos_backward <- left_join(combos_backward, forward %>% select(block, treat, visit, theta), 
                             by = c("blk"="block", "trt"="treat", "vst"="visit")) 


cat(nrow(combos_backward))