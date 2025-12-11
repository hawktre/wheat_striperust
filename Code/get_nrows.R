library(dplyr)
library(here)

forward <- readRDS(here("DataProcessed/results/forward_model/forward_fits.rds"))
mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat_arrays.rds"))

blocks <- dimnames(mod_dat$intensity)[["block"]]
treats <- dimnames(mod_dat$intensity)[["treat"]]
visits <- dimnames(mod_dat$intensity)[["visit"]]
configs <- dimnames(mod_dat$groups)[["config"]]

combos_backward <- expand.grid(config = configs, blk = blocks, trt = treats, vst = visits[-1], stringsAsFactors = FALSE)
combos_backward <- left_join(combos_backward, forward %>% select(block, treat, visit, theta), 
                             by = c("blk"="block", "trt"="treat", "vst"="visit")) |> 
  filter(!(config == "64" & trt %in% c(2,4)))

cat(nrow(combos_backward))