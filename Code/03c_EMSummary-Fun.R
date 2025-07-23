## ---------------------------
##
## Script name: 03c_EMSummary-Fun
##
## Purpose of script: Write Summaries into functions
##
## Author: Trent VanHawkins
##
## Date Created: 2025-07-22
##
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## ---------------------------

## view outputs in non-scientific notation

options(scipen = 6, digits = 4) 

## ---------------------------

## load up the packages we will need:  (uncomment as required)

library(tidyverse)
library(data.table)
library(here)

mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat_arrays.rds"))
backward_fits <- readRDS(here("DataProcessed/results/backward_model/backwards_fits.rds"))

source_pred <- function(config, blk, trt, vst, p_mat, mod_dat) {

  # Ensure groups is a factor with levels = 1:n_groups
  groups <- as.factor(mod_dat$groups[,config])
  levels(groups) <- sort(unique(groups))
  
  # 1. Compute p̄ matrix (average p_mat rows by destination group)
  group_ids <- sort(unique(groups))
  n_groups <- length(group_ids)
  
  p_bar <- matrix(NA, n_groups, n_groups)
  rownames(p_bar) <- colnames(p_bar) <- as.character(group_ids)
  
  for (g in group_ids) {
    rows_in_group <- which(groups == g)
    p_bar[as.character(g), ] <- colMeans(p_mat[rows_in_group, , drop = FALSE])
  }
  
  # 2. Predict source group for each destination group (row-wise argmax)
  predicted_source <- which(p_bar == max(p_bar), arr.ind = TRUE)[[2]]
  
  # 3. Compare to ground truth
  true_source <- mod_dat$truth[blk, config]  # vector of length n_groups
  
  # 4. Compute distance-weighted accuracy
  dist <- mod_dat$grid_dist[[config]]
  error <- dist[predicted_source, true_source]
  dist_acc <- 1 - (error / max(dist[predicted_source,]))  # normalized distance accuracy
  
  result <- list(
    config = config,
    block = blk,
    treat = trt,
    visit = vst,
    predicted_source = predicted_source,
    true_source = true_source,
    dist_acc = dist_acc
  )
  # 5. Compute other metrics 
  return(as.data.table(result))
}

backward_preds <- pmap(backward_fits[treat == 1], ~source_pred(config = ..1,
                                 blk = ..2,
                                 trt = ..3,
                                 vst = ..4,
                                 p_mat = ..9,
                                 mod_dat)) %>% rbindlist()
