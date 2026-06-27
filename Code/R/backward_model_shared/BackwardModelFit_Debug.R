## ---------------------------
##
## Script name: 03d_BackwardModelFit_Local_DEBUG.R
##
## Purpose of script: Debug version — test one row at a time
##                    Updated for log-scale initialization parsing.
##
## Author: Trent VanHawkins
##
## ---------------------------
options(scipen = 6, digits = 4)
library(dplyr)
library(here)
library(data.table)

## Read in necessary functions
source(here("Code/R/backward_model_shared/03a_BackwardGradFunShared.R"))
source(here("Code/R/backward_model_shared/03b_BackwardModelFunShared.R"))

## Read in forward fits and experimental data
forward <- readRDS(here("DataProcessed/results/forward_model/forward_fits.rds"))
mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat_arrays.rds"))

# Set up Indices ----------------------------------------------------------
blocks <- dimnames(mod_dat$intensity)[["block"]]
treats <- dimnames(mod_dat$intensity)[["treat"]]
visits <- dimnames(mod_dat$intensity)[["visit"]]
configs <- dimnames(mod_dat$groups)[["config"]]
n_src <- c(1, 2, 3, 4)

# Create full combinations
combos_backward <- expand.grid(
  config = configs,
  block = blocks,
  treat = treats,
  n_src = n_src,
  visit = visits[-1],
  stringsAsFactors = FALSE
) |>
  filter(
    !(config == "4" & n_src == 4),
    !(config == "64" & n_src > 3),
    !(config == "64" & treat == "4")
  )

combos_backward <- left_join(
  combos_backward,
  forward %>% select(block, treat, visit, theta),
  by = c("block", "treat", "visit")
)

prediction_cols <- c(
  "mean_error",
  "n_correct",
  "acc",
  "dist_acc",
  "component_dist_acc",
  "predicted_source",
  "true_source",
  "naive_source",
  "naive_mean_error",
  "naive_n_correct",
  "naive_acc",
  "naive_dist_acc",
  "naive_component_dist_acc"
)

# ── DEBUG CONTROLS ────────────────────────────────────────────────────────────
# Set the row index you want to test, then step through the blocks below.
DEBUG_ROW <- 1 # <-- change this to any row in combos_backward
SHOW_COMBO <- TRUE # print the combo being tested
# ─────────────────────────────────────────────────────────────────────────────

# Step 1: Inspect the combo & process initial values -----------------------
combo <- combos_backward[DEBUG_ROW, ]

# Extract and structurally transform initial theta parameters to log-scale
init_pars <- combo$theta[[1]]
log_pars <- c("delta", "gamma", "kappa", "phi")
init_pars[log_pars] <- log(pmax(init_pars[log_pars], 1e-5)) 
names(init_pars) <- c("beta", "log_delta", "log_gamma", "log_kappa", "log_phi")

if (SHOW_COMBO) {
  message(">>> Testing row ", DEBUG_ROW, " of ", nrow(combos_backward))
  print(combo[, c("config", "block", "treat", "visit", "n_src")])
  message("Forward input theta (natural): ", paste(round(combo$theta[[1]], 4), collapse = ", "))
  message("EM transformed input (log-scale): ", paste(round(init_pars, 4), collapse = ", "))
}


# Step 2: Run backward_fit ------------------------------------------------
backward_result <- backward_fit(
  config = combo$config,
  blk = combo$block,
  trt = combo$treat,
  vst = combo$visit,
  n_src = combo$n_src,
  inits = init_pars, # Pass structural log-parameters
  mod_dat = mod_dat,
  tol = 1e-4,
  max_iter = 200
)

message("Converged: ", backward_result$converged)
print(backward_result)


# Step 3: Run source_pred (only if converged and p_mat exists) ------------
predictions <- NULL

if (backward_result$converged && !is.null(backward_result$p_mat)) {
  predictions <- source_pred(
    config = combo$config,
    blk = combo$block,
    trt = combo$treat,
    vst = combo$visit,
    n_src = combo$n_src,
    p_mat = backward_result$p_mat[[1]],
    mod_dat = mod_dat
  )
  message("Predictions:")
  print(predictions)
} else {
  message("Skipping source_pred — model did not converge or p_mat is NULL.")
}


# Step 4: Merge results ---------------------------------------------------
backward_result_clean <- backward_result[, !c("p_mat")]

row_result <- if (!is.null(predictions)) {
  merge(
    backward_result_clean,
    predictions,
    by = c("config", "block", "treat", "visit", "n_src")
  )
} else {
  backward_result_clean[, (prediction_cols) := NA]
  backward_result_clean
}

message("Final merged result for row ", DEBUG_ROW, ":")
print(row_result)
