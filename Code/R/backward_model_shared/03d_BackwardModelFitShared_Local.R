## ---------------------------
##
## Script name: 03d_BackwardModelFit_Local.R
##
## Purpose of script: Fit the backward source-prediction model (local parallel version)
##
## Author: Trent VanHawkins
##
## Date Created: 2025-08-20
##
##
## ---------------------------
options(scipen = 6, digits = 4)
library(dplyr)
library(here)
library(parallel)
library(data.table)

## Read in necessary functions
source(here("Code/R/backward_model_shared/03a_BackwardGradFunShared.R"))
source(here("Code/R/backward_model_shared/03b_BackwardModelFunShared.R"))

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

# Fit all combinations in parallel ----------------------------------------
start <- Sys.time()

results <- mclapply(
  seq_len(nrow(combos_backward)),
  function(i) {
    combo <- combos_backward[i, ]

    tryCatch(
      {
        backward_result <- backward_fit(
          config = combo$config,
          blk = combo$block,
          trt = combo$treat,
          vst = combo$visit,
          n_src = combo$n_src,
          inits = combo$theta[[1]],
          mod_dat = mod_dat,
          tol = 1e-4,
          max_iter = 200
        )

        predictions <- if (
          backward_result$converged && !is.null(backward_result$p_mat)
        ) {
          source_pred(
            config = combo$config,
            blk = combo$block,
            trt = combo$treat,
            vst = combo$visit,
            n_src = combo$n_src,
            p_mat = backward_result$p_mat[[1]],
            mod_dat = mod_dat
          )
        } else {
          NULL
        }

        backward_result <- backward_result[, !c("p_mat")]

        result <- if (!is.null(predictions)) {
          merge(
            backward_result,
            predictions,
            by = c("config", "block", "treat", "visit", "n_src")
          )
        } else {
          backward_result[, (prediction_cols) := NA]
          backward_result
        }

        # Report progress in batches of 50
        if (i %% 50 == 0) {
          message("Completed ", i, " of ", nrow(combos_backward), " rows")
        }

        result
      },
      error = function(e) {
        message("Task ", i, " failed: ", conditionMessage(e))
        NULL
      }
    )
  },
  mc.cores = parallel::detectCores() - 1
)

runtime <- difftime(Sys.time(), start, units = "mins")
message("Total runtime: ", round(runtime, 2), " minutes")

# Report failed tasks
failed <- which(sapply(results, is.null))
if (length(failed) > 0) {
  message("Failed tasks: ", paste(failed, collapse = ", "))
  message("Failed combinations:")
  print(combos_backward[
    failed,
    c("config", "block", "treat", "visit", "n_src")
  ])
}

# Combine and save
final_results <- rbindlist(results[!sapply(results, is.null)])

output_dir <- here("DataProcessed/results/backward_model/")
saveRDS(final_results, file.path(output_dir, "backward_fits_shared.rds"))

message("Done. ", nrow(final_results), " combinations saved successfully.")
