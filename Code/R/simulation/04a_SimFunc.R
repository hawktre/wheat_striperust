## ---------------------------
##
## Script name: 04a_SimFunc.R
##
## Purpose of script: Compile all necessary functions for simulation.
##
## Author: Trent VanHawkins
##
## Date Created: 2025-07-22
##
##
## ---------------------------
options(scipen = 6, digits = 4)

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

# Function to simulate the data -------------------------------------------
disease_sim <- function(pars, mu, alpha) {
  phi <- exp(pars[["phi"]])

  a <- mu * phi
  b <- phi * (1 - mu)

  min_detectable <- 0.0001
  max_detectable <- 1 - 0.0001

  sim_dat <- map2_dbl(
    a,
    b,
    ~ {
      val <- rbeta(1, shape1 = .x, shape2 = .y)
      pmin(pmax(val, min_detectable), max_detectable)
    }
  )

  zeros <- rbinom(length(sim_dat), size = 1, prob = 1 - alpha)
  return(sim_dat * zeros)
}

# Wrapper Function for the simulation -------------------------------------
single_sim <- function(
  sim_id,
  dat,
  forward_mod,
  kappa_try,
  n_src = c(1, 2, 3, 4),
  output_dir = here("DataProcessed/results/simulation")
) {
  base_seed <- 404
  set.seed(base_seed + sim_id)

  tryCatch(
    {
      # Set up indices
      blocks <- dimnames(dat$intensity)[["block"]]
      treats <- dimnames(dat$intensity)[["treat"]]
      visits <- dimnames(dat$intensity)[["visit"]]
      configs <- dimnames(dat$groups)[["config"]]

      # 1. Simulate new intensity values
      for (blk in blocks) {
        for (trt in treats) {
          for (vst in visits[-1]) {
            fit <- forward_mod[
              block == blk & treat == as.numeric(trt) & visit == as.numeric(vst)
            ]
            pars <- fit[["theta"]][[1]]
            alpha <- fit[["alpha"]]
            fitted <- fit[["fitted"]][[1]]
            dat$intensity[, blk, trt, vst] <- disease_sim(pars, fitted, alpha)
          }
        }
      }

      # 2. Fit Forward Model
      combos_forward <- expand.grid(
        block = blocks,
        treat = treats,
        visit = visits[-1],
        stringsAsFactors = FALSE
      )

      forward <- pmap(
        combos_forward,
        ~ forward_fit(..1, ..2, ..3, dat, kappa_try)
      ) |>
        rbindlist()

      # 3. Fit backward model with all n_src values
      combos_backward <- expand.grid(
        config = configs,
        block = blocks,
        treat = treats,
        n_src = n_src,
        visit = visits[-1],
        stringsAsFactors = FALSE
      ) |>
        filter(
          !(config == "64" & treat > 2),
          !(config == "64" & n_src > 2),
          !(config == "4" & n_src == "4")
        ) |>
        left_join(
          forward |> select(block, treat, visit, theta),
          by = c("block", "treat", "visit")
        )

      # Process each backward combination
      results_list <- lapply(seq_len(nrow(combos_backward)), function(i) {
        combo <- combos_backward[i, ]

        backward_result <- backward_fit(
          config = combo$config,
          blk = combo$block,
          trt = combo$treat,
          vst = combo$visit,
          n_src = combo$n_src,
          inits = combo$theta[[1]],
          mod_dat = dat,
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
            mod_dat = dat
          )
        } else {
          NULL
        }

        backward_result <- backward_result[, !c("p_mat")]

        if (!is.null(predictions)) {
          merge(
            backward_result,
            predictions,
            by = c("config", "block", "treat", "visit", "n_src")
          )
        } else {
          backward_result[, (prediction_cols) := NA]
          backward_result
        }
      })

      # Combine all results
      backward_results <- rbindlist(results_list)

      # Merge with forward results and add simulation ID
      final_results <- merge(
        backward_results,
        forward |> mutate(visit = as.numeric(visit)),
        by = c("block", "treat", "visit"),
        suffix = c(".backward", ".forward")
      ) |>
        mutate(sim = sim_id) |>
        select(sim, everything())

      return(final_results)
    },
    error = function(e) {
      err_file <- file.path(output_dir, sprintf("sim_error_%s.txt", sim_id))
      cat(
        "Simulation failed:\n",
        "sim_id: ",
        sim_id,
        "\n",
        "message: ",
        conditionMessage(e),
        "\n",
        "class: ",
        paste(class(e), collapse = ", "),
        "\n",
        "call: ",
        deparse(conditionCall(e)),
        "\n",
        file = err_file,
        append = FALSE,
        sep = ""
      )

      dump_name <- file.path(output_dir, sprintf("sim_dump_%s", sim_id))
      dump.frames(dump_name, to.file = TRUE)

      cat(
        "\nDumped frames to: ",
        dump_name,
        ".rda\n",
        file = err_file,
        append = TRUE
      )

      stop(e)
    }
  )
}
