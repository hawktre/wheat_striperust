## ---------------------------
##
## Purpose of script: Fit forward epidemic model via multi-start optimization,
##                    extract fitted values and deviance residuals, and
##                    generate comprehensive diagnostic plots.
##
## Author: Trent VanHawkins
##
## ---------------------------

library(here)
library(tidyverse)
library(cmdstanr)
library(data.table)

# Read in the data --------------------------------------------------------
mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat_arrays.rds"))

# Set up indices ----------------------------------------------------------
blocks <- dimnames(mod_dat$intensity)[["block"]]
treats <- dimnames(mod_dat$intensity)[["treat"]]
visits <- dimnames(mod_dat$intensity)[["visit"]]

combos <- expand.grid(
  block = blocks,
  treat = treats,
  visit = visits[2:length(visits)],
  stringsAsFactors = FALSE
)

# Stan Model Setup --------------------------------------------------------
mod <- cmdstan_model(
  here("Code/R/forward_model/stripe_rust_mle.stan")
)

# Function to run multi-start optimization for a single transition ----------
fit_stan_transition_multistart <- function(
  blk,
  trt,
  vst,
  data_list,
  n_starts = 5,
  lp_tol = 0.01,
  base_seed = 12345
) {
  # Extract data for this transition
  y_current <- data_list$intensity[, blk, trt, vst]
  y_prev <- data_list$intensity[, blk, trt, as.numeric(vst) - 1]
  wind_mat <- data_list$wind[,, blk, trt, vst]
  dist_mat <- data_list$dist

  N <- length(y_current)

  stan_data <- list(
    N = N,
    y = y_current,
    y_lag = y_prev,
    dist_mat = dist_mat,
    wind_mat = wind_mat,

    # Vague Hyperparameters for priors
    beta_mu = 0,
    beta_sigma = 2.5,
    delta_mu = 0,
    delta_sigma = 2.5,
    gamma_mu = 0,
    gamma_sigma = 2.5,
    kappa_mu = 1,
    kappa_sigma = 1,
    phi_alpha = 3,
    phi_beta = 0.1
  )

  # Storage objects for looping
  best_fit <- NULL
  best_lp <- -Inf
  all_lps <- numeric(n_starts)

  # Execute optimization loop using shifting random seeds
  for (i in 1:n_starts) {
    run_seed <- base_seed + i

    fit_opt <- tryCatch(
      {
        capture.output(
          {
            fit <- mod$optimize(
              data = stan_data,
              seed = run_seed,
              refresh = 0
            )
          },
          type = "output"
        )
        fit
      },
      error = function(e) {
        NULL
      }
    )

    if (!is.null(fit_opt)) {
      current_lp <- fit_opt$lp()
      all_lps[i] <- current_lp

      # Retain track of global maximum
      if (current_lp > best_lp) {
        best_lp <- current_lp
        best_fit <- fit_opt
      }
    } else {
      all_lps[i] <- -Inf
    }
  }

  # Evaluate Convergence: Count hits on the highest log-posterior peak
  valid_lps <- all_lps[is.finite(all_lps)]

  if (length(valid_lps) > 0) {
    runs_matching_best <- sum(abs(valid_lps - best_lp) <= lp_tol)
    global_peak_converged <- (runs_matching_best > 1)
  } else {
    runs_matching_best <- 0
    global_peak_converged <- FALSE
  }

  # Process and pull parameters from the optimal run
  if (!is.null(best_fit)) {
    capture.output(
      {
        # Safely pull exact global scalar MLE parameters via direct vector extraction
        mle_estimates <- best_fit$mle(c(
          "beta",
          "delta",
          "gamma",
          "kappa",
          "phi",
          "alpha"
        ))

        # Pull numeric MLE mu and deviance residual vectors directly by name
        mu_vectors <- best_fit$mle("mu")
        dev_res_vectors <- best_fit$mle("deviance_residual")
      },
      type = "output"
    )

    # Build summary data.table explicitly from vectors
    summ <- data.table(
      variable = names(mle_estimates),
      estimate = as.numeric(mle_estimates),
      block = blk,
      treat = trt,
      visit = vst,
      best_lp = best_lp,
      starts_matching = runs_matching_best,
      is_reliable_peak = global_peak_converged
    )

    # Map out the long node-level predictions and residuals explicitly
    vals <- data.table(
      node = 1:N,
      block = blk,
      treat = trt,
      visit = vst,
      true_y = y_current,
      fitted_mu = as.numeric(mu_vectors),
      dev_res = as.numeric(dev_res_vectors)
    )
  } else {
    # Fallback structure if all seeds crashed
    summ <- data.table(
      variable = c("beta", "delta", "gamma", "kappa", "phi", "alpha"),
      estimate = NA_real_,
      block = blk,
      treat = trt,
      visit = vst,
      best_lp = -Inf,
      starts_matching = 0,
      is_reliable_peak = FALSE
    )
    vals <- data.table(
      node = 1:N,
      block = blk,
      treat = trt,
      visit = vst,
      true_y = y_current,
      fitted_mu = NA_real_,
      dev_res = NA_real_
    )
  }

  return(list(summary = summ, values = vals))
}

# Run fits sequentially with purrr's native progress bar --------------------
message("Processing transition combinations...")

results_raw <- pmap(
  combos,
  ~ fit_stan_transition_multistart(..1, ..2, ..3, mod_dat, n_starts = 5),
  .progress = TRUE
)

# Unpack outputs cleanly into active session objects ------------------------
stan_summaries <- rbindlist(lapply(results_raw, `[[`, "summary"), fill = TRUE)
fitted_values <- rbindlist(lapply(results_raw, `[[`, "values"), fill = TRUE)

message(
  "All optimization runs completed. Objects 'stan_summaries' and 'fitted_values' are ready in memory."
)

saveRDS(
  stan_summaries,
  here(
    "DataProcessed/results/forward_model/stan_forward_fits_multistart.rds"
  )
)

# Plotting Function 1: Fitted vs True ---------------------------------------
plot_fitted_vs_true <- function(pred_data, target_block = NULL) {
  plot_df <- if (!is.null(target_block)) {
    pred_data[block == target_block]
  } else {
    pred_data
  }

  g <- ggplot(plot_df, aes(x = true_y, y = fitted_mu)) +
    geom_point(alpha = 0.4, color = "seagreen") +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = "firebrick",
      linewidth = 0.8
    ) +
    facet_grid(visit ~ treat, labeller = label_both) +
    labs(
      title = paste(
        "Fitted vs. True Disease Intensities",
        if (!is.null(target_block)) paste("- Block", target_block) else ""
      ),
      subtitle = "Dashed line indicating ideal 1:1 perfect prediction matching profile",
      x = "Observed Disease Proportion (True Y)",
      y = "Expected Conditional Infection Mean (Fitted Mu)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "grey95", color = NA),
      panel.spacing = unit(1, "lines")
    )
  return(g)
}


# Plotting Function 2: Residuals vs Fitted -----------------------------------
plot_deviance_residuals <- function(pred_data, target_block = NULL) {
  plot_df <- if (!is.null(target_block)) {
    pred_data[block == target_block]
  } else {
    pred_data
  }

  g <- ggplot(plot_df, aes(x = fitted_mu, y = dev_res)) +
    geom_point(alpha = 0.4, color = "steelblue") +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "firebrick",
      linewidth = 0.8
    ) +
    facet_grid(visit ~ treat, labeller = label_both) +
    labs(
      title = paste(
        "Deviance Residuals vs. Fitted Values",
        if (!is.null(target_block)) paste("- Block", target_block) else ""
      ),
      subtitle = "Horizontal line at 0 indicates perfect fit; look for unwanted trends or trumpet shapes",
      x = "Fitted Expected Value (Mu)",
      y = "Sign-Adjusted Deviance Residual"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "grey95", color = NA),
      panel.spacing = unit(1, "lines")
    )
  return(g)
}


# Display Plots -----------------------------------------------------------
# Filter to Block A to keep the grid clean and highly scannable
plot_fitted_vs_true(fitted_values, target_block = "A")
plot_deviance_residuals(fitted_values, target_block = "A")
plot_fitted_vs_true(fitted_values, target_block = "B")
plot_deviance_residuals(fitted_values, target_block = "B")
plot_fitted_vs_true(fitted_values, target_block = "C")
plot_deviance_residuals(fitted_values, target_block = "C")
plot_fitted_vs_true(fitted_values, target_block = "D")
plot_deviance_residuals(fitted_values, target_block = "D")

stan_summaries |> 
  ggplot(aes(x = visit, y = estimate)) +
  geom_point(aes(color = block)) +
  facet_grid(treat ~ variable, scales = "free_y")
