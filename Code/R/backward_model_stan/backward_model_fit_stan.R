## ---------------------------
##
## Script name: Local_BackwardModelFit_ForwardInits.
##
## Purpose of script: Fit the backward source-prediction model locally
##                    using multi-core parallelization initialized with
##                    global MLEs from the forward model.
##
## Author: Trent VanHawkins
##
## ---------------------------

options(scipen = 6, digits = 4)
library(dplyr)
library(here)
library(data.table)
library(cmdstanr)
library(furrr)
library(purrr)

# Read in forward fits and experimental baseline matrices -------------------
forward <- readRDS(here(
  "DataProcessed/results/forward_model/stan_forward_fits_multistart.rds"
))
mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat_arrays.rds"))

# Reshape forward estimates from long format into a nested parameter-list format
forward_theta <- forward[,
  .(theta = list(structure(estimate, names = variable))),
  by = .(block, treat, visit)
]

# Set up Grid Combinations --------------------------------------------------
blocks <- dimnames(mod_dat$intensity)[["block"]]
treats <- dimnames(mod_dat$intensity)[["treat"]]
visits <- dimnames(mod_dat$intensity)[["visit"]]
configs <- dimnames(mod_dat$groups)[["config"]]

combos_backward <- expand.grid(
  config = configs,
  blk = blocks,
  trt = treats,
  vst = visits[-1],
  stringsAsFactors = FALSE
)
combos_backward <- left_join(
  combos_backward,
  forward_theta,
  by = c("blk" = "block", "trt" = "treat", "vst" = "visit")
) %>%
  filter(trt == 1)

# Compile Stan Backward Engine ---------------------------------------------
mod_backward <- cmdstan_model(here(
  "Code/R/backward_model_stan/stripe_rust_backwards.stan"
))

# Core Execution Wrapper Function -------------------------------------------
run_local_transition <- function(
  config,
  blk,
  trt,
  vst,
  theta_forward,
  data_list
) {
  # Extract Target Matrix Proportions
  y_current <- data_list$intensity[, blk, trt, vst]
  y_prev <- data_list$intensity[, blk, trt, as.numeric(vst) - 1]
  wind_mat <- data_list$wind[,, blk, trt, vst]
  dist_mat <- data_list$dist
  group_membership <- data_list$groups[, config]

  # Establish Latent Configurations Matrix (K)
  n_src <- as.numeric(trt)
  unique_groups <- sort(unique(group_membership))
  latent_combos <- combn(unique_groups, n_src)
  K_components <- ncol(latent_combos)

  # Build the binary indicator mapping matrix required by Stan (N x K)
  N_nodes <- length(y_current)
  source_in_comp_matrix <- matrix(0, nrow = N_nodes, ncol = K_components)
  for (k in 1:K_components) {
    source_in_comp_matrix[, k] <- as.numeric(
      group_membership %in% latent_combos[, k]
    )
  }

  # Calculate Reference Scaling Factors from Baseline
  sd_auto_ref <- sd(y_prev * (1 - y_prev))
  kernel_ref <- (dist_mat + 0.01)^(-2.0)
  diag(kernel_ref) <- 0
  sd_disp_ref <- sd((wind_mat * kernel_ref) %*% y_prev)

  # Bundle Data List for Stan
  stan_data <- list(
    N = N_nodes,
    K = K_components,
    y = y_current,
    y_lag = y_prev,
    dist_mat = dist_mat,
    wind_mat = wind_mat,
    source_in_component = source_in_comp_matrix,
    sd_auto = sd_auto_ref,
    sd_disp_fixed = sd_disp_ref
  )

  # Construct initialization list mapping forward MLEs to backward names
  init_list <- list(
    beta = as.numeric(theta_forward["beta"]),
    delta = as.numeric(theta_forward["delta"]),
    gamma = as.numeric(theta_forward["gamma"]),
    kappa = as.numeric(theta_forward["kappa"]),
    phi = as.numeric(theta_forward["phi"]),
    alpha = as.numeric(theta_forward["alpha"]),
    pi_vec = as.numeric(rep(1 / K_components, K_components))
  )

  # Target Optimization
  fit_opt <- tryCatch(
    {
      # pass native flags to optimize() to completely suppress logs without capturing stdout
      mod_backward$optimize(
        data = stan_data,
        init = list(init_list),
        seed = 12345,
        refresh = 0,
        show_exceptions = FALSE
      )
    },
    error = function(e) {
      NULL
    }
  )

  # Check if fit_opt is a valid object before pulling parameters
  if (!is.null(fit_opt) && inherits(fit_opt, "CmdStanMLE")) {
    p_mat_extracted <- matrix(
      fit_opt$mle("p_mat"),
      nrow = N_nodes,
      ncol = K_components
    )
    pi_estimates <- fit_opt$mle("pi_vec")
    theta_estimates <- fit_opt$mle(c(
      "beta",
      "delta",
      "gamma",
      "kappa",
      "phi",
      "alpha"
    ))

    return(data.table(
      config = config,
      block = blk,
      treat = trt,
      visit = as.numeric(vst),
      em_iters = NA,
      converged = TRUE,
      Q_final = fit_opt$lp(),
      theta = list(theta_estimates),
      p_mat = list(p_mat_extracted),
      pi = list(pi_estimates)
    ))
  } else {
    return(data.table(
      config = config,
      block = blk,
      treat = trt,
      visit = as.numeric(vst),
      em_iters = NA,
      converged = FALSE,
      Q_final = -Inf,
      theta = list(NA),
      p_mat = list(NA),
      pi = list(NA)
    ))
  }
}

# ==============================================================================
# SINGLE-TRANSITION ISOLATED TESTING BLOCK
# ==============================================================================
message("Running an isolated single-transition test...")

test_row_idx <- 44
test_combo <- combos_backward[test_row_idx, ]

message(sprintf(
  "Testing - Config: %s | Block: %s | Treat: %s | Visit: %s",
  test_combo$config,
  test_combo$blk,
  test_combo$trt,
  test_combo$vst
))

test_start_time <- Sys.time()

single_test_result <- run_local_transition(
  config = test_combo$config,
  blk = test_combo$blk,
  trt = test_combo$trt,
  vst = test_combo$vst,
  theta_forward = test_combo$theta[[1]],
  data_list = mod_dat
)

test_end_time <- Sys.time()
message(
  "Test transition completed in ",
  round(difftime(test_end_time, test_start_time, units = "secs"), 2),
  " seconds."
)
print(single_test_result)

if (single_test_result$converged) {
  test_p_mat <- single_test_result$p_mat[[1]]
  test_pi <- single_test_result$pi[[1]]
  test_theta <- single_test_result$theta[[1]]

  message("\n--- Diagnostic Structural Sanity Checks ---")
  message(
    "Dimensions of Posterior Matrix (p_mat) [Should be N x K]: ",
    paste(dim(test_p_mat), collapse = " x ")
  )
  message(
    "Number of Mixture Weights (pi) [Should match K components]: ",
    length(test_pi)
  )
  message(
    "Sum of Mixture Weights (pi) [Should equal exactly 1.0]: ",
    sum(test_pi)
  )

  message("\nOptimized Backward Parameter Targets:")
  print(test_theta)
} else {
  warning("The single test optimization path failed to converge.")
}
# ==============================================================================

# Local Sequential Execution Setup ------------------------------------------
message(
  "\nStarting local sequential backward modeling initialized via forward MLEs..."
)
message(
  "Transitions will process one by one. Internal Stan optimization tracking enabled."
)

start_time <- Sys.time()

# Standard pmap keeps everything on your main R thread with a native progress bar
results_list <- pmap(
  list(
    combos_backward$config,
    combos_backward$blk,
    combos_backward$trt,
    combos_backward$vst,
    combos_backward$theta
  ),
  ~ run_local_transition(
    config = ..1,
    blk = ..2,
    trt = ..3,
    vst = ..4,
    theta_forward = ..5,
    data_list = mod_dat
  ),
  .progress = TRUE
)

end_time <- Sys.time()
message(
  "All transitions completed in ",
  round(difftime(end_time, start_time, units = "mins"), 2),
  " minutes."
)

# Bind results cleanly into an active session master data.table
backward_master_results <- rbindlist(results_list, fill = TRUE)

# Print a final execution completion check summary
print(backward_master_results[, .(
  total_runs = .N,
  success_count = sum(converged)
)])
