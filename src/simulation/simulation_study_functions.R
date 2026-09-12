# Functions for the whole-epidemic source-detection simulation study.

source(here::here(
  "src", "source_detection", "source_detection_functions.R"
))
source(here::here(
  "src", "simulation", "whole_epidemic_pilot_functions.R"
))

simulation_configs <- function(treat, available_configs) {
  requested <- if (as.integer(treat) == 1L) {
    c("4", "8h", "8v", "16", "64")
  } else {
    c("4", "8h", "8v", "16")
  }
  intersect(requested, available_configs)
}

summarize_forward_simulation_fit <- function(
  forward_fit, generating_fit, block, treat, simulation
) {
  estimates <- forward_fit$summary
  truth <- generating_fit$summary |>
    dplyr::select(
      visit,
      true_beta = beta,
      true_delta = delta,
      true_gamma = gamma,
      true_kappa = kappa,
      true_phi = phi
    )

  optimizer <- dplyr::bind_rows(lapply(seq_along(forward_fit$fits), function(i) {
    selected <- forward_fit$fits[[i]]
    fit <- selected$fit
    data.frame(
      visit = rownames(forward_fit$theta)[i],
      function_evaluations = unname(fit$counts[["function"]]),
      gradient_evaluations = unname(fit$counts[["gradient"]]),
      optimizer_message = if (is.null(fit$message)) NA_character_ else fit$message
    )
  }))

  estimates |>
    dplyr::mutate(visit = as.character(visit)) |>
    dplyr::left_join(
      truth |> dplyr::mutate(visit = as.character(visit)),
      by = "visit"
    ) |>
    dplyr::left_join(optimizer, by = "visit") |>
    dplyr::mutate(
      simulation = simulation,
      block = block,
      treat = as.character(treat),
      beta_bias = beta - true_beta,
      delta_bias = delta - true_delta,
      gamma_bias = gamma - true_gamma,
      kappa_bias = kappa - true_kappa,
      phi_bias = phi - true_phi,
      log_kappa_bias = log(kappa) - log(true_kappa),
      log_phi_bias = log(phi) - log(true_phi),
      beta_wald_covered = abs(beta_bias) <= 1.96 * se_beta,
      delta_wald_covered = abs(delta_bias) <= 1.96 * se_delta,
      gamma_wald_covered = abs(gamma_bias) <= 1.96 * se_gamma,
      kappa_wald_covered = abs(log_kappa_bias) <= 1.96 * se_log_kappa,
      phi_wald_covered = abs(log_phi_bias) <= 1.96 * se_log_phi,
      beta_standardized_bias = beta_bias / se_beta,
      delta_standardized_bias = delta_bias / se_delta,
      gamma_standardized_bias = gamma_bias / se_gamma,
      log_kappa_standardized_bias = log_kappa_bias / se_log_kappa,
      log_phi_standardized_bias = log_phi_bias / se_log_phi,
      .before = 1
    )
}

summarize_backward_fit <- function(fit, simulation) {
  final_history <- fit$history[nrow(fit$history), , drop = FALSE]
  candidate <- fit$candidate_summary
  true_row <- candidate[candidate$true_candidate, , drop = FALSE]
  winner <- candidate[which.max(candidate$posterior), , drop = FALSE]

  diagnostics <- data.frame(
    simulation = simulation,
    block = winner$block,
    treat = as.character(winner$treat),
    config = winner$config,
    n_sources = winner$n_sources,
    through_visit = winner$through_visit,
    fit_scope = winner$fit_scope,
    status = fit$parameter_summary$status[1],
    iterations = final_history$iteration,
    converged = isTRUE(final_history$converged),
    observed_loglik = final_history$observed_loglik,
    final_loglik_change = final_history$loglik_change,
    final_parameter_change = final_history$max_parameter_change,
    final_responsibility_change = final_history$max_responsibility_change,
    final_m_step_convergence = final_history$m_step_convergence,
    all_iterations_monotonic = all(
      fit$history$monotonic[!is.na(fit$history$monotonic)]
    ),
    any_m_step_failure = any(
      fit$history$m_step_convergence > 0,
      na.rm = TRUE
    ),
    predicted_candidate_id = winner$candidate_id,
    true_candidate_id = winner$true_candidate_id,
    true_source_set_selected = winner$candidate_id == winner$true_candidate_id,
    winning_posterior = winner$posterior,
    true_posterior = if (nrow(true_row)) true_row$posterior else NA_real_,
    true_rank = if (nrow(true_row)) true_row$rank else NA_real_,
    posterior_entropy = -sum(ifelse(
      candidate$posterior > 0,
      candidate$posterior * log(candidate$posterior),
      0
    ))
  )
  diagnostics$effective_candidates <- exp(diagnostics$posterior_entropy)

  accuracy <- fit$accuracy_summary |>
    dplyr::mutate(
      simulation = simulation,
      treat = as.character(treat),
      .before = 1
    )
  parameters <- fit$parameter_summary |>
    dplyr::mutate(
      simulation = simulation,
      treat = as.character(treat),
      .before = 1
    )

  list(diagnostics = diagnostics, accuracy = accuracy, parameters = parameters)
}

run_backward_simulation_fits <- function(
  simulated_mod_dat,
  forward_fit,
  block,
  treat,
  simulation,
  configs
) {
  visits <- dimnames(simulated_mod_dat$intensity)[["visit"]][-1L]
  fit_grid <- rbind(
    expand.grid(
      config = configs,
      endpoint_visit = visits,
      analysis_mode = "transition",
      stringsAsFactors = FALSE
    ),
    expand.grid(
      config = configs,
      endpoint_visit = visits,
      analysis_mode = "cumulative",
      stringsAsFactors = FALSE
    )
  )

  fits <- vector("list", nrow(fit_grid))
  status <- vector("list", nrow(fit_grid))
  for (i in seq_len(nrow(fit_grid))) {
    spec <- fit_grid[i, ]
    config <- spec$config
    endpoint <- spec$endpoint_visit
    mode <- spec$analysis_mode
    n_sources <- length(get_true_source_groups(
      simulated_mod_dat, block, treat, config
    ))

    started <- proc.time()[["elapsed"]]
    result <- tryCatch(
      fit_source_detection_experiment(
        block = block,
        treat = treat,
        config = config,
        n_sources = n_sources,
        mod_dat = simulated_mod_dat,
        theta_initial = forward_fit$theta,
        through_visit = if (mode == "cumulative") endpoint else NULL,
        transition_visit = if (mode == "transition") endpoint else NULL,
        include_naive = TRUE
      ),
      error = function(e) e
    )
    elapsed <- proc.time()[["elapsed"]] - started
    success <- !inherits(result, "error")
    status[[i]] <- data.frame(
      simulation = simulation,
      block = block,
      treat = as.character(treat),
      config = config,
      n_sources = n_sources,
      endpoint_visit = endpoint,
      analysis_mode = mode,
      success = success,
      elapsed_seconds = elapsed,
      error = if (success) NA_character_ else conditionMessage(result)
    )
    if (success) fits[[i]] <- result
  }

  successful <- Filter(Negate(is.null), fits)
  summaries <- lapply(successful, summarize_backward_fit, simulation = simulation)
  list(
    fit_status = dplyr::bind_rows(status),
    diagnostics = dplyr::bind_rows(lapply(summaries, `[[`, "diagnostics")),
    accuracy = dplyr::bind_rows(lapply(summaries, `[[`, "accuracy")),
    parameters = dplyr::bind_rows(lapply(summaries, `[[`, "parameters")),
    fits = fits
  )
}

run_simulation_scenario <- function(
  simulation,
  block,
  treat,
  mod_dat,
  generating_fit,
  kappa_grid = exp(seq(log(0.5), log(2), length.out = 5)),
  save_fits = FALSE
) {
  scenario_seed <- 100000L * as.integer(simulation) +
    1000L * match(block, dimnames(mod_dat$intensity)[["block"]]) +
    10L * match(as.character(treat), dimnames(mod_dat$intensity)[["treat"]])
  simulated <- simulate_whole_epidemic(
    mod_dat = mod_dat,
    forward_fit = generating_fit,
    block = block,
    treat = treat,
    mechanism = "secondary_spread",
    seed = scenario_seed
  )
  simulated_mod_dat <- replace_experiment_intensity(
    mod_dat, block, treat, simulated
  )

  forward_started <- proc.time()[["elapsed"]]
  forward_result <- tryCatch(
    fit_forward_experiment(
      block, treat, simulated_mod_dat, kappa_grid = kappa_grid
    ),
    error = function(e) e
  )
  forward_elapsed <- proc.time()[["elapsed"]] - forward_started
  if (inherits(forward_result, "error")) {
    return(list(
      scenario_status = data.frame(
        simulation = simulation, block = block, treat = as.character(treat),
        forward_success = FALSE, forward_elapsed_seconds = forward_elapsed,
        error = conditionMessage(forward_result)
      ),
      epidemic_summary = epidemic_summary(
        simulated, block, treat, simulation, "secondary_spread", "simulated"
      ),
      progression_diagnostics = progression_diagnostics(
        simulated, block, treat, simulation
      )
    ))
  }

  configs <- simulation_configs(
    treat, dimnames(simulated_mod_dat$groups)[["config"]]
  )
  backward <- run_backward_simulation_fits(
    simulated_mod_dat, forward_result, block, treat, simulation, configs
  )
  result <- list(
    scenario_status = data.frame(
      simulation = simulation, block = block, treat = as.character(treat),
      forward_success = TRUE, forward_elapsed_seconds = forward_elapsed,
      error = NA_character_
    ),
    epidemic_summary = epidemic_summary(
      simulated, block, treat, simulation, "secondary_spread", "simulated"
    ),
    progression_diagnostics = progression_diagnostics(
      simulated, block, treat, simulation
    ),
    forward_metrics = summarize_forward_simulation_fit(
      forward_result, generating_fit, block, treat, simulation
    ),
    backward_fit_status = backward$fit_status,
    backward_diagnostics = backward$diagnostics,
    backward_accuracy = backward$accuracy,
    backward_parameters = backward$parameters
  )
  if (save_fits) {
    result$simulated_intensity <- simulated
    result$forward_fit <- forward_result
    result$backward_fits <- backward$fits
  }
  result
}

combine_simulation_results <- function(results) {
  component <- function(name) {
    dplyr::bind_rows(lapply(results, function(x) x[[name]]))
  }
  list(
    scenario_status = component("scenario_status"),
    epidemic_summary = component("epidemic_summary"),
    progression_diagnostics = component("progression_diagnostics"),
    forward_metrics = component("forward_metrics"),
    backward_fit_status = component("backward_fit_status"),
    backward_diagnostics = component("backward_diagnostics"),
    backward_accuracy = component("backward_accuracy"),
    backward_parameters = component("backward_parameters")
  )
}
