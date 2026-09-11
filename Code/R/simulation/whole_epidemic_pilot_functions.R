# Functions for a small whole-epidemic simulation pilot.
#
# These functions deliberately reuse the production transition equation but
# keep simulation and diagnostic calculations separate from the main analysis.

simulate_zib_transition <- function(mu, phi, alpha) {
  # This is an exact draw from the fitted conditional ZIB likelihood. The
  # likelihood does not impose monotone severity, so simulated trajectories
  # may decrease or return to zero; those events are recorded as diagnostics.
  is_zero <- rbinom(length(mu), size = 1L, prob = alpha) == 1L
  y <- numeric(length(mu))
  positive <- !is_zero
  y[positive] <- rbeta(
    sum(positive),
    shape1 = mu[positive] * phi,
    shape2 = (1 - mu[positive]) * phi
  )
  y
}

simulate_whole_epidemic <- function(
  mod_dat,
  forward_fit,
  block,
  treat,
  mechanism = c("secondary_spread", "persistent_sources"),
  seed = 1L,
  d0 = 0.01
) {
  mechanism <- match.arg(mechanism)
  visits <- dimnames(mod_dat$intensity)[["visit"]]
  transition_visits <- visits[-1L]
  n_plants <- dim(mod_dat$intensity)[1L]
  simulated <- matrix(
    0,
    nrow = n_plants,
    ncol = length(visits),
    dimnames = list(plant = dimnames(mod_dat$intensity)[["plant"]], visit = visits)
  )

  # The observed first visit supplies the experimentally imposed initial state.
  simulated[, 1L] <- mod_dat$intensity[, block, treat, visits[1L]]
  source_plants <- get_true_source_groups(mod_dat, block, treat, "64")
  source_indicator <- as.numeric(seq_len(n_plants) %in% source_plants)
  observed_alpha <- vapply(transition_visits, function(visit) {
    mean(mod_dat$intensity[, block, treat, visit] == 0)
  }, numeric(1))

  set.seed(seed)
  for (t in seq_along(transition_visits)) {
    visit <- transition_visits[t]
    transmitter_indicator <- if (mechanism == "secondary_spread") {
      rep(1, n_plants)
    } else {
      source_indicator
    }
    components <- transition_components(
      theta = forward_fit$theta[visit, ],
      y_previous = simulated[, t],
      wind_matrix = mod_dat$wind[, , block, treat, visit],
      distance_matrix = mod_dat$dist,
      candidate_indicator = matrix(transmitter_indicator, ncol = 1L),
      d0 = d0
    )
    simulated[, t + 1L] <- simulate_zib_transition(
      mu = components$mu[, 1L],
      phi = exp(forward_fit$theta[visit, "log_phi"]),
      alpha = observed_alpha[t]
    )
  }

  simulated
}

progression_diagnostics <- function(y, block, treat, simulation) {
  do.call(rbind, lapply(seq.int(2L, ncol(y)), function(v) {
    previous <- y[, v - 1L]
    current <- y[, v]
    change <- current - previous
    decreased <- change < 0
    data.frame(
      simulation = simulation,
      block = block,
      treat = as.character(treat),
      visit = colnames(y)[v],
      n_decreased = sum(decreased),
      proportion_decreased = mean(decreased),
      mean_decrease = if (any(decreased)) mean(-change[decreased]) else 0,
      maximum_decrease = if (any(decreased)) max(-change[decreased]) else 0,
      n_positive_to_zero = sum(previous > 0 & current == 0),
      proportion_positive_to_zero = mean(previous > 0 & current == 0),
      any_decrease = any(decreased),
      any_positive_to_zero = any(previous > 0 & current == 0)
    )
  }))
}

replace_experiment_intensity <- function(mod_dat, block, treat, simulated) {
  simulated_data <- mod_dat
  simulated_data$intensity[, block, treat, ] <- simulated
  simulated_data
}

epidemic_summary <- function(y, block, treat, replicate, mechanism, data_type) {
  do.call(rbind, lapply(seq_len(ncol(y)), function(v) {
    values <- y[, v]
    data.frame(
      block = block,
      treat = treat,
      replicate = replicate,
      mechanism = mechanism,
      data_type = data_type,
      visit = colnames(y)[v],
      zero_fraction = mean(values == 0),
      mean_severity = mean(values),
      median_severity = median(values),
      q90_severity = unname(quantile(values, 0.9)),
      maximum_severity = max(values)
    )
  }))
}

plant_level_summary <- function(
  y, coordinates, source_plants, block, treat, replicate, mechanism, data_type
) {
  result <- expand.grid(
    plant = seq_len(nrow(y)),
    visit = colnames(y),
    stringsAsFactors = FALSE
  )
  result$severity <- as.vector(y)
  result$east <- coordinates$east[result$plant]
  result$north <- coordinates$north[result$plant]
  result$true_source <- result$plant %in% source_plants
  result$block <- block
  result$treat <- treat
  result$replicate <- replicate
  result$mechanism <- mechanism
  result$data_type <- data_type
  result
}

oracle_source_diagnostics <- function(
  simulated_mod_dat,
  forward_fit,
  block,
  treat,
  config,
  replicate,
  mechanism,
  d0 = 0.01
) {
  experiment <- extract_experiment(simulated_mod_dat, block, treat)
  true_sources <- get_true_source_groups(simulated_mod_dat, block, treat, config)
  candidates <- make_candidate_sets(
    simulated_mod_dat$groups[, config], length(true_sources)
  )
  loglik <- candidate_loglik_matrix(
    forward_fit$theta,
    experiment,
    candidates$indicator,
    d0
  )
  log_prior <- rep(-log(length(candidates$ids)), length(candidates$ids))
  true_id <- paste(true_sources, collapse = "_")

  summaries <- make_posterior_summaries(loglik, log_prior, candidates$ids)
  do.call(rbind, lapply(split(summaries, interaction(
    summaries$evidence, summaries$visit, drop = TRUE
  )), function(x) {
    true_row <- x[x$candidate_id == true_id, , drop = FALSE]
    winner <- x[which.max(x$posterior), , drop = FALSE]
    data.frame(
      block = block,
      treat = treat,
      replicate = replicate,
      mechanism = mechanism,
      config = config,
      n_sources = length(true_sources),
      evidence = winner$evidence,
      visit = winner$visit,
      true_candidate_id = true_id,
      predicted_candidate_id = winner$candidate_id,
      correct = winner$candidate_id == true_id,
      true_rank = true_row$rank,
      true_posterior = true_row$posterior,
      winner_posterior = winner$posterior,
      effective_candidates = exp(-sum(ifelse(
        x$posterior > 0, x$posterior * log(x$posterior), 0
      )))
    )
  }))
}
