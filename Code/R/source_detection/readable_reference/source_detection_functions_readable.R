# Readable reference implementation of the source-detection model.
#
# This file mirrors the production analysis but favors explicit calculations,
# loops, and named intermediate objects over speed. It intentionally uses the
# finite-difference gradient supplied by optim() so it can serve as an
# independent check on the optimized analytic-gradient implementation.

log_sum_exp_readable <- function(values) {
  largest_value <- max(values)
  largest_value + log(sum(exp(values - largest_value)))
}

make_candidate_sets_readable <- function(group_id, number_of_sources) {
  possible_groups <- sort(unique(as.numeric(group_id)))
  candidate_matrix <- combn(possible_groups, number_of_sources)
  candidate_names <- character(ncol(candidate_matrix))

  for (candidate in seq_len(ncol(candidate_matrix))) {
    candidate_names[candidate] <- paste(
      candidate_matrix[, candidate], collapse = "_"
    )
  }
  colnames(candidate_matrix) <- candidate_names

  list(matrix = candidate_matrix, names = candidate_names)
}

get_true_sources_readable <- function(model_data, block, treatment, resolution) {
  recorded_truth <- model_data$truth[block, treatment, , resolution]
  sort(unique(as.numeric(na.omit(recorded_truth))))
}

extract_study_window_readable <- function(
  model_data,
  block,
  treatment,
  through_visit = NULL,
  single_transition = NULL
) {
  study_visits <- dimnames(model_data$intensity)[["visit"]]

  if (!is.null(single_transition)) {
    transition_position <- match(as.character(single_transition), study_visits)
    included_visits <- study_visits[c(
      transition_position - 1L,
      transition_position
    )]
  } else {
    if (is.null(through_visit)) through_visit <- tail(study_visits, 1L)
    endpoint_position <- match(as.character(through_visit), study_visits)
    included_visits <- study_visits[seq_len(endpoint_position)]
  }

  transition_visits <- included_visits[-1L]
  severity_by_visit <- vector("list", length(included_visits))
  wind_by_transition <- vector("list", length(transition_visits))

  for (visit_index in seq_along(included_visits)) {
    visit <- included_visits[visit_index]
    severity_by_visit[[visit_index]] <-
      model_data$intensity[, block, treatment, visit]
  }

  for (transition_index in seq_along(transition_visits)) {
    visit <- transition_visits[transition_index]
    wind_by_transition[[transition_index]] <-
      model_data$wind[, , block, treatment, visit]
  }

  list(
    included_visits = included_visits,
    transition_visits = transition_visits,
    through_visit = tail(included_visits, 1L),
    final_study_visit = tail(study_visits, 1L),
    severity = severity_by_visit,
    wind = wind_by_transition,
    distance = model_data$dist
  )
}

candidate_dispersal_readable <- function(
  previous_severity,
  wind_matrix,
  distance_matrix,
  group_id,
  candidate_sources,
  kappa,
  distance_offset = 0.01
) {
  number_of_plants <- length(previous_severity)
  dispersal <- numeric(number_of_plants)

  for (target in seq_len(number_of_plants)) {
    target_dispersal <- 0

    for (source in seq_len(number_of_plants)) {
      source_is_allowed <- group_id[source] %in% candidate_sources
      source_is_target <- source == target

      if (source_is_allowed && !source_is_target) {
        distance_weight <-
          (distance_matrix[target, source] + distance_offset)^(-kappa)
        source_contribution <-
          previous_severity[source] *
          wind_matrix[target, source] *
          distance_weight
        target_dispersal <- target_dispersal + source_contribution
      }
    }

    dispersal[target] <- target_dispersal
  }

  dispersal
}

zib_loglikelihood_readable <- function(
  theta,
  current_severity,
  previous_severity,
  wind_matrix,
  distance_matrix,
  group_id,
  candidate_sources,
  distance_offset = 0.01
) {
  beta <- theta[["beta"]]
  delta <- theta[["delta"]]
  gamma <- theta[["gamma"]]
  kappa <- exp(theta[["log_kappa"]])
  phi <- exp(theta[["log_phi"]])

  dispersal <- candidate_dispersal_readable(
    previous_severity = previous_severity,
    wind_matrix = wind_matrix,
    distance_matrix = distance_matrix,
    group_id = group_id,
    candidate_sources = candidate_sources,
    kappa = kappa,
    distance_offset = distance_offset
  )

  autoinfection <- previous_severity * (1 - previous_severity)
  linear_predictor <- beta + delta * autoinfection + gamma * dispersal
  conditional_mean <- plogis(linear_predictor)
  conditional_mean <- pmin(pmax(conditional_mean, 1e-8), 1 - 1e-8)
  zero_probability <- mean(current_severity == 0)

  contributions <- numeric(length(current_severity))
  for (plant in seq_along(current_severity)) {
    observed_severity <- current_severity[plant]

    if (observed_severity == 0) {
      contributions[plant] <- log(zero_probability)
    } else {
      beta_log_density <- dbeta(
        observed_severity,
        shape1 = conditional_mean[plant] * phi,
        shape2 = (1 - conditional_mean[plant]) * phi,
        log = TRUE
      )
      contributions[plant] <-
        log1p(-zero_probability) + beta_log_density
    }
  }

  sum(contributions)
}

negative_zib_loglikelihood_readable <- function(...) {
  -zib_loglikelihood_readable(...)
}

all_candidate_loglikelihoods_readable <- function(
  theta_by_transition,
  experiment,
  group_id,
  candidate_matrix,
  distance_offset = 0.01
) {
  number_of_candidates <- ncol(candidate_matrix)
  number_of_transitions <- length(experiment$transition_visits)
  loglikelihoods <- matrix(
    NA_real_,
    nrow = number_of_candidates,
    ncol = number_of_transitions,
    dimnames = list(
      candidate = colnames(candidate_matrix),
      transition = experiment$transition_visits
    )
  )

  for (candidate in seq_len(number_of_candidates)) {
    candidate_sources <- candidate_matrix[, candidate]

    for (transition in seq_len(number_of_transitions)) {
      loglikelihoods[candidate, transition] <- zib_loglikelihood_readable(
        theta = theta_by_transition[transition, ],
        current_severity = experiment$severity[[transition + 1L]],
        previous_severity = experiment$severity[[transition]],
        wind_matrix = experiment$wind[[transition]],
        distance_matrix = experiment$distance,
        group_id = group_id,
        candidate_sources = candidate_sources,
        distance_offset = distance_offset
      )
    }
  }

  loglikelihoods
}

candidate_responsibilities_readable <- function(
  candidate_loglikelihood,
  log_prior
) {
  log_weights <- log_prior + candidate_loglikelihood
  exp(log_weights - log_sum_exp_readable(log_weights))
}

negative_q_for_transition_readable <- function(
  theta,
  transition,
  experiment,
  group_id,
  candidate_matrix,
  responsibilities,
  distance_offset = 0.01
) {
  expected_loglikelihood <- 0

  for (candidate in seq_len(ncol(candidate_matrix))) {
    candidate_loglikelihood <- zib_loglikelihood_readable(
      theta = theta,
      current_severity = experiment$severity[[transition + 1L]],
      previous_severity = experiment$severity[[transition]],
      wind_matrix = experiment$wind[[transition]],
      distance_matrix = experiment$distance,
      group_id = group_id,
      candidate_sources = candidate_matrix[, candidate],
      distance_offset = distance_offset
    )

    expected_loglikelihood <- expected_loglikelihood +
      responsibilities[candidate] * candidate_loglikelihood
  }

  -expected_loglikelihood
}

initialize_transition_readable <- function(
  current_severity,
  previous_severity,
  wind_matrix,
  distance_matrix,
  initial_kappa,
  distance_offset = 0.01
) {
  full_forward_dispersal <- candidate_dispersal_readable(
    previous_severity = previous_severity,
    wind_matrix = wind_matrix,
    distance_matrix = distance_matrix,
    group_id = seq_along(previous_severity),
    candidate_sources = seq_along(previous_severity),
    kappa = initial_kappa,
    distance_offset = distance_offset
  )
  positive_plants <- which(current_severity > 0)
  transformed_severity <- qlogis(pmin(
    pmax(current_severity[positive_plants], 1e-5),
    1 - 1e-5
  ))
  autoinfection <-
    previous_severity[positive_plants] *
    (1 - previous_severity[positive_plants])

  moment_regression <- lm(
    transformed_severity ~
      autoinfection + full_forward_dispersal[positive_plants]
  )
  regression_coefficients <- coef(moment_regression)
  fitted_conditional_mean <- plogis(fitted(moment_regression))
  response_scale_residual_variance <- mean(
    residuals(moment_regression)^2 *
      (fitted_conditional_mean * (1 - fitted_conditional_mean))^2
  )
  initial_phi <-
    mean(fitted_conditional_mean * (1 - fitted_conditional_mean)) /
    max(response_scale_residual_variance, 1e-6) - 1
  initial_phi <- max(initial_phi, 1e-3)

  c(
    beta = unname(regression_coefficients[1]),
    delta = unname(regression_coefficients[2]),
    gamma = unname(regression_coefficients[3]),
    log_kappa = log(initial_kappa),
    log_phi = log(initial_phi)
  )
}

fit_forward_model_readable <- function(
  block,
  treatment,
  model_data,
  initial_kappa_grid = exp(seq(log(0.5), log(2), length.out = 5)),
  distance_offset = 0.01
) {
  experiment <- extract_study_window_readable(model_data, block, treatment)
  number_of_transitions <- length(experiment$transition_visits)
  theta_by_transition <- matrix(
    NA_real_,
    nrow = number_of_transitions,
    ncol = 5,
    dimnames = list(
      experiment$transition_visits,
      c("beta", "delta", "gamma", "log_kappa", "log_phi")
    )
  )

  for (transition in seq_len(number_of_transitions)) {
    fits_for_transition <- vector("list", length(initial_kappa_grid))

    for (start in seq_along(initial_kappa_grid)) {
      initial_kappa <- initial_kappa_grid[start]
      initial_theta <- initialize_transition_readable(
        current_severity = experiment$severity[[transition + 1L]],
        previous_severity = experiment$severity[[transition]],
        wind_matrix = experiment$wind[[transition]],
        distance_matrix = experiment$distance,
        initial_kappa = initial_kappa,
        distance_offset = distance_offset
      )

      every_plant_is_a_source <- seq_len(nrow(experiment$distance))
      fits_for_transition[[start]] <- optim(
        par = initial_theta,
        fn = negative_zib_loglikelihood_readable,
        method = "BFGS",
        control = list(maxit = 5000, reltol = 1e-8),
        current_severity = experiment$severity[[transition + 1L]],
        previous_severity = experiment$severity[[transition]],
        wind_matrix = experiment$wind[[transition]],
        distance_matrix = experiment$distance,
        group_id = seq_len(nrow(experiment$distance)),
        candidate_sources = every_plant_is_a_source,
        distance_offset = distance_offset
      )
    }

    optimized_loglikelihoods <- numeric(length(fits_for_transition))
    for (start in seq_along(fits_for_transition)) {
      optimized_loglikelihoods[start] <- -fits_for_transition[[start]]$value
    }
    best_start <- which.max(optimized_loglikelihoods)
    theta_by_transition[transition, ] <-
      fits_for_transition[[best_start]]$par
  }

  theta_by_transition
}

fit_source_detection_readable <- function(
  block,
  treatment,
  resolution,
  number_of_sources,
  model_data,
  initial_theta,
  through_visit = NULL,
  single_transition = NULL,
  maximum_em_iterations = 200,
  loglikelihood_tolerance = 1e-3,
  responsibility_tolerance = 1e-3,
  distance_offset = 0.01
) {
  experiment <- extract_study_window_readable(
    model_data = model_data,
    block = block,
    treatment = treatment,
    through_visit = through_visit,
    single_transition = single_transition
  )
  group_id <- model_data$groups[, resolution]
  candidates <- make_candidate_sets_readable(
    group_id,
    number_of_sources
  )
  number_of_candidates <- ncol(candidates$matrix)
  log_prior <- rep(-log(number_of_candidates), number_of_candidates)
  theta <- initial_theta[experiment$transition_visits, , drop = FALSE]

  transition_loglikelihood <- all_candidate_loglikelihoods_readable(
    theta,
    experiment,
    group_id,
    candidates$matrix,
    distance_offset
  )
  candidate_loglikelihood <- rowSums(transition_loglikelihood)
  responsibilities <- candidate_responsibilities_readable(
    candidate_loglikelihood,
    log_prior
  )
  observed_loglikelihood <- log_sum_exp_readable(
    log_prior + candidate_loglikelihood
  )

  history <- data.frame(
    iteration = 0,
    observed_loglikelihood = observed_loglikelihood,
    maximum_responsibility_change = NA_real_
  )

  for (em_iteration in seq_len(maximum_em_iterations)) {
    new_theta <- theta

    for (transition in seq_along(experiment$transition_visits)) {
      transition_fit <- optim(
        par = theta[transition, ],
        fn = negative_q_for_transition_readable,
        method = "BFGS",
        control = list(maxit = 1000, reltol = 1e-8),
        transition = transition,
        experiment = experiment,
        group_id = group_id,
        candidate_matrix = candidates$matrix,
        responsibilities = responsibilities,
        distance_offset = distance_offset
      )
      new_theta[transition, ] <- transition_fit$par
    }

    new_transition_loglikelihood <- all_candidate_loglikelihoods_readable(
      new_theta,
      experiment,
      group_id,
      candidates$matrix,
      distance_offset
    )
    new_candidate_loglikelihood <- rowSums(new_transition_loglikelihood)
    new_responsibilities <- candidate_responsibilities_readable(
      new_candidate_loglikelihood,
      log_prior
    )
    new_observed_loglikelihood <- log_sum_exp_readable(
      log_prior + new_candidate_loglikelihood
    )

    loglikelihood_change <-
      new_observed_loglikelihood - observed_loglikelihood
    maximum_responsibility_change <- max(abs(
      new_responsibilities - responsibilities
    ))

    history <- rbind(history, data.frame(
      iteration = em_iteration,
      observed_loglikelihood = new_observed_loglikelihood,
      maximum_responsibility_change = maximum_responsibility_change
    ))

    theta <- new_theta
    transition_loglikelihood <- new_transition_loglikelihood
    candidate_loglikelihood <- new_candidate_loglikelihood
    responsibilities <- new_responsibilities
    observed_loglikelihood <- new_observed_loglikelihood

    if (
      abs(loglikelihood_change) <= loglikelihood_tolerance &&
      maximum_responsibility_change <= responsibility_tolerance
    ) {
      break
    }
  }

  fit_scope <- if (!is.null(single_transition)) {
    "transition_refit"
  } else if (experiment$through_visit == experiment$final_study_visit) {
    "whole_study"
  } else {
    "cumulative_refit"
  }
  true_sources <- get_true_sources_readable(
    model_data,
    block,
    treatment,
    resolution
  )
  winning_candidate <- which.max(responsibilities)
  predicted_sources <- candidates$matrix[, winning_candidate]
  naive_sources <- naive_prediction_readable(
    experiment,
    group_id,
    number_of_sources
  )
  model_accuracy <- source_accuracy_readable(
    predicted_sources,
    true_sources,
    model_data$grid_dist[[resolution]]
  )
  naive_accuracy <- source_accuracy_readable(
    naive_sources,
    true_sources,
    model_data$grid_dist[[resolution]]
  )
  candidate_summary <- data.frame(
    block = block,
    treat = treatment,
    config = resolution,
    n_sources = number_of_sources,
    true_n_sources = length(true_sources),
    through_visit = experiment$through_visit,
    fit_scope = fit_scope,
    candidate_id = candidates$names,
    loglik = candidate_loglikelihood,
    posterior = responsibilities,
    rank = rank(-responsibilities, ties.method = "min"),
    predicted = seq_len(number_of_candidates) == winning_candidate,
    true_candidate_id = paste(true_sources, collapse = "_"),
    true_candidate = candidates$names == paste(true_sources, collapse = "_")
  )

  parameter_summary <- data.frame(
    block = block,
    treat = treatment,
    config = resolution,
    n_sources = number_of_sources,
    true_n_sources = length(true_sources),
    through_visit = experiment$through_visit,
    fit_scope = fit_scope,
    visit = experiment$transition_visits,
    beta = theta[, "beta"],
    delta = theta[, "delta"],
    gamma = theta[, "gamma"],
    kappa = exp(theta[, "log_kappa"]),
    phi = exp(theta[, "log_phi"]),
    observed_loglik = observed_loglikelihood
  )

  accuracy_summary <- data.frame(
    block = block,
    treat = treatment,
    config = resolution,
    n_sources = number_of_sources,
    true_n_sources = length(true_sources),
    through_visit = experiment$through_visit,
    fit_scope = fit_scope,
    visit = experiment$through_visit,
    predicted_candidate_id = paste(predicted_sources, collapse = "_"),
    naive_candidate_id = paste(naive_sources, collapse = "_"),
    source_count_correct = number_of_sources == length(true_sources),
    standard_accuracy = model_accuracy$standard_accuracy,
    precision = model_accuracy$precision,
    recall = model_accuracy$recall,
    f1 = model_accuracy$f1,
    distance_weighted_accuracy =
      model_accuracy$distance_weighted_accuracy,
    naive_standard_accuracy = naive_accuracy$standard_accuracy,
    naive_precision = naive_accuracy$precision,
    naive_recall = naive_accuracy$recall,
    naive_f1 = naive_accuracy$f1,
    naive_distance_weighted_accuracy =
      naive_accuracy$distance_weighted_accuracy
  )

  list(
    candidate_summary = candidate_summary,
    parameter_summary = parameter_summary,
    accuracy_summary = accuracy_summary,
    history = history,
    theta = theta,
    transition_loglikelihood = transition_loglikelihood
  )
}

naive_prediction_readable <- function(
  experiment,
  group_id,
  number_of_sources
) {
  endpoint_severity <- tail(experiment$severity, 1L)[[1]]
  group_names <- sort(unique(as.numeric(group_id)))
  maximum_severity <- numeric(length(group_names))

  for (group in seq_along(group_names)) {
    plants_in_group <- which(group_id == group_names[group])
    maximum_severity[group] <- max(endpoint_severity[plants_in_group])
  }

  ranking <- order(-maximum_severity, group_names)
  sort(group_names[ranking[seq_len(number_of_sources)]])
}

source_accuracy_readable <- function(
  predicted_sources,
  true_sources,
  distance_matrix
) {
  predicted_sources <- as.numeric(predicted_sources)
  true_sources <- as.numeric(true_sources)
  number_predicted <- length(predicted_sources)
  number_true <- length(true_sources)
  assignment_size <- max(number_predicted, number_true)
  unmatched_cost <- max(distance_matrix)

  assignment_cost <- matrix(
    unmatched_cost,
    nrow = assignment_size,
    ncol = assignment_size
  )
  assignment_cost[seq_len(number_predicted), seq_len(number_true)] <-
    distance_matrix[predicted_sources, true_sources, drop = FALSE]
  assignment <- RcppHungarian::HungarianSolver(assignment_cost)$pairs

  spatial_accuracy <- numeric(assignment_size)
  for (pair in seq_len(nrow(assignment))) {
    predicted_index <- assignment[pair, 1]
    true_index <- assignment[pair, 2]

    if (predicted_index <= number_predicted && true_index <= number_true) {
      predicted_group <- predicted_sources[predicted_index]
      true_group <- true_sources[true_index]
      prediction_distance <- distance_matrix[predicted_group, true_group]
      maximum_distance <- max(distance_matrix[, true_group])
      spatial_accuracy[pair] <-
        1 - prediction_distance / maximum_distance
    } else {
      spatial_accuracy[pair] <- 0
    }
  }

  number_correct <- sum(predicted_sources %in% true_sources)
  precision <- number_correct / number_predicted
  recall <- number_correct / number_true
  f1 <- if (precision + recall == 0) {
    0
  } else {
    2 * precision * recall / (precision + recall)
  }

  list(
    standard_accuracy = number_correct / assignment_size,
    precision = precision,
    recall = recall,
    f1 = f1,
    distance_weighted_accuracy = mean(spatial_accuracy)
  )
}
