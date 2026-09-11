# Functions for the main forward-backward source-detection analysis.
#
# A source set is fixed for an entire block x treatment experiment. Model
# parameters are shared by all candidate sets but may differ by transition.
# beta, delta, and gamma are unconstrained; kappa and phi are optimized on
# their log scales.

log_sum_exp <- function(x) {
  x_max <- max(x)
  x_max + log(sum(exp(x - x_max)))
}

inv_logit <- function(x) plogis(x)

logit <- function(p) qlogis(p)

make_candidate_sets <- function(group_id, n_sources) {
  group_levels <- sort(unique(as.numeric(group_id)))
  if (n_sources < 1L || n_sources > length(group_levels)) {
    stop("n_sources must be between 1 and the number of spatial groups.")
  }

  candidate_matrix <- combn(group_levels, n_sources)
  candidate_ids <- apply(candidate_matrix, 2, paste, collapse = "_")
  colnames(candidate_matrix) <- candidate_ids

  indicator <- vapply(seq_len(ncol(candidate_matrix)), function(k) {
    as.numeric(group_id %in% candidate_matrix[, k])
  }, numeric(length(group_id)))
  if (is.null(dim(indicator))) indicator <- matrix(indicator, ncol = 1L)
  colnames(indicator) <- candidate_ids

  list(
    matrix = candidate_matrix,
    ids = candidate_ids,
    indicator = indicator
  )
}

get_true_source_groups <- function(mod_dat, block, treat, config) {
  sort(unique(na.omit(as.numeric(mod_dat$truth[block, treat, , config]))))
}

# Pair predicted and true source cells by the minimum total spatial distance.
source_accuracy <- function(predicted_sources, true_sources, distance_matrix) {
  predicted_sources <- as.numeric(predicted_sources)
  true_sources <- as.numeric(true_sources)
  if (length(predicted_sources) != length(true_sources)) {
    stop("Predicted and true source sets must have the same length.")
  }

  error_matrix <- distance_matrix[predicted_sources, true_sources, drop = FALSE]
  assignment <- RcppHungarian::HungarianSolver(error_matrix)$pairs
  assigned_predictions <- predicted_sources[assignment[, 1]]
  assigned_truth <- true_sources[assignment[, 2]]
  assigned_distance <- error_matrix[assignment]

  # The denominator is specific to each true source: its distance to the
  # furthest candidate-cell centroid in the study area.
  maximum_distance <- vapply(assigned_truth, function(true_source) {
    max(distance_matrix[, true_source])
  }, numeric(1))
  individual_spatial_accuracy <- ifelse(
    maximum_distance > 0,
    1 - assigned_distance / maximum_distance,
    as.numeric(assigned_distance == 0)
  )

  list(
    standard_accuracy = mean(predicted_sources %in% true_sources),
    n_correct = sum(predicted_sources %in% true_sources),
    mean_distance = mean(assigned_distance),
    distance_weighted_accuracy = mean(individual_spatial_accuracy),
    assignment = data.frame(
      true_source = assigned_truth,
      predicted_source = assigned_predictions,
      distance = assigned_distance,
      maximum_distance = maximum_distance,
      spatial_accuracy = individual_spatial_accuracy
    )
  )
}

# Accuracy when the predicted and true source sets may have different sizes.
# Real source groups are paired by minimum total Euclidean distance after the
# distance matrix is padded with dummy sources. An unmatched predicted or true
# source receives zero spatial accuracy. For equal-sized sets this reduces to
# the original Hungarian assignment and accuracy definitions above.
source_accuracy_general <- function(
  predicted_sources,
  true_sources,
  distance_matrix
) {
  predicted_sources <- as.numeric(predicted_sources)
  true_sources <- as.numeric(true_sources)
  n_predicted <- length(predicted_sources)
  n_true <- length(true_sources)
  assignment_size <- max(n_predicted, n_true)
  maximum_study_distance <- max(distance_matrix)

  padded_error <- matrix(
    maximum_study_distance,
    nrow = assignment_size,
    ncol = assignment_size
  )
  padded_error[seq_len(n_predicted), seq_len(n_true)] <-
    distance_matrix[predicted_sources, true_sources, drop = FALSE]
  if (n_predicted < assignment_size && n_true < assignment_size) {
    padded_error[
      (n_predicted + 1L):assignment_size,
      (n_true + 1L):assignment_size
    ] <- 0
  }

  assignment <- RcppHungarian::HungarianSolver(padded_error)$pairs
  assignment_details <- do.call(rbind, lapply(seq_len(nrow(assignment)), function(i) {
    predicted_index <- assignment[i, 1]
    true_index <- assignment[i, 2]
    matched_prediction <- predicted_index <= n_predicted
    matched_truth <- true_index <= n_true
    predicted_source <- if (matched_prediction) {
      predicted_sources[predicted_index]
    } else {
      NA_real_
    }
    true_source <- if (matched_truth) true_sources[true_index] else NA_real_

    if (matched_prediction && matched_truth) {
      distance <- distance_matrix[predicted_source, true_source]
      maximum_distance <- max(distance_matrix[, true_source])
      spatial_accuracy <- if (maximum_distance > 0) {
        1 - distance / maximum_distance
      } else {
        as.numeric(distance == 0)
      }
    } else {
      distance <- NA_real_
      maximum_distance <- if (matched_truth) {
        max(distance_matrix[, true_source])
      } else {
        NA_real_
      }
      spatial_accuracy <- 0
    }

    data.frame(
      true_source = true_source,
      predicted_source = predicted_source,
      distance = distance,
      maximum_distance = maximum_distance,
      spatial_accuracy = spatial_accuracy
    )
  }))

  n_correct <- sum(predicted_sources %in% true_sources)
  precision <- n_correct / n_predicted
  recall <- n_correct / n_true
  f1 <- if (precision + recall > 0) {
    2 * precision * recall / (precision + recall)
  } else {
    0
  }

  list(
    source_count_correct = n_predicted == n_true,
    n_predicted = n_predicted,
    n_true = n_true,
    n_correct = n_correct,
    # This symmetric extension equals the previous exact accuracy when the
    # source-set sizes agree and penalizes missing or extra sources otherwise.
    standard_accuracy = n_correct / assignment_size,
    precision = precision,
    recall = recall,
    f1 = f1,
    mean_distance = mean(
      assignment_details$distance,
      na.rm = TRUE
    ),
    distance_weighted_accuracy = mean(
      assignment_details$spatial_accuracy
    ),
    assignment = assignment_details
  )
}

# Select the groups containing the greatest observed severity at the endpoint
# of a fitted window. Group ID resolves exact ties deterministically.
naive_source_prediction <- function(experiment, group_id, n_sources) {
  endpoint_severity <- tail(experiment$y, 1L)[[1]]
  group_severity <- tapply(endpoint_severity, group_id, max, na.rm = TRUE)
  group_ids <- as.numeric(names(group_severity))
  ordering <- order(-as.numeric(group_severity), group_ids)
  sort(group_ids[ordering[seq_len(n_sources)]])
}

dispersal_by_candidate <- function(
  y_previous,
  wind_matrix,
  distance_matrix,
  candidate_indicator,
  kappa,
  d0 = 0.01,
  derivative = FALSE
) {
  shifted_distance <- distance_matrix + d0
  kernel <- shifted_distance^(-kappa)
  if (derivative) kernel <- -kernel * log(shifted_distance)

  transport <- wind_matrix * kernel
  diag(transport) <- 0
  transport %*% (y_previous * candidate_indicator)
}

transition_components <- function(
  theta,
  y_previous,
  wind_matrix,
  distance_matrix,
  candidate_indicator,
  d0 = 0.01,
  include_kappa_derivative = FALSE
) {
  kappa <- exp(theta[["log_kappa"]])
  dispersal <- dispersal_by_candidate(
    y_previous,
    wind_matrix,
    distance_matrix,
    candidate_indicator,
    kappa,
    d0
  )
  auto <- y_previous * (1 - y_previous)
  eta <- theta[["beta"]] +
    theta[["delta"]] * auto +
    theta[["gamma"]] * dispersal
  mu <- pmin(pmax(inv_logit(eta), 1e-8), 1 - 1e-8)

  result <- list(mu = mu, dispersal = dispersal, auto = auto)
  if (include_kappa_derivative) {
    result$dispersal_kappa_derivative <- dispersal_by_candidate(
      y_previous,
      wind_matrix,
      distance_matrix,
      candidate_indicator,
      kappa,
      d0,
      derivative = TRUE
    )
  }
  result
}

zib_loglik_matrix <- function(y, mu, phi) {
  if (is.null(dim(mu))) mu <- matrix(mu, ncol = 1L)
  alpha <- mean(y == 0)
  result <- matrix(0, nrow = length(y), ncol = ncol(mu))
  zero <- y == 0
  positive <- !zero

  if (any(zero)) result[zero, ] <- log(alpha)
  if (any(positive)) {
    result[positive, ] <- log1p(-alpha) + dbeta(
      y[positive],
      shape1 = mu[positive, , drop = FALSE] * phi,
      shape2 = (1 - mu[positive, , drop = FALSE]) * phi,
      log = TRUE
    )
  }
  result
}

transition_logliks <- function(
  theta,
  y_current,
  y_previous,
  wind_matrix,
  distance_matrix,
  candidate_indicator,
  d0 = 0.01
) {
  components <- transition_components(
    theta,
    y_previous,
    wind_matrix,
    distance_matrix,
    candidate_indicator,
    d0
  )
  colSums(zib_loglik_matrix(
    y_current,
    components$mu,
    exp(theta[["log_phi"]])
  ))
}

neg_expected_transition_loglik <- function(
  theta,
  y_current,
  y_previous,
  wind_matrix,
  distance_matrix,
  candidate_indicator,
  responsibilities,
  d0 = 0.01
) {
  -sum(responsibilities * transition_logliks(
    theta,
    y_current,
    y_previous,
    wind_matrix,
    distance_matrix,
    candidate_indicator,
    d0
  ))
}

neg_expected_transition_gradient <- function(
  theta,
  y_current,
  y_previous,
  wind_matrix,
  distance_matrix,
  candidate_indicator,
  responsibilities,
  d0 = 0.01
) {
  phi <- exp(theta[["log_phi"]])
  kappa <- exp(theta[["log_kappa"]])
  components <- transition_components(
    theta,
    y_previous,
    wind_matrix,
    distance_matrix,
    candidate_indicator,
    d0,
    include_kappa_derivative = TRUE
  )

  positive <- which(y_current > 0)
  mu <- components$mu[positive, , drop = FALSE]
  dispersal <- components$dispersal[positive, , drop = FALSE]
  dispersal_derivative <-
    components$dispersal_kappa_derivative[positive, , drop = FALSE]
  auto <- components$auto[positive]
  y <- matrix(y_current[positive], nrow = length(positive), ncol = ncol(mu))
  weights <- matrix(
    responsibilities,
    nrow = length(positive),
    ncol = length(responsibilities),
    byrow = TRUE
  )

  score_eta <- phi *
    (logit(y) - digamma(mu * phi) + digamma((1 - mu) * phi)) *
    mu * (1 - mu)
  weighted_score <- weights * score_eta

  d_beta <- sum(weighted_score)
  d_delta <- sum(weighted_score * auto)
  d_gamma <- sum(weighted_score * dispersal)
  d_log_kappa <- sum(
    weighted_score * theta[["gamma"]] * dispersal_derivative
  ) * kappa

  score_phi <- digamma(phi) -
    mu * digamma(mu * phi) -
    (1 - mu) * digamma((1 - mu) * phi) +
    mu * log(y) + (1 - mu) * log1p(-y)
  d_log_phi <- sum(weights * score_phi) * phi

  -c(
    beta = d_beta,
    delta = d_delta,
    gamma = d_gamma,
    log_kappa = d_log_kappa,
    log_phi = d_log_phi
  )
}

initialize_transition_theta <- function(
  y_current,
  y_previous,
  wind_matrix,
  distance_matrix,
  initial_kappa,
  d0 = 0.01
) {
  positive <- which(y_current > 0)
  if (length(positive) < 4L) {
    stop("At least four positive observations are needed for initialization.")
  }

  full_indicator <- matrix(1, nrow = length(y_previous), ncol = 1L)
  spatial <- dispersal_by_candidate(
    y_previous,
    wind_matrix,
    distance_matrix,
    full_indicator,
    initial_kappa,
    d0
  )[, 1]
  initial_fit <- lm(
    logit(pmin(pmax(y_current[positive], 1e-5), 1 - 1e-5)) ~
      I(y_previous[positive] * (1 - y_previous[positive])) +
      spatial[positive]
  )
  coefficients <- coef(initial_fit)
  if (any(!is.finite(coefficients))) {
    stop("Method-of-moments regression produced non-finite coefficients.")
  }

  mu <- inv_logit(fitted(initial_fit))
  residual_variance <- mean(resid(initial_fit)^2 * (mu * (1 - mu))^2)
  phi <- mean(mu * (1 - mu)) / max(residual_variance, 1e-6) - 1
  phi <- max(phi, 1e-3)

  c(
    beta = unname(coefficients[1]),
    delta = unname(coefficients[2]),
    gamma = unname(coefficients[3]),
    log_kappa = log(initial_kappa),
    log_phi = log(phi)
  )
}

# Construct transition-specific method-of-moments starting values without
# fitting the full forward model. This is used by the initialization
# sensitivity analysis so the subsequent EM fit begins directly from the
# specified kappa value.
initialize_experiment_theta <- function(
  block,
  treat,
  mod_dat,
  initial_kappa,
  through_visit = NULL,
  transition_visit = NULL,
  d0 = 0.01
) {
  experiment <- extract_experiment(
    mod_dat, block, treat, through_visit, transition_visit
  )

  theta <- do.call(rbind, lapply(
    seq_along(experiment$transition_visits),
    function(t) {
      initialize_transition_theta(
        y_current = experiment$y[[t + 1L]],
        y_previous = experiment$y[[t]],
        wind_matrix = experiment$wind[[t]],
        distance_matrix = experiment$distance,
        initial_kappa = initial_kappa,
        d0 = d0
      )
    }
  ))
  rownames(theta) <- experiment$transition_visits
  theta
}

# Numerically evaluate the observed information at a fitted forward-model
# optimum. The Hessian is computed from the objective rather than the custom
# gradient so uncertainty diagnostics do not depend on the analytic gradient.
forward_hessian_diagnostics <- function(
  theta,
  y_current,
  y_previous,
  wind_matrix,
  distance_matrix,
  d0 = 0.01
) {
  full_indicator <- matrix(1, nrow = length(y_previous), ncol = 1L)
  hessian <- tryCatch(
    optimHess(
      par = theta,
      fn = neg_expected_transition_loglik,
      y_current = y_current,
      y_previous = y_previous,
      wind_matrix = wind_matrix,
      distance_matrix = distance_matrix,
      candidate_indicator = full_indicator,
      responsibilities = 1,
      d0 = d0
    ),
    error = function(e) NULL
  )

  empty_result <- list(
    hessian = NULL,
    vcov = NULL,
    rank = NA_integer_,
    minimum_eigenvalue = NA_real_,
    maximum_eigenvalue = NA_real_,
    condition_number = NA_real_,
    positive_definite = FALSE,
    inversion_succeeded = FALSE,
    usable_vcov = FALSE,
    standard_errors = setNames(rep(NA_real_, length(theta)), names(theta)),
    gamma_log_kappa_correlation = NA_real_
  )
  if (is.null(hessian) || any(!is.finite(hessian))) return(empty_result)

  hessian <- (hessian + t(hessian)) / 2
  dimnames(hessian) <- list(names(theta), names(theta))
  eigenvalues <- tryCatch(
    eigen(hessian, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) rep(NA_real_, nrow(hessian))
  )
  scale <- max(abs(eigenvalues), na.rm = TRUE)
  eigen_tolerance <- sqrt(.Machine$double.eps) * max(scale, 1)
  positive_definite <- all(is.finite(eigenvalues)) &&
    min(eigenvalues) > eigen_tolerance
  hessian_rank <- sum(abs(eigenvalues) > eigen_tolerance)
  condition_number <- tryCatch(
    kappa(hessian, exact = TRUE),
    error = function(e) NA_real_
  )

  covariance <- tryCatch(
    chol2inv(chol(hessian)),
    error = function(e) NULL
  )
  if (!is.null(covariance)) {
    dimnames(covariance) <- list(names(theta), names(theta))
  }
  inversion_succeeded <- !is.null(covariance)
  usable_vcov <- inversion_succeeded &&
    positive_definite &&
    all(is.finite(covariance)) &&
    all(diag(covariance) > 0)
  standard_errors <- setNames(rep(NA_real_, length(theta)), names(theta))
  gamma_log_kappa_correlation <- NA_real_
  if (usable_vcov) {
    standard_errors <- sqrt(diag(covariance))
    gamma_log_kappa_correlation <- covariance["gamma", "log_kappa"] /
      (standard_errors["gamma"] * standard_errors["log_kappa"])
  }

  list(
    hessian = hessian,
    vcov = covariance,
    rank = hessian_rank,
    minimum_eigenvalue = min(eigenvalues),
    maximum_eigenvalue = max(eigenvalues),
    condition_number = condition_number,
    positive_definite = positive_definite,
    inversion_succeeded = inversion_succeeded,
    usable_vcov = usable_vcov,
    standard_errors = standard_errors,
    gamma_log_kappa_correlation = gamma_log_kappa_correlation
  )
}

extract_experiment <- function(
  mod_dat,
  block,
  treat,
  through_visit = NULL,
  transition_visit = NULL
) {
  all_visits <- dimnames(mod_dat$intensity)[["visit"]]
  if (!is.null(through_visit) && !is.null(transition_visit)) {
    stop("Specify through_visit or transition_visit, not both.")
  }
  if (!is.null(transition_visit)) {
    endpoint <- match(as.character(transition_visit), all_visits)
    if (is.na(endpoint) || endpoint < 2L) {
      stop("transition_visit must identify a visit after the first visit.")
    }
    visits <- all_visits[c(endpoint - 1L, endpoint)]
    transition_visits <- all_visits[endpoint]
    return(list(
      through_visit = as.character(transition_visit),
      final_study_visit = tail(all_visits, 1L),
      visits = visits,
      transition_visits = transition_visits,
      y = lapply(visits, function(v) mod_dat$intensity[, block, treat, v]),
      wind = list(mod_dat$wind[, , block, treat, transition_visits]),
      distance = mod_dat$dist
    ))
  }
  if (is.null(through_visit)) through_visit <- tail(all_visits, 1L)
  endpoint <- match(as.character(through_visit), all_visits)
  if (is.na(endpoint) || endpoint < 2L) {
    stop("through_visit must identify a visit after the first visit.")
  }
  visits <- all_visits[seq_len(endpoint)]
  transition_visits <- visits[-1]
  list(
    through_visit = as.character(through_visit),
    final_study_visit = tail(all_visits, 1L),
    visits = visits,
    transition_visits = transition_visits,
    y = lapply(visits, function(v) mod_dat$intensity[, block, treat, v]),
    wind = lapply(transition_visits, function(v) {
      mod_dat$wind[, , block, treat, v]
    }),
    distance = mod_dat$dist
  )
}

add_forward_hessian_diagnostics <- function(
  forward_fit,
  block,
  treat,
  mod_dat,
  d0 = 0.01
) {
  experiment <- extract_experiment(mod_dat, block, treat)
  diagnostics <- lapply(seq_along(experiment$transition_visits), function(t) {
    forward_hessian_diagnostics(
      theta = forward_fit$theta[t, ],
      y_current = experiment$y[[t + 1L]],
      y_previous = experiment$y[[t]],
      wind_matrix = experiment$wind[[t]],
      distance_matrix = experiment$distance,
      d0 = d0
    )
  })

  for (t in seq_along(diagnostics)) {
    forward_fit$fits[[t]]$hessian_diagnostics <- diagnostics[[t]]
  }
  standard_errors <- do.call(rbind, lapply(
    diagnostics, `[[`, "standard_errors"
  ))
  forward_fit$summary$hessian_rank <- vapply(
    diagnostics, `[[`, integer(1), "rank"
  )
  forward_fit$summary$hessian_minimum_eigenvalue <- vapply(
    diagnostics, `[[`, numeric(1), "minimum_eigenvalue"
  )
  forward_fit$summary$hessian_maximum_eigenvalue <- vapply(
    diagnostics, `[[`, numeric(1), "maximum_eigenvalue"
  )
  forward_fit$summary$hessian_condition_number <- vapply(
    diagnostics, `[[`, numeric(1), "condition_number"
  )
  forward_fit$summary$hessian_positive_definite <- vapply(
    diagnostics, `[[`, logical(1), "positive_definite"
  )
  forward_fit$summary$hessian_inversion_succeeded <- vapply(
    diagnostics, `[[`, logical(1), "inversion_succeeded"
  )
  forward_fit$summary$hessian_usable_vcov <- vapply(
    diagnostics, `[[`, logical(1), "usable_vcov"
  )
  forward_fit$summary$se_beta <- standard_errors[, "beta"]
  forward_fit$summary$se_delta <- standard_errors[, "delta"]
  forward_fit$summary$se_gamma <- standard_errors[, "gamma"]
  forward_fit$summary$se_log_kappa <- standard_errors[, "log_kappa"]
  forward_fit$summary$se_log_phi <- standard_errors[, "log_phi"]
  forward_fit$summary$se_kappa <- forward_fit$summary$kappa *
    forward_fit$summary$se_log_kappa
  forward_fit$summary$se_phi <- forward_fit$summary$phi *
    forward_fit$summary$se_log_phi
  forward_fit$summary$gamma_log_kappa_correlation <- vapply(
    diagnostics, `[[`, numeric(1), "gamma_log_kappa_correlation"
  )
  forward_fit$hessian_diagnostics <- diagnostics
  forward_fit
}

fit_forward_experiment <- function(
  block,
  treat,
  mod_dat,
  kappa_grid = exp(seq(log(0.5), log(2), length.out = 5)),
  d0 = 0.01
) {
  experiment <- extract_experiment(mod_dat, block, treat)
  full_indicator <- matrix(1, nrow = nrow(mod_dat$dist), ncol = 1L)

  fits <- lapply(seq_along(experiment$transition_visits), function(t) {
    attempts <- lapply(kappa_grid, function(initial_kappa) {
      initial <- tryCatch(
        initialize_transition_theta(
          experiment$y[[t + 1L]], experiment$y[[t]], experiment$wind[[t]],
          experiment$distance, initial_kappa, d0
        ),
        error = function(e) NULL
      )
      if (is.null(initial)) return(NULL)

      fit <- tryCatch(optim(
        par = initial,
        fn = neg_expected_transition_loglik,
        gr = neg_expected_transition_gradient,
        method = "BFGS",
        control = list(maxit = 5000, reltol = 1e-8),
        y_current = experiment$y[[t + 1L]],
        y_previous = experiment$y[[t]],
        wind_matrix = experiment$wind[[t]],
        distance_matrix = experiment$distance,
        candidate_indicator = full_indicator,
        responsibilities = 1,
        d0 = d0
      ), error = function(e) NULL)
      if (is.null(fit) || !is.finite(fit$value)) return(NULL)
      list(fit = fit, initial_kappa = initial_kappa)
    })
    attempts <- Filter(Negate(is.null), attempts)
    converged <- Filter(function(x) x$fit$convergence == 0L, attempts)
    eligible <- if (length(converged)) converged else attempts
    if (!length(eligible)) {
      stop("No finite forward fit for block ", block, ", treatment ", treat,
           ", transition ending at visit ", experiment$transition_visits[t], ".")
    }
    eligible[[which.min(vapply(eligible, function(x) x$fit$value, numeric(1)))]]
  })

  theta <- do.call(rbind, lapply(fits, function(x) x$fit$par))
  rownames(theta) <- experiment$transition_visits
  summary <- data.frame(
    block = block,
    treat = treat,
    visit = experiment$transition_visits,
    beta = theta[, "beta"],
    delta = theta[, "delta"],
    gamma = theta[, "gamma"],
    kappa = exp(theta[, "log_kappa"]),
    phi = exp(theta[, "log_phi"]),
    neg_loglik = vapply(fits, function(x) x$fit$value, numeric(1)),
    initial_kappa = vapply(fits, function(x) x$initial_kappa, numeric(1)),
    convergence = vapply(fits, function(x) x$fit$convergence, integer(1))
  )

  add_forward_hessian_diagnostics(
    list(theta = theta, summary = summary, fits = fits),
    block = block,
    treat = treat,
    mod_dat = mod_dat,
    d0 = d0
  )
}

candidate_loglik_matrix <- function(
  theta_by_transition,
  experiment,
  candidate_indicator,
  d0 = 0.01
) {
  result <- matrix(
    NA_real_,
    nrow = ncol(candidate_indicator),
    ncol = length(experiment$transition_visits),
    dimnames = list(
      candidate = colnames(candidate_indicator),
      transition = experiment$transition_visits
    )
  )
  for (t in seq_along(experiment$transition_visits)) {
    result[, t] <- transition_logliks(
      theta_by_transition[t, ],
      experiment$y[[t + 1L]],
      experiment$y[[t]],
      experiment$wind[[t]],
      experiment$distance,
      candidate_indicator,
      d0
    )
  }
  result
}

responsibilities_from_loglik <- function(loglik, log_prior) {
  values <- log_prior + loglik
  exp(values - log_sum_exp(values))
}

make_posterior_summaries <- function(loglik_matrix, log_prior, candidate_ids) {
  transition <- do.call(rbind, lapply(seq_len(ncol(loglik_matrix)), function(t) {
    posterior <- responsibilities_from_loglik(loglik_matrix[, t], log_prior)
    data.frame(
      visit = colnames(loglik_matrix)[t], candidate_id = candidate_ids,
      loglik = loglik_matrix[, t], posterior = posterior,
      rank = rank(-posterior, ties.method = "min"),
      predicted = posterior == max(posterior),
      evidence = "transition"
    )
  }))
  cumulative <- do.call(rbind, lapply(seq_len(ncol(loglik_matrix)), function(t) {
    cumulative_loglik <- rowSums(loglik_matrix[, seq_len(t), drop = FALSE])
    posterior <- responsibilities_from_loglik(cumulative_loglik, log_prior)
    data.frame(
      visit = colnames(loglik_matrix)[t], candidate_id = candidate_ids,
      loglik = cumulative_loglik, posterior = posterior,
      rank = rank(-posterior, ties.method = "min"),
      predicted = posterior == max(posterior),
      evidence = "cumulative"
    )
  }))
  rbind(transition, cumulative)
}

fit_source_detection_experiment <- function(
  block,
  treat,
  config,
  n_sources,
  mod_dat,
  theta_initial,
  through_visit = NULL,
  transition_visit = NULL,
  max_iterations = 200L,
  loglik_tolerance = 1e-3,
  responsibility_tolerance = 1e-3,
  monotonicity_tolerance = 1e-8,
  d0 = 0.01,
  include_naive = TRUE
) {
  experiment <- extract_experiment(
    mod_dat, block, treat, through_visit, transition_visit
  )
  fit_scope <- if (!is.null(transition_visit)) {
    "transition_refit"
  } else if (experiment$through_visit == experiment$final_study_visit) {
    "whole_study"
  } else {
    "cumulative_refit"
  }
  candidates <- make_candidate_sets(mod_dat$groups[, config], n_sources)
  true_sources <- get_true_source_groups(mod_dat, block, treat, config)
  true_n_sources <- length(true_sources)
  true_candidate_id <- paste(true_sources, collapse = "_")
  log_prior <- rep(-log(length(candidates$ids)), length(candidates$ids))
  theta_current <- theta_initial[experiment$transition_visits, , drop = FALSE]

  loglik_matrix <- candidate_loglik_matrix(
    theta_current, experiment, candidates$indicator, d0
  )
  candidate_loglik <- rowSums(loglik_matrix)
  responsibilities <- responsibilities_from_loglik(candidate_loglik, log_prior)
  observed_loglik <- log_sum_exp(log_prior + candidate_loglik)

  history <- data.frame(
    iteration = 0L, observed_loglik = observed_loglik,
    loglik_change = NA_real_, max_parameter_change = NA_real_,
    max_responsibility_change = NA_real_, m_step_convergence = NA_integer_,
    monotonic = NA, converged = FALSE
  )
  responsibility_history <- data.frame(
    iteration = 0L, candidate_id = candidates$ids,
    loglik = candidate_loglik, responsibility = responsibilities,
    rank = rank(-responsibilities, ties.method = "min"),
    predicted = responsibilities == max(responsibilities)
  )
  status <- "maximum_iterations"

  for (iteration in seq_len(max_iterations)) {
    fits <- lapply(seq_along(experiment$transition_visits), function(t) {
      optim(
        par = theta_current[t, ],
        fn = neg_expected_transition_loglik,
        gr = neg_expected_transition_gradient,
        method = "BFGS",
        control = list(maxit = 1000, reltol = 1e-8),
        y_current = experiment$y[[t + 1L]],
        y_previous = experiment$y[[t]],
        wind_matrix = experiment$wind[[t]],
        distance_matrix = experiment$distance,
        candidate_indicator = candidates$indicator,
        responsibilities = responsibilities,
        d0 = d0
      )
    })
    theta_new <- do.call(rbind, lapply(fits, `[[`, "par"))
    rownames(theta_new) <- experiment$transition_visits
    m_convergence <- max(vapply(fits, `[[`, integer(1), "convergence"))

    loglik_matrix_new <- candidate_loglik_matrix(
      theta_new, experiment, candidates$indicator, d0
    )
    candidate_loglik_new <- rowSums(loglik_matrix_new)
    responsibilities_new <- responsibilities_from_loglik(
      candidate_loglik_new, log_prior
    )
    observed_loglik_new <- log_sum_exp(log_prior + candidate_loglik_new)
    loglik_change <- observed_loglik_new - observed_loglik
    parameter_change <- max(abs(theta_new - theta_current))
    responsibility_change <- max(abs(responsibilities_new - responsibilities))
    monotonic <- loglik_change >= -monotonicity_tolerance
    converged <- abs(loglik_change) <= loglik_tolerance &&
      responsibility_change <= responsibility_tolerance &&
      m_convergence == 0L

    history <- rbind(history, data.frame(
      iteration = iteration, observed_loglik = observed_loglik_new,
      loglik_change = loglik_change,
      max_parameter_change = parameter_change,
      max_responsibility_change = responsibility_change,
      m_step_convergence = m_convergence,
      monotonic = monotonic, converged = converged
    ))
    responsibility_history <- rbind(responsibility_history, data.frame(
      iteration = iteration, candidate_id = candidates$ids,
      loglik = candidate_loglik_new, responsibility = responsibilities_new,
      rank = rank(-responsibilities_new, ties.method = "min"),
      predicted = responsibilities_new == max(responsibilities_new)
    ))

    theta_current <- theta_new
    loglik_matrix <- loglik_matrix_new
    candidate_loglik <- candidate_loglik_new
    responsibilities <- responsibilities_new
    observed_loglik <- observed_loglik_new

    if (!monotonic) {
      status <- "observed_loglik_decreased"
      break
    }
    if (converged) {
      status <- "converged"
      break
    }
  }

  source_columns <- as.data.frame(t(candidates$matrix))
  names(source_columns) <- paste0("source", seq_len(n_sources))
  candidate_summary <- cbind(data.frame(
    block = block, treat = treat, config = config, n_sources = n_sources,
    through_visit = experiment$through_visit, fit_scope = fit_scope,
    candidate_id = candidates$ids, loglik = candidate_loglik,
    posterior = responsibilities,
    rank = rank(-responsibilities, ties.method = "min"),
    predicted = responsibilities == max(responsibilities),
    true_candidate_id = true_candidate_id,
    true_candidate = candidates$ids == true_candidate_id
  ), source_columns)
  parameter_summary <- data.frame(
    block = block, treat = treat, config = config, n_sources = n_sources,
    through_visit = experiment$through_visit, fit_scope = fit_scope,
    visit = experiment$transition_visits,
    beta = theta_current[, "beta"], delta = theta_current[, "delta"],
    gamma = theta_current[, "gamma"],
    kappa = exp(theta_current[, "log_kappa"]),
    phi = exp(theta_current[, "log_phi"]),
    status = status, iterations = max(history$iteration),
    observed_loglik = observed_loglik
  )
  diagnostic_posteriors <- make_posterior_summaries(
    loglik_matrix, log_prior, candidates$ids
  )
  diagnostic_posteriors$block <- block
  diagnostic_posteriors$treat <- treat
  diagnostic_posteriors$config <- config
  diagnostic_posteriors$n_sources <- n_sources
  diagnostic_posteriors$through_visit <- experiment$through_visit
  diagnostic_posteriors$fit_scope <- fit_scope

  predicted_sources <- candidates$matrix[, which.max(responsibilities)]
  fitted_accuracy <- source_accuracy_general(
    predicted_sources,
    true_sources,
    mod_dat$grid_dist[[config]]
  )
  if (include_naive) {
    naive_sources <- naive_source_prediction(
      experiment,
      mod_dat$groups[, config],
      n_sources
    )
    naive_accuracy <- source_accuracy_general(
      naive_sources,
      true_sources,
      mod_dat$grid_dist[[config]]
    )
  } else {
    naive_sources <- numeric(0)
    naive_accuracy <- list(
      source_count_correct = NA,
      standard_accuracy = NA_real_,
      n_correct = NA_integer_,
      precision = NA_real_,
      recall = NA_real_,
      f1 = NA_real_,
      mean_distance = NA_real_,
      distance_weighted_accuracy = NA_real_,
      assignment = NULL
    )
  }
  accuracy_summary <- data.frame(
    block = block, treat = treat, config = config, n_sources = n_sources,
    true_n_sources = true_n_sources,
    through_visit = experiment$through_visit, fit_scope = fit_scope,
    evidence = fit_scope, visit = experiment$through_visit,
    predicted_candidate_id = candidates$ids[which.max(responsibilities)],
    true_candidate_id = true_candidate_id,
    naive_candidate_id = if (include_naive) {
      paste(naive_sources, collapse = "_")
    } else {
      NA_character_
    },
    source_count_correct = fitted_accuracy$source_count_correct,
    standard_accuracy = fitted_accuracy$standard_accuracy,
    n_correct = fitted_accuracy$n_correct,
    precision = fitted_accuracy$precision,
    recall = fitted_accuracy$recall,
    f1 = fitted_accuracy$f1,
    mean_distance = fitted_accuracy$mean_distance,
    distance_weighted_accuracy = fitted_accuracy$distance_weighted_accuracy,
    naive_source_count_correct = naive_accuracy$source_count_correct,
    naive_standard_accuracy = naive_accuracy$standard_accuracy,
    naive_n_correct = naive_accuracy$n_correct,
    naive_precision = naive_accuracy$precision,
    naive_recall = naive_accuracy$recall,
    naive_f1 = naive_accuracy$f1,
    naive_mean_distance = naive_accuracy$mean_distance,
    naive_distance_weighted_accuracy =
      naive_accuracy$distance_weighted_accuracy,
    observed_loglik = observed_loglik
  )
  diagnostic_groups <- split(
    diagnostic_posteriors,
    interaction(
      diagnostic_posteriors$evidence,
      diagnostic_posteriors$visit,
      drop = TRUE
    )
  )
  diagnostic_accuracy <- do.call(rbind, lapply(diagnostic_groups, function(x) {
    winner <- x[which.max(x$posterior), ]
    winner_index <- match(winner$candidate_id, candidates$ids)
    metrics <- source_accuracy_general(
      candidates$matrix[, winner_index],
      true_sources,
      mod_dat$grid_dist[[config]]
    )
    if (include_naive) {
      diagnostic_experiment <- extract_experiment(
        mod_dat,
        block,
        treat,
        through_visit = if (winner$evidence == "cumulative") {
          winner$visit
        } else {
          NULL
        },
        transition_visit = if (winner$evidence == "transition") {
          winner$visit
        } else {
          NULL
        }
      )
      naive_diagnostic_sources <- naive_source_prediction(
        diagnostic_experiment,
        mod_dat$groups[, config],
        n_sources
      )
      naive_metrics <- source_accuracy_general(
        naive_diagnostic_sources,
        true_sources,
        mod_dat$grid_dist[[config]]
      )
    } else {
      naive_diagnostic_sources <- numeric(0)
      naive_metrics <- naive_accuracy
    }
    data.frame(
      block = block, treat = treat, config = config, n_sources = n_sources,
      true_n_sources = true_n_sources,
      through_visit = experiment$through_visit, fit_scope = fit_scope,
      evidence = winner$evidence, visit = as.character(winner$visit),
      predicted_candidate_id = winner$candidate_id,
      true_candidate_id = true_candidate_id,
      naive_candidate_id = if (include_naive) {
        paste(naive_diagnostic_sources, collapse = "_")
      } else {
        NA_character_
      },
      source_count_correct = metrics$source_count_correct,
      standard_accuracy = metrics$standard_accuracy,
      n_correct = metrics$n_correct,
      precision = metrics$precision,
      recall = metrics$recall,
      f1 = metrics$f1,
      mean_distance = metrics$mean_distance,
      distance_weighted_accuracy = metrics$distance_weighted_accuracy,
      naive_source_count_correct = naive_metrics$source_count_correct,
      naive_standard_accuracy = naive_metrics$standard_accuracy,
      naive_n_correct = naive_metrics$n_correct,
      naive_precision = naive_metrics$precision,
      naive_recall = naive_metrics$recall,
      naive_f1 = naive_metrics$f1,
      naive_mean_distance = naive_metrics$mean_distance,
      naive_distance_weighted_accuracy =
        naive_metrics$distance_weighted_accuracy,
      observed_loglik = observed_loglik
    )
  }))
  accuracy_summary <- rbind(accuracy_summary, diagnostic_accuracy)
  assignment_details <- fitted_accuracy$assignment
  assignment_details$method <- "model"
  if (include_naive) {
    naive_assignment_details <- naive_accuracy$assignment
    naive_assignment_details$method <- "naive_maximum_severity"
    assignment_details <- rbind(assignment_details, naive_assignment_details)
  }
  assignment_details$block <- block
  assignment_details$treat <- treat
  assignment_details$config <- config
  assignment_details$n_sources <- n_sources
  assignment_details$true_n_sources <- true_n_sources
  assignment_details$through_visit <- experiment$through_visit
  assignment_details$fit_scope <- fit_scope

  parameter_summary$true_n_sources <- true_n_sources

  list(
    candidate_summary = candidate_summary,
    parameter_summary = parameter_summary,
    responsibility_history = transform(
      responsibility_history,
      block = block, treat = treat, config = config, n_sources = n_sources,
      through_visit = experiment$through_visit, fit_scope = fit_scope
    ),
    diagnostic_posteriors = diagnostic_posteriors,
    accuracy_summary = accuracy_summary,
    whole_study_assignment = assignment_details,
    history = transform(
      history,
      block = block, treat = treat, config = config, n_sources = n_sources,
      through_visit = experiment$through_visit, fit_scope = fit_scope
    ),
    theta = theta_current,
    candidate_sets = candidates$matrix,
    transition_loglik = loglik_matrix
  )
}
