#' Diagnose Selection Conditioning Sensitivity
#'
#' @description
#' Describes how integrating the contexts retained by a fitted selection model
#' would reweight their pre-selection Gaussian distribution. This compares
#' reporting models; it is not a numerical approximation error bound.
#'
#' @param object A fitted normal-outcome model with a bound selection model.
#' @param max_posterior_samples Maximum number of posterior draws to evaluate.
#'   Draws are deterministically thinned over posterior row order.
#' @param latent_samples Number of fresh pre-selection retained-context draws
#'   per posterior draw. Defaults to `512`.
#' @param seed Integer seed for the diagnostic simulation. The calling R
#'   session's random-number state is restored before returning.
#' @param integration_control Optional settings from
#'   [set_selection_likelihood_control()] for diagnostic event integration.
#'   `NULL` uses the fitted model's settings.
#'
#' @details
#' At each posterior draw the fitted model defines a conditional Gaussian
#' response law, retained contexts, and any best-rule publication groups. The
#' diagnostic evaluates each Gaussian selection unit's normalizer at fresh
#' retained contexts, preserving the fitted product or best weight rule and any
#' best-rule publication groups. Product-only models require no grouping.
#' Within a connected source unit, their product gives the
#' relative reweighting when the retained contexts are integrated before
#' selection normalization.
#'
#' `ess_fraction` estimates `mean(W)^2 / mean(W^2)`;
#' `total_variation` estimates the distance between the pre-selection context
#' law and its normalized `W`-reweighted law. The fitted and reference
#' conditioning specifications, publication identities, and source roles are
#' retained in the `target` attribute. With no retained context this comparison
#' is inapplicable and its reweighting is identically constant.
#'
#' The fitted weight rule is held fixed. These summaries do not compare
#' arbitrary weight rules, identify the correct reporting mechanism, or reveal
#' parameter regions absent from the fitted posterior.
#'
#' Simulation budgets control diagnostic Monte Carlo precision separately from
#' the fitted Gaussian event integration controls. The largest unit-level
#' posterior-median total variation distance is silent up to `0.05`, produces
#' a message above `0.05` and below `0.10`, and a warning from `0.10`.
#' These notifications describe sensitivity rather than statistical adequacy.
#'
#' This optional comparison is computed only when requested; fitting does not
#' compute it. Calls without simulation arguments reuse an explicitly stored
#' result, and explicit settings reuse it only when they match. Use
#' `add_selection_sensitivity_diagnostics()` to compute and store a result.
#' Chain extension silently refreshes only stored diagnostics using their stored
#' simulation settings; label-only updates retain them. Failed refreshes remain
#' stored as unavailable without discarding the fitted posterior. Call
#' `selection_sensitivity_diagnostics()` to inspect the comparison or its error.
#' Sensitivity notifications are emitted only by these explicit helper calls,
#' not by fitting, extension, or ordinary model printing.
#'
#' @return A data frame with one row per connected source unit. The `q05`,
#'   `median`, and `q95` columns summarize ESS fractions, total variation, and
#'   the interquartile range of log weights across posterior draws. Conditional
#'   Monte Carlo standard errors and simulation budgets are also reported.
#'   `settings` and `target` attributes record simulation settings and the
#'   compared model specifications. The add function returns the fitted object.
#' @seealso [bselmodel.mv()], [set_selection_likelihood_control()]
#' @export
selection_sensitivity_diagnostics <- function(
    object, max_posterior_samples = 256L, latent_samples = 512L,
    seed = 1L, integration_control = NULL) {

  if (!inherits(object, "brma") || is.null(.data_selection_model(object[["data"]]))) {
    stop("'object' must be a fitted model with a bound selection model.", call. = FALSE)
  }
  if (is.null(object[["fit"]])) {
    stop("'object' must contain a fitted posterior.", call. = FALSE)
  }

  cached <- object[["selection_sensitivity_diagnostics"]]
  if (!is.null(cached) && missing(max_posterior_samples) &&
      missing(latent_samples) && missing(seed) && missing(integration_control)) {
    if (inherits(cached, "error")) {
      stop(cached)
    }
    .selection_sensitivity_notify_result(cached)
    return(cached)
  }

  max_posterior_samples <- .normalize_max_samples(
    max_posterior_samples,
    argument = "max_posterior_samples"
  )
  latent_samples <- .normalize_max_samples(
    latent_samples,
    argument = "latent_samples",
    minimum = 2L
  )
  if (is.infinite(latent_samples)) {
    stop("'latent_samples' must be a finite positive integer.", call. = FALSE)
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed) || seed < 0 || seed > .Machine$integer.max ||
      seed != floor(seed)) {
    stop("'seed' must be a single non-negative integer.", call. = FALSE)
  }
  seed <- as.integer(seed)
  if (!is.null(integration_control)) {
    integration_control <- .check_selection_likelihood_control(integration_control)
  }
  settings <- list(
    max_posterior_samples = max_posterior_samples,
    latent_samples        = latent_samples,
    seed                  = seed,
    integration_control   = integration_control
  )
  if (!is.null(cached) && identical(attr(cached, "settings"), settings)) {
    if (inherits(cached, "error")) {
      stop(cached)
    }
    .selection_sensitivity_notify_result(cached)
    return(cached)
  }

  result <- .compute_selection_sensitivity_diagnostics(
    object                = object,
    max_posterior_samples = max_posterior_samples,
    latent_samples        = latent_samples,
    seed                  = seed,
    integration_control   = integration_control
  )
  .selection_sensitivity_notify_result(result)
  return(result)
}


#' @rdname selection_sensitivity_diagnostics
#' @param ... Simulation arguments passed to
#'   `selection_sensitivity_diagnostics()`.
#' @export
add_selection_sensitivity_diagnostics <- function(object, ...) {

  object[["selection_sensitivity_diagnostics"]] <- NULL
  object[["selection_sensitivity_diagnostics"]] <-
    selection_sensitivity_diagnostics(object, ...)
  return(object)
}


.refresh_selection_sensitivity_diagnostics <- function(object) {

  if (is.null(object[["selection_sensitivity_diagnostics"]]) ||
      is.null(.data_selection_model(object[["data"]]))) {
    return(object)
  }

  settings <- attr(object[["selection_sensitivity_diagnostics"]], "settings")
  object[["selection_sensitivity_diagnostics"]] <- NULL
  result <- tryCatch(
    do.call(
      .compute_selection_sensitivity_diagnostics,
      c(list(object = object), settings)
    ),
    error = function(error) error
  )
  if (inherits(result, "error")) {
    attr(result, "settings") <- settings
  }
  object[["selection_sensitivity_diagnostics"]] <- result
  return(object)
}


.selection_sensitivity_notify_result <- function(result) {

  if (is.null(result)) {
    return(invisible(NULL))
  }
  if (inherits(result, "error")) {
    stop(result)
  }
  .selection_sensitivity_notify(result[["total_variation_median"]])
  invisible(NULL)
}


.compute_selection_sensitivity_diagnostics <- function(
    object, max_posterior_samples = 256L, latent_samples = 512L, seed = 1L,
    integration_control = NULL) {

  rng_kind <- RNGkind()
  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (has_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    do.call(RNGkind, as.list(rng_kind))
    if (has_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)

  sample_info <- .known_v_diagnostic_posterior_samples(
    object = object, max_samples = max_posterior_samples,
    caller = "selection_sensitivity_diagnostics()", warn = FALSE
  )
  posterior_samples <- sample_info[["posterior_samples"]]
  data <- object[["data"]]
  priors <- object[["priors"]]
  model <- .data_selection_model(data)
  if (is.null(model)) {
    stop("Selection sensitivity diagnostics are unavailable without a bound selection model.",
      call. = FALSE)
  }
  means <- .evaluate.brma.mu(
    fit = object[["fit"]], outcome_data = data[["outcome"]],
    mods_data = data[["mods"]],
    mods_formula = if (.is_mods(object)) .create_fit_formula_list(data, "mods") else NULL,
    mods_priors = if (.is_random(object)) priors[["location"]] else priors[["mods"]],
    priors = priors, is_mods = .is_mods(object),
    is_PET = .is_PET(object), is_PEESE = .is_PEESE(object),
    effect_direction = .effect_direction(object), bias_adjusted = FALSE,
    K = nrow(data[["outcome"]]), posterior_samples = posterior_samples
  )
  context <- .selection_context(object, posterior_samples)
  .selection_require_step_evaluable(context, "selection_sensitivity_diagnostics()")
  scales <- if (.is_data_random(data)) {
    list(tau_within = matrix(0, nrow(means), ncol(means)),
         tau_between = matrix(0, nrow(means), ncol(means)))
  } else {
    .zplot_tau_samples(object, posterior_samples)
  }
  parts <- .predict_joint_selection_gaussian_parts(
    object = object, data = data, posterior_samples = posterior_samples,
    fixed_mu = means, within = scales[["tau_within"]], between = scales[["tau_between"]],
    draw_context = FALSE
  )
  execution_plan <- .data_selection_execution_plan(data)
  if (!is.null(integration_control)) {
    integration_control <- .check_selection_likelihood_control(integration_control)
    execution_plan <- .selection_joint_execution_plan_with_control(execution_plan, integration_control)
  }
  publication_groups <- model[["groups"]][["row_blocks"]]
  adjacency <- diag(TRUE, ncol(means))
  for (rows in c(parts[["full_dependency_blocks"]], publication_groups)) {
    adjacency[rows, rows] <- TRUE
  }
  row_blocks <- .known_v_block_indices(adjacency * 1)
  result <- .selection_sensitivity_run(
    means = parts[["means"]], covariance = parts[["covariance"]],
    context_covariance = parts[["context_covariance"]],
    selection_sei = .outcome_data_sei(object), selection_context = context,
    normalization_units = execution_plan[["row_blocks"]], row_blocks = row_blocks,
    execution_plan = execution_plan, latent_samples = latent_samples,
    random_covariance = parts[["random_covariance"]]
  )
  reference <- model
  for (mode in c("estimate_random_effects", "other_random_effects", "known_sampling_variance")) {
    reference[[mode]] <- "integrate"
    for (branch in reference[["active_branches"]]) {
      reference[["branches"]][[branch]][[mode]] <- "integrate"
    }
  }
  for (source in seq_along(reference[["sources"]][["random"]])) {
    reference[["sources"]][["random"]][[source]][["retained"]] <- FALSE
  }
  retained_random <- any(vapply(model[["sources"]][["random"]],
    function(source) isTRUE(source[["retained"]]), logical(1)))
  retained_sampling <- .selection_retains_sampling(data) &&
    .selection_sampling_structure(data)[["rank"]] > 0L
  reference_data <- data
  attr(reference_data, "selection_model") <- reference

  result[["posterior_samples"]] <- sample_info[["n_used"]]
  result[["total_posterior_samples"]] <- sample_info[["n_total"]]
  result[["latent_samples"]] <- latent_samples
  attr(result, "settings") <- list(
    max_posterior_samples = max_posterior_samples, latent_samples = latent_samples,
    seed = seed, integration_control = integration_control
  )
  attr(result, "target") <- list(
    comparison = "retained_context_reweighting",
    fitted = model, reference = reference,
    applicable = retained_random || retained_sampling,
    source_roles = list(random = model[["sources"]][["random"]],
      sampling = if (.is_data_known_v(data)) .known_v_selection_metadata(.data_known_v_data(data)) else NULL),
    sampling_structure = .selection_postfit_target_metadata(data)[["sampling_structure"]],
    reference_sampling_structure = .selection_postfit_target_metadata(reference_data)[["sampling_structure"]],
    publication_groups = model[["groups"]],
    likelihood_unit = "connected_source_and_publication_event",
    row_blocks = row_blocks,
    integration_control = execution_plan[c("points_per_scramble", "max_points_per_scramble",
      "scrambles", "seed", "relative_tolerance")]
  )
  class(result) <- c("selection_sensitivity_diagnostics", "data.frame")
  result
}


.selection_sensitivity_notify <- function(total_variation) {

  maximum_tv <- max(total_variation)
  tv_percent <- sprintf("%.1f", 100 * maximum_tv)
  if (maximum_tv >= 0.10) {
    warning(
      "Selection conditioning sensitivity was substantial: the largest unit-level ",
      "posterior-median total variation distance was ", tv_percent, "%. ",
      "Retaining and integrating the specified contexts define different reporting models. ",
      "Inspect 'selection_sensitivity_diagnostics()' for the fitted and reference specifications.",
      call. = FALSE
    )
  } else if (maximum_tv > 0.05) {
    message(
      "Selection conditioning sensitivity: the largest unit-level ",
      "posterior-median total variation distance was ", tv_percent, "%. ",
      "Inspect 'selection_sensitivity_diagnostics()' for unit-level ",
      "results and the fitted and reference specifications."
    )
  }
  invisible(NULL)
}


.selection_sensitivity_run <- function(
    means, covariance, context_covariance, selection_sei, selection_context,
    normalization_units, row_blocks, execution_plan, latent_samples,
    random_covariance = NULL) {

  means <- as.matrix(means)
  S <- nrow(means)
  K <- ncol(means)
  B <- length(row_blocks)
  if (!identical(.block_covariance_dim(covariance), c(S, K, K)) ||
      !identical(.block_covariance_dim(context_covariance), c(S, K, K)) ||
      length(selection_sei) != K || !B ||
      !identical(sort(as.integer(unlist(row_blocks, use.names = FALSE))), seq_len(K)) ||
      !identical(sort(as.integer(unlist(normalization_units, use.names = FALSE))), seq_len(K))) {
    stop("Selection sensitivity diagnostic inputs are inconsistent.", call. = FALSE)
  }
  metric_names <- c("ess_fraction", "ess_fraction_mcse", "total_variation",
    "total_variation_mcse", "log_weight_iqr")
  metrics <- stats::setNames(lapply(metric_names, function(x) matrix(NA_real_, S, B)), metric_names)
  integration_error <- matrix(0, S, B)
  integration_quadrature_error <- matrix(0, S, B)
  zero_context <- .block_covariance_zero_draws(context_covariance)
  for (draw in seq_len(S)) {
    if (isTRUE(selection_context[["use_normal"]][[draw]]) ||
        zero_context[[draw]]) {
      constant <- .selection_sensitivity_weight_metrics(rep(0, latent_samples))
      for (metric in metric_names) metrics[[metric]][draw, ] <- constant[[metric]]
      next
    }
    # The assembled K x K matrix, not the blocks: this spectral factorization
    # of the whole matrix and of its blocks differ in the last bits.
    latent <- .selection_sensitivity_covariance_draws(
      .block_covariance_dense(context_covariance, draw), latent_samples
    )
    conditional_means <- sweep(latent, 2L, means[draw, ], "+")
    draw_context <- BayesTools::selection_context_subset_rows(
      selection_context, rep.int(draw, latent_samples)
    )
    log_normalizers <- matrix(0, latent_samples, length(normalization_units))
    event_error <- numeric(length(normalization_units))
    event_quadrature_error <- numeric(length(normalization_units))
    for (event in seq_along(normalization_units)) {
      observations <- normalization_units[[event]]
      n <- length(observations)
      sigma <- .block_covariance_sub(covariance, draw, observations)
      rank_one <- if (!is.null(random_covariance)) {
        .selection_joint_declared_rank_one_loading(
          execution_plan[["sampling_factor_blocks"]][[event]],
          .block_covariance_sub(random_covariance, draw, observations))
      } else NULL
      lower <- matrix(sigma[lower.tri(sigma, diag = TRUE)], latent_samples,
        n * (n + 1L) / 2L, byrow = TRUE)
      event_context <- BayesTools::selection_context_subset_observations(draw_context, observations)
      mass <- .selection_joint_checked_event(
        compute = function(plan, rows) {
          if (is.null(rows)) rows <- seq_len(latent_samples)
          .selection_gaussian_event_mass(
            mean = conditional_means[rows, observations, drop = FALSE],
            covariance_lower = if (is.null(rank_one)) lower[rows, , drop = FALSE] else NULL,
            selection_se = selection_sei[observations],
            selection_context = BayesTools::selection_context_subset_rows(event_context, rows),
            execution_plan = plan,
            rank_one_loading = if (!is.null(rank_one)) matrix(rank_one,
              length(rows), n, byrow = TRUE) else NULL
          )
        },
        execution_plan = execution_plan,
        subject = "Selection sensitivity event integration",
        remedy = paste0("Increase 'max_points_per_scramble' in ",
          "'integration_control = set_selection_likelihood_control()'.")
      )
      if (any(!is.finite(mass[["log_mass"]]))) {
        stop("Selection sensitivity normalizers must be finite and positive.", call. = FALSE)
      }
      log_normalizers[, event] <- mass[["log_mass"]]
      event_error[[event]] <- max(mass[["relative_mcse"]])
      event_quadrature_error[[event]] <- max(mass[["relative_quadrature_error"]])
    }
    for (block in seq_along(row_blocks)) {
      events <- which(vapply(normalization_units, function(rows) {
        all(rows %in% row_blocks[[block]])
      }, logical(1)))
      block_metrics <- .selection_sensitivity_weight_metrics(
        rowSums(log_normalizers[, events, drop = FALSE])
      )
      for (metric in metric_names) metrics[[metric]][draw, block] <- block_metrics[[metric]]
      integration_error[draw, block] <- max(event_error[events])
      integration_quadrature_error[draw, block] <- max(event_quadrature_error[events])
    }
  }
  result <- .selection_sensitivity_summarize(metrics, row_blocks)
  result[["integration_relative_mcse_max"]] <- apply(integration_error, 2L, max)
  result[["integration_relative_quadrature_error_max"]] <-
    apply(integration_quadrature_error, 2L, max)
  result
}


.selection_sensitivity_covariance_draws <- function(
    covariance, latent_samples) {

  if (all(covariance == 0)) {
    return(matrix(0, nrow = latent_samples, ncol = nrow(covariance)))
  }
  factor <- .covariance_sampling_factor(
    .covariance_factorization(covariance)
  )
  if (is.null(factor)) {
    stop(
      "Selection sensitivity latent covariance is not positive semidefinite.",
      call. = FALSE
    )
  }
  matrix(
    stats::rnorm(latent_samples * nrow(covariance)),
    nrow = latent_samples,
    ncol = nrow(covariance)
  ) %*% factor
}


.selection_sensitivity_weight_metrics <- function(log_weights) {

  if (!is.numeric(log_weights) || length(log_weights) < 2L ||
      any(!is.finite(log_weights))) {
    stop(
      "Selection sensitivity log weights are invalid.",
      call. = FALSE
    )
  }
  weights <- exp(log_weights - max(log_weights))
  mean_weight <- mean(weights)
  mean_square <- mean(weights^2)
  normalized  <- weights / mean_weight

  ess_fraction <- mean_weight^2 / mean_square
  ess_influence <-
    (2 * mean_weight / mean_square) * (weights - mean_weight) -
    (mean_weight^2 / mean_square^2) * (weights^2 - mean_square)
  ess_mcse <- stats::sd(ess_influence) / sqrt(length(weights))

  absolute_difference <- abs(normalized - 1)
  total_variation     <- 0.5 * mean(absolute_difference)
  derivative <- -0.5 * mean(
    weights * sign(normalized - 1)
  ) / mean_weight^2
  tv_influence <- 0.5 *
    (absolute_difference - mean(absolute_difference)) +
    derivative * (weights - mean_weight)
  tv_mcse <- stats::sd(tv_influence) / sqrt(length(weights))

  list(
    ess_fraction        = ess_fraction,
    ess_fraction_mcse   = ess_mcse,
    total_variation     = total_variation,
    total_variation_mcse = tv_mcse,
    log_weight_iqr      = unname(stats::IQR(log_weights, type = 8))
  )
}


.selection_sensitivity_summarize <- function(metrics, row_blocks) {

  probabilities <- c(0.05, 0.5, 0.95)
  summarize <- function(values) {
    unname(stats::quantile(values, probabilities, names = FALSE, type = 8))
  }
  B <- length(row_blocks)
  out <- data.frame(
    block       = seq_len(B),
    rows        = vapply(
      row_blocks,
      function(rows) paste(rows, collapse = ","),
      character(1L)
    ),
    n_estimates = vapply(row_blocks, length, integer(1L)),
    stringsAsFactors = FALSE
  )
  for (metric in c("ess_fraction", "total_variation", "log_weight_iqr")) {
    values <- t(vapply(
      seq_len(B),
      function(block) summarize(metrics[[metric]][, block]),
      numeric(3L)
    ))
    out[[paste0(metric, "_q05")]]    <- values[, 1L]
    out[[paste0(metric, "_median")]] <- values[, 2L]
    out[[paste0(metric, "_q95")]]    <- values[, 3L]
  }
  for (metric in c("ess_fraction_mcse", "total_variation_mcse")) {
    out[[paste0(metric, "_median")]] <- apply(
      metrics[[metric]],
      2L,
      stats::median
    )
    out[[paste0(metric, "_max")]] <- apply(
      metrics[[metric]],
      2L,
      max
    )
  }
  out
}
