# Exact Gaussian retained-intercept conditioning for qCMDE. This changes the
# conditioning coordinates, not the marginal target or its scalar prior.

.iwmde_retained_location_dynamic_prior <- function(value) {

  if (is.expression(value) || is.call(value) || is.function(value) || is.environment(value)) {
    return(TRUE)
  }
  if (is.list(value)) {
    return(any(vapply(value, .iwmde_retained_location_dynamic_prior, logical(1L))))
  }
  FALSE
}

.iwmde_retained_location_plan <- function(context, parameter, parameter_spec) {

  data <- context[["data"]]
  fit <- context[["object"]][["fit"]]
  if (!.is_data_joint_selection(data) || !.is_data_random(data) ||
      .data_outcome_type(data) != "norm" || .selection_retains_sampling(data) ||
      .is_data_scale(data) || .is_data_weights(data) ||
      !identical(parameter_spec[["type"]], "primitive") ||
      !inherits(fit, "BayesTools_fit") || length(context[["indicator_names"]]) ||
      .iwmde_retained_location_dynamic_prior(context[["flat_prior_list"]]) ||
      any(nzchar(.random_allocation_inclusion_indicators(context[["flat_prior_list"]]))) ||
      .iwmde_parameter_controls_sampled_random_sd(context, parameter)) {
    return(NULL)
  }
  coordinates <- BayesTools::parameter_coordinates(fit)
  coordinate <- coordinates[coordinates[["coordinate_name"]] == parameter, , drop = FALSE]
  if (nrow(coordinate) != 1L || coordinate[["role"]] != "fixed_coefficient" ||
      coordinate[["internal"]] || coordinate[["formula_parameter"]] != "mu" ||
      coordinate[["term"]] != "intercept") {
    return(NULL)
  }
  dependencies <- BayesTools::JAGS_formula_coordinate_dependencies(fit, parameter)
  if (any(dependencies[["formula_parameter"]] != "mu") ||
      any(dependencies[["dependency_type"]] != "coefficient")) {
    return(NULL)
  }
  samples <- context[["posterior_samples"]]
  if (!parameter %in% colnames(samples) || nrow(samples) == 0L) return(NULL)
  prior <- .iwmde_focal_prior(context, parameter, samples[1L, ])
  if (is.null(prior) || !BayesTools::is.prior.simple(prior) ||
      BayesTools::is.prior.factor(prior) || BayesTools::is.prior.point(prior) ||
      BayesTools::is.prior.discrete(prior) ||
      !.iwmde_can_use_focal_prior_delta(prior) ||
      !is.null(attr(prior, "multiply_by", exact = TRUE))) {
    return(NULL)
  }
  # The shared accessor certifies affine dependence, including rejection of
  # logged intercepts, expressions and random-scale dependencies.
  basis <- BayesTools::JAGS_formula_predictor_basis(fit,
    directions = stats::setNames(1, parameter), posterior_samples = samples[1L, , drop = FALSE])
  K <- nrow(data[["outcome"]])
  if (!identical(basis[["status"]], "affine") ||
      !identical(basis[["parameter"]], "mu") ||
      !identical(dim(basis[["basis"]]), c(1L, K)) ||
      any(basis[["basis"]] != 1)) {
    return(NULL)
  }
  model <- .data_selection_model(data)
  sources <- model[["sources"]][["random"]]
  retained <- vapply(sources, function(source) isTRUE(source[["retained"]]), logical(1L))
  retained_names <- vapply(sources[retained], `[[`, character(1L), "name")
  design <- .fitted_formula_design(context[["object"]], "mu", required = TRUE)
  terms <- design[["random_effects"]]
  eligible <- vapply(terms, function(term) {
    term[["block_name"]] %in% retained_names &&
      identical(.random_effect_term_compile_mode(term), "sampled") &&
      term[["structure"]] %in% c("id", "diag", "us") &&
      .brma_mv_random_term_is_pure_intercept(term) &&
      is.matrix(term[["model_matrix"]]) &&
      identical(dim(term[["model_matrix"]]), c(K, 1L)) &&
      all(term[["model_matrix"]] == 1) && length(term[["sd_parameter_names"]]) == 1L
  }, logical(1L))
  if (sum(eligible) != 1L) return(NULL)
  term <- terms[[which(eligible)]]
  groups <- term[["group_map"]]
  J <- term[["n_groups"]]
  if (!is.numeric(groups) || length(groups) != K || anyNA(groups) ||
      any(groups != as.integer(groups)) || !is.numeric(J) || length(J) != 1L ||
      !is.finite(J) || J < 1L || J != as.integer(J) ||
      !identical(sort(unique(as.integer(groups))), seq_len(J))) {
    return(NULL)
  }
  list(parameter = parameter, block = term[["block_name"]], term = term,
    group_map = as.integer(groups), group_rows = match(seq_len(J), groups),
    n_groups = as.integer(J), prior = prior)
}

.iwmde_retained_location_statistics <- function(coefficients, sd) {

  if (!is.matrix(coefficients) || !is.numeric(coefficients) ||
      !nrow(coefficients) || !ncol(coefficients) || any(!is.finite(coefficients)) ||
      !is.numeric(sd) || length(sd) != nrow(coefficients) ||
      any(!is.finite(sd)) || any(sd < 0)) {
    stop("Retained-location Gaussian conditional inputs are invalid.", call. = FALSE)
  }
  if (any(sd == 0)) return(NULL)
  S <- nrow(coefficients)
  J <- ncol(coefficients)
  change <- .iwmde_normal_location_likelihood_change_estimate(
    yi = numeric(J), mu = -coefficients, mu_basis = matrix(1, S, J),
    tau_within = matrix(sd, S, J), sei = numeric(J), data_weights = NULL)
  if (is.null(change)) return(NULL)
  if (any(!is.finite(change[["linear"]])) ||
      any(!is.finite(change[["quadratic"]])) || any(change[["quadratic"]] <= 0)) {
    stop("Retained-location Gaussian conditional statistics are invalid.", call. = FALSE)
  }
  change
}

.iwmde_retained_location_row_states <- function(context, plan, rows) {

  if (!identical(plan[["method"]], "q_grid_cmde") || !length(rows)) return(NULL)
  parameter <- plan[["target"]][["parameter"]]
  parameter_spec <- plan[["execution_spec"]]
  if (!is.null(parameter_spec[["direction"]])) return(NULL)
  coordinate <- parameter
  coordinate_spec <- parameter_spec
  if (identical(parameter_spec[["type"]], "linear")) {
    weights <- .iwmde_linear_weights(parameter_spec[["weights"]])
    if (length(weights) != 1L || !identical(unname(weights), 1)) return(NULL)
    coordinate <- names(weights)[[1L]]
    coordinate_spec <- list(type = "primitive", parameter = coordinate)
  }
  # An exact unit-weight alias has the same retained-location conditional law.
  # The original target label/spec remain authoritative for output and fallback.
  metadata <- .iwmde_retained_location_plan(context, coordinate, coordinate_spec)
  if (is.null(metadata)) return(NULL)
  data <- context[["data"]]
  samples <- context[["posterior_samples"]][rows, , drop = FALSE]
  S <- length(rows)
  K <- nrow(data[["outcome"]])
  # Fixed, target-independent sampling information. A singular sampling law
  # has no such finite bound in this first route and keeps ordinary qCMDE.
  sampling <- .selection_joint_sampling_block(
    .data_selection_execution_plan(data)[["sampling"]], seq_len(K))
  root <- tryCatch(chol(sampling), error = function(e) NULL)
  if (is.null(root)) return(NULL)
  sampling_information <- sum(forwardsolve(t(root), rep(1, K))^2)
  information_cap <- S * sampling_information
  if (!is.finite(information_cap) || information_cap <= 0) return(NULL)
  sd <- .random_effect_term_sd_samples(metadata[["term"]], samples, K)
  sd <- .expand_brma_mv_heterogeneity_samples(sd, S, K)
  if (any(!is.finite(sd)) || any(sd < 0)) {
    stop("Retained-location prior standard deviations are invalid.", call. = FALSE)
  }
  if (any(sd != matrix(sd[, 1L], S, K))) return(NULL)
  sd <- sd[, 1L]
  information <- (sqrt(metadata[["n_groups"]]) / sd)^2
  use_transformed <- sd > 0 & is.finite(information) & information <= information_cap
  if (!any(use_transformed)) return(NULL)
  selected <- which(use_transformed)
  selected_samples <- samples[selected, , drop = FALSE]
  # Conditional formula prediction materializes actual sampled coefficients
  # (including noncentered/mean-centered reconstructions), never BLUP means.
  contributions <- .evaluate.brma.random_effects(
    fit = context[["object"]][["fit"]], data = data, priors = context[["priors"]],
    posterior_samples = selected_samples, same_data = TRUE, required = TRUE,
    formula_target = "conditional", blocks = metadata[["block"]], object = context[["object"]])
  if (!is.matrix(contributions) ||
      !identical(dim(contributions), c(length(selected), K)) ||
      any(!is.finite(contributions))) {
    stop("Retained-location sampled coefficients are invalid.", call. = FALSE)
  }
  coefficients <- contributions[, metadata[["group_rows"]], drop = FALSE]
  if (any(contributions != coefficients[, metadata[["group_map"]], drop = FALSE])) {
    stop("Retained-location sampled coefficients disagree within groups.", call. = FALSE)
  }
  change <- .iwmde_retained_location_statistics(coefficients, sd[selected])
  if (is.null(change)) return(NULL)
  focal_prior <- .iwmde_focal_log_prior_values(metadata[["prior"]],
    selected_samples[, coordinate], coordinate)
  if (any(!is.finite(focal_prior))) {
    .iwmde_stop_construction_failure("q_grid_cmde", parameter,
      rows[selected[!is.finite(focal_prior)]],
      stage = "retained-location conditional prior evaluation",
      detail = "a retained conditional prior value was not finite")
  }
  transform <- list(type = "retained_gaussian_intercept", parameter = parameter,
    block = metadata[["block"]], n_groups = metadata[["n_groups"]])
  states <- vector("list", S)
  for (i in seq_along(selected)) {
    position <- selected[[i]]
    states[[position]] <- .iwmde_new_row_state(list(
      row_index = rows[[position]], row = selected_samples[i, ], active_key = "all",
      current = selected_samples[i, coordinate],
      baseline_log_q = focal_prior[[i]], baseline_log_lik = 0,
      baseline_log_prior = focal_prior[[i]], baseline_focal_log_prior = focal_prior[[i]],
      focal_prior = metadata[["prior"]], use_focal_prior_delta = TRUE,
      likelihood_mode = "conditional", state_scope = "local",
      gaussian_change = list(linear = change[["linear"]][[i]], quadratic = change[["quadratic"]][[i]]),
      conditioning_transform = transform))
  }
  ordinary <- which(!use_transformed)
  if (length(ordinary)) {
    ordinary_states <- .iwmde_row_states_grouped_marginal(context, rows[ordinary],
      parameter, parameter_spec, estimator = plan[["method"]])
    if (is.null(ordinary_states)) {
      ordinary_states <- .iwmde_row_states(context, rows[ordinary], parameter,
        parameter_spec, estimator = plan[["method"]])
    }
    states[ordinary] <- ordinary_states
  }
  attr(states, "conditioning_policy") <- c(transform, list(
    information_multiplier = S, sampling_information = sampling_information,
    information_cap = information_cap, transformed_rows = length(selected),
    ordinary_rows = length(ordinary),
    classification_hash = .iwmde_hash("retained_location_charts", list(rows, use_transformed))))
  states
}

.iwmde_log_q_grid_retained_location <- function(context, parameter, values,
                                                 row_states, replacement) {

  transforms <- lapply(row_states, `[[`, "conditioning_transform")
  if (!length(transforms) || !all(vapply(transforms, function(transform) {
    is.list(transform) && identical(transform[["type"]], "retained_gaussian_intercept") &&
      identical(transform[["parameter"]], parameter)
  }, logical(1L)))) {
    stop("Retained-location conditional row metadata are inconsistent.", call. = FALSE)
  }
  current <- vapply(row_states, `[[`, numeric(1L), "current")
  change <- list(
    linear = vapply(row_states, function(state) state[["gaussian_change"]][["linear"]], numeric(1L)),
    quadratic = vapply(row_states, function(state) state[["gaussian_change"]][["quadratic"]], numeric(1L)))
  prior <- row_states[[1L]][["focal_prior"]]
  if (!all(vapply(row_states, function(state) identical(state[["focal_prior"]], prior), logical(1L)))) {
    stop("Retained-location conditional rows have different focal priors.", call. = FALSE)
  }
  # Constants from other nuisance priors cancel in each conditional. Evaluate
  # the original focal prior directly, avoiding large add/subtract constants.
  log_prior <- rep(.iwmde_focal_log_prior_values(prior, values, parameter), length(row_states))
  .iwmde_normal_location_log_q_grid(values, current, numeric(length(row_states)),
    change, log_prior)
}

.iwmde_retained_location_normalizer <- function(state) {

  prior <- state[["focal_prior"]]
  current <- state[["current"]]
  linear <- state[["gaussian_change"]][["linear"]]
  precision <- state[["gaussian_change"]][["quadratic"]]
  if (!is.finite(current) || !is.finite(linear) || !is.finite(precision) || precision <= 0) {
    stop("Retained-location normalizer inputs are invalid.", call. = FALSE)
  }
  scale <- 1 / sqrt(precision)
  center <- current + linear / precision
  if (!is.finite(center) || !is.finite(scale) || scale <= 0) {
    stop("Retained-location Gaussian kernel cannot be represented.", call. = FALSE)
  }
  support <- .iwmde_prior_support(prior)
  if (identical(prior[["distribution"]], "normal")) {
    prior_mean <- prior[["parameters"]][["mean"]]
    prior_sd <- prior[["parameters"]][["sd"]]
    total_sd <- .root_sum_squares(prior_sd, scale)
    small <- min(prior_sd, scale)
    large <- max(prior_sd, scale)
    posterior_sd <- small / sqrt(1 + (small / large)^2)
    posterior_mean <- prior_mean * (scale / total_sd)^2 + center * (prior_sd / total_sd)^2
    log_mass <- .iwmde_log_normal_interval_prob(support[1L], support[2L],
      posterior_mean, posterior_sd)
    log_normalizer <- .iwmde_focal_log_prior(prior, current) -
      stats::dnorm(current, posterior_mean, posterior_sd, log = TRUE) + log_mass
    if (!is.finite(log_normalizer)) {
      stop("Retained-location normal-prior normalization was not finite.", call. = FALSE)
    }
    return(list(log_normalizer = log_normalizer, relative_error = 0,
      method = "normal_product", posterior_mean = posterior_mean, posterior_sd = posterior_sd))
  }

  # Standardize to this row's Gaussian width. Prior quantiles split narrow
  # prior features; support endpoints retain their original exact meaning.
  prior_points <- as.numeric(BayesTools::quant(prior,
    c(.001, .01, .1, .25, .5, .75, .9, .99, .999)))
  prior_points <- prior_points[is.finite(prior_points)]
  support_z <- (support - center) / scale
  anchors <- sort(unique(c(0, (current - center) / scale,
    (prior_points - center) / scale)))
  anchors <- anchors[is.finite(anchors)]
  target <- .iwmde_qcmde_refinement_target()
  tolerance <- min(target / 4, .Machine$double.eps^.25)
  radius <- max(4, if (support_z[1L] > 0) support_z[1L] + 4 else 0,
    if (support_z[2L] < 0) -support_z[2L] + 4 else 0)
  log_kernel <- function(z) {
    .iwmde_focal_log_prior_values(prior, center + scale * z) + log(scale) - .5 * z^2
  }
  for (iteration in seq_len(10L)) {
    lower <- max(-radius, support_z[1L])
    upper <- min(radius, support_z[2L])
    if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
      stop("Retained-location prior support cannot be integrated on the Gaussian scale.", call. = FALSE)
    }
    cuts <- sort(unique(c(lower, upper, anchors[anchors > lower & anchors < upper])))
    log_values <- log_errors <- rep(-Inf, length(cuts) - 1L)
    for (part in seq_len(length(cuts) - 1L)) {
      left <- cuts[[part]]
      right <- cuts[[part + 1L]]
      probes <- c(left, left + (right - left) * c(.25, .5, .75), right)
      heights <- log_kernel(probes)
      shift <- max(heights[is.finite(heights)], -Inf)
      if (!is.finite(shift)) {
        stop("Retained-location integration could not find a finite positive interior density.", call. = FALSE)
      }
      integral <- stats::integrate(function(z) exp(log_kernel(z) - shift),
        left, right, rel.tol = tolerance, abs.tol = 0, stop.on.error = TRUE)
      if (!is.finite(integral$value) || integral$value <= 0 ||
          !is.finite(integral$abs.error) || integral$abs.error < 0) {
        stop("Retained-location conditional integration returned invalid mass or error.", call. = FALSE)
      }
      log_values[[part]] <- log(integral$value) + shift
      if (integral$abs.error > 0) log_errors[[part]] <- log(integral$abs.error) + shift
    }
    log_mass <- .log_sum_exp(log_values)
    log_error <- if (any(is.finite(log_errors))) .log_sum_exp(log_errors) else -Inf
    # For any proper prior, the omitted Gaussian-weighted prior integral is
    # at most exp(-radius^2/2). Keep this bound in log scale, not underflowed 0.
    log_tail <- if (lower == support_z[1L] && upper == support_z[2L]) -Inf else -.5 * radius^2
    log_error_total <- if (is.finite(log_error) || is.finite(log_tail)) {
      .log_sum_exp(c(log_error, log_tail))
    } else {
      -Inf
    }
    log_lower_mass <- .iwmde_logspace_subtract(log_mass, log_error)
    relative_error <- if (is.finite(log_lower_mass)) exp(log_error_total - log_lower_mass) else Inf
    if (is.finite(relative_error) && relative_error <= target) {
      log_normalizer <- .5 * (linear / sqrt(precision))^2 + log_mass
      if (!is.finite(log_normalizer)) {
        stop("Retained-location conditional normalizer was not finite.", call. = FALSE)
      }
      return(list(log_normalizer = log_normalizer, relative_error = relative_error,
        method = "scaled_gaussian_integral"))
    }
    radius <- max(2 * radius, sqrt(max(0, -2 * (log_mass + log(target / 4)))))
  }
  stop("Retained-location conditional normalization was rejected by diagnostics: relative integration error was ",
    format(relative_error, digits = 4), ". Inspect the scalar prior and conditional-kernel scale.", call. = FALSE)
}
