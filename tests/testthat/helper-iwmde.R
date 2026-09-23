# ============================================================================ #
# helper-iwmde.R
# ============================================================================ #

if (!exists("skip_if_missing_fits", mode = "function", inherits = FALSE)) {
  source(testthat::test_path("common-functions.R"), local = TRUE)
}

.bridge_shared_prior_semantics <- function(object, focal) {

  priors <- object[["priors"]]
  if (is.null(priors[["outcome"]][["mu"]]) &&
      !is.null(priors[["mods"]][["intercept"]])) {
    priors[["outcome"]][["mu"]] <- priors[["mods"]][["intercept"]]
  }
  priors[["mods"]][["intercept"]] <- NULL
  priors[["outcome"]][[focal]]     <- NULL
  priors[["mods"]][[focal]]        <- NULL
  priors <- lapply(priors, function(group) {
    if (is.null(group)) {
      return(NULL)
    }
    group <- lapply(group, function(prior) {
      list(
        distribution = prior[["distribution"]],
        parameters   = prior[["parameters"]],
        truncation   = prior[["truncation"]]
      )
    })
    group[sort(names(group))]
  })
  priors <- priors[vapply(priors, length, integer(1)) > 0L]

  return(priors[sort(names(priors))])
}


.expect_bridge_nesting <- function(null, full, focal) {

  expect_equal(null[["data"]][["outcome"]], full[["data"]][["outcome"]])
  null_data_attributes           <- attributes(null[["data"]])
  full_data_attributes           <- attributes(full[["data"]])
  null_data_attributes[["mods"]] <- NULL
  full_data_attributes[["mods"]] <- NULL
  expect_equal(null_data_attributes, full_data_attributes)
  expect_equal(null[["likelihood"]], full[["likelihood"]])
  expect_equal(
    .bridge_shared_prior_semantics(null, focal),
    .bridge_shared_prior_semantics(full, focal)
  )
}

.expect_iwmde_ok <- function(out, parameters, estimator = "q_grid_cmde",
                             mass_tolerance = NULL) {

  expect_named(out, parameters)

  for (parameter in parameters) {
    diagnostic <- out[[parameter]]

    expect_equal(diagnostic[["status"]], "ok")
    expect_length(diagnostic[["iwmde"]][["x"]], 20L)
    expect_length(diagnostic[["iwmde"]][["y"]], 20L)
    expect_true(all(is.finite(diagnostic[["iwmde"]][["y"]])))
    expect_true(all(diagnostic[["iwmde"]][["y"]] >= 0))
    expect_true(is.finite(diagnostic[["diagnostics"]][["integral"]]))
    expect_gt(diagnostic[["diagnostics"]][["integral"]], 0)
    expect_gt(diagnostic[["diagnostics"]][["display_mass_fraction"]], 0)
    expect_true(is.finite(diagnostic[["diagnostics"]][["display_mass_fraction"]]))
    expect_lt(diagnostic[["diagnostics"]][["display_mass_fraction"]], 10)
    expect_true(is.finite(diagnostic[["diagnostics"]][["min_ess"]]))
    expect_gt(diagnostic[["diagnostics"]][["min_ess"]], 0)
    expect_true(
      diagnostic[["diagnostics"]][["max_weight_share"]] >= 0 &&
        diagnostic[["diagnostics"]][["max_weight_share"]] <= 1
    )
    expect_true(is.finite(diagnostic[["diagnostics"]][["normalization_points"]]))
    expect_gt(diagnostic[["diagnostics"]][["normalization_points"]], 0)
    expect_true(is.finite(diagnostic[["diagnostics"]][["n_batches"]]))
    expect_equal(diagnostic[["diagnostics"]][["estimator"]], estimator)
    .expect_iwmde_mass_contract(
      diagnostic     = diagnostic,
      estimator      = estimator,
      mass_tolerance = mass_tolerance
    )
    expect_equal(
      diagnostic[["active_mass"]] +
        sum(diagnostic[["point_masses"]][["mass"]]),
      1,
      tolerance = 1e-8
    )
  }
}


.skip_if_missing_raw_fits <- function(names) {

  skip_if_missing_fits(names)
}


.mock_iwmde_empty_point_masses <- function() {

  data.frame(x = numeric(), mass = numeric())
}


.mock_iwmde_estimator_samples <- function(context, parameter,
                                          parameter_spec = NULL) {

  n <- 500L
  if (is.list(context) && !is.null(context[["posterior_samples"]])) {
    n_context <- nrow(context[["posterior_samples"]])
    if (is.finite(n_context) && n_context > 1L) {
      n <- n_context
    }
  }

  values <- try({
    spec <- .iwmde_parameter_spec(context, parameter, parameter_spec)
    if (!identical(spec[["status"]], "ok")) {
      stop("unsupported mock parameter", call. = FALSE)
    }
    values <- .iwmde_parameter_values(context, parameter, spec)
    condition_rows <- .iwmde_parameter_condition_rows(
      context        = context,
      parameter_spec = spec
    )
    if (length(condition_rows) == length(values)) {
      values <- values[condition_rows]
    }
    values
  }, silent = TRUE)

  if (!inherits(values, "try-error")) {
    values <- as.numeric(values)
    if (sum(is.finite(values)) >= 2L &&
        diff(range(values[is.finite(values)])) > sqrt(.Machine$double.eps)) {
      return(values)
    }
  }

  return(stats::qnorm(stats::ppoints(n)))
}


.mock_iwmde_density_grid <- function(samples, n_points) {

  n_points <- as.integer(n_points)[1]
  if (!is.finite(n_points) || n_points < 2L) {
    n_points <- 20L
  }

  finite <- as.numeric(samples)
  finite <- finite[is.finite(finite)]
  if (length(finite) < 2L ||
      diff(range(finite)) <= sqrt(.Machine$double.eps)) {
    finite <- stats::qnorm(stats::ppoints(max(20L, n_points)))
  }

  x <- seq(min(finite), max(finite), length.out = n_points)
  y <- stats::dnorm(x, mean = mean(finite), sd = max(stats::sd(finite), .1))

  return(list(x = x, y = y))
}


.mock_iwmde_good_diagnostics <- function(estimator, rows, value = NA_real_,
                                         include_bf = FALSE) {

  rows <- as.integer(rows)[1]
  if (!is.finite(rows) || rows < 500L) {
    rows <- 500L
  }
  ess <- max(100, rows / 2)

  list(
    integral                          = 1,
    plot_integral                     = 1,
    point_mass_total                  = 0,
    plot_total_mass                   = 1,
    display_mass_fraction             = 1,
    pilot_normalization_integral      = 1,
    final_normalization_integral      = 1,
    support_grid_normalization_integral = 1,
    normalization_relative_error      = 0,
    normalization_points              = 100L,
    normalization_range               = c(-1, 1),
    normalization_initial_points      = 50L,
    normalization_initial_range       = c(-1, 1),
    normalization_mass_ratio          = 1,
    max_ordinate_relative_change      = 0,
    max_normalizer_relative_change    = 0,
    p95_normalizer_relative_change    = 0,
    median_normalizer_relative_change = 0,
    max_quadrature_relative_change    = 0,
    normalization_refined_points      = 100L,
    normalization_refined_range       = c(-1, 1),
    n_refinement_steps                = 1L,
    active_mass                       = 1,
    n_candidate_rows                  = rows,
    n_denominator_rows                = rows,
    n_estimator_rows                  = rows,
    n_evaluated_rows                  = rows,
    weight_partitions                 = NULL,
    n_active                          = rows,
    n_active_state_keys               = 1L,
    min_active_state_rows             = rows,
    n_total                           = rows,
    min_finite_terms                  = rows,
    finite_terms                      = rows,
    min_ess                           = ess,
    ess                               = ess,
    max_weight_share                  = .01,
    max_mcse                          = .001,
    max_relative_mcse                 = .01,
    relative_mcse                     = .01,
    plot_scale_relative_mcse          = .001,
    bulk_probability_range            = c(.05, .95),
    bulk_x_range                      = c(-1, 1),
    bulk_max_relative_mcse            = .01,
    bulk_min_ess                      = ess,
    bulk_max_weight_share             = .01,
    tail_probabilities                = c(.05, .95),
    tail_target_x                     = c(-1, 1),
    tail_evaluation_x                 = c(-1, 1),
    tail_relative_mcse                = c(.01, .01),
    tail_ess                          = c(ess, ess),
    tail_max_weight_share             = c(.01, .01),
    plot_integral_mcse                = .001,
    plot_integral_relative_mcse       = .001,
    batch_size                        = max(1L, floor(rows / 20L)),
    n_batches                         = 20L,
    estimator                         = estimator,
    weight_method                     = if (identical(estimator, "q_grid_cmde")) {
      "conditional_grid"
    } else {
      "moment_match"
    },
    bf_included                       = isTRUE(include_bf),
    bf_value                          = value,
    bf_evaluation_value               = value,
    bf_ordinate                       = if (isTRUE(include_bf)) .4 else NA_real_,
    bf_mcse                           = .001,
    bf_relative_mcse                  = .01,
    bf_error_percent                  = 1,
    bf_finite_terms                   = rows,
    bf_ess                            = ess,
    bf_max_weight_share               = .01,
    bf_max_log_ratio                  = 0,
    bf_pilot_ordinate                 = .4,
    bf_validation_ordinate            = .4,
    bf_ordinate_relative_change       = 0,
    bf_ordinate_log_change            = 0,
    bf_pilot_ordinate_relative_change = 0,
    bf_pilot_ordinate_log_change      = 0
  )
}


.mock_iwmde_estimate_success <- function(omit_ordinate_levels = character()) {

  omit_ordinate_levels <- as.character(omit_ordinate_levels)
  force(omit_ordinate_levels)

  function(context, parameter, density_method, density_control,
           outputs = c("density", "ordinate"), values = NULL,
           parameter_spec = NULL, metadata = NULL, cache = NULL) {

    if (is.null(metadata)) {
      metadata <- list()
    }
    if (is.null(density_control)) {
      density_control <- list()
    }

    density_method <- .density_method_normalize_precomputed(density_method)
    estimator      <- .density_method_iwmde_estimator(density_method)
    samples        <- .mock_iwmde_estimator_samples(
      context        = context,
      parameter      = parameter,
      parameter_spec = parameter_spec
    )
    grid <- .mock_iwmde_density_grid(
      samples  = samples,
      n_points = density_control[["n_points"]]
    )
    rows <- max(500L, length(samples))
    base_diagnostic <- list(
      status       = "ok",
      samples      = samples,
      xlim         = range(grid[["x"]]),
      active_mass  = 1,
      point_masses = .mock_iwmde_empty_point_masses(),
      iwmde        = list(
        x         = grid[["x"]],
        y         = grid[["y"]],
        estimator = estimator
      )
    )

    out <- list(
      target             = list(parameter = parameter, parameter_spec = parameter_spec),
      diagnostics        = list(density = NULL, ordinate = NULL),
      posterior_density  = NULL,
      posterior_ordinate = NULL
    )

    if ("density" %in% outputs) {
      density_diagnostic <- base_diagnostic
      density_diagnostic[["diagnostics"]] <- .mock_iwmde_good_diagnostics(
        estimator = estimator,
        rows      = rows
      )
      out[["diagnostics"]][["density"]] <- density_diagnostic
      out[["posterior_density"]] <- .iwmde_posterior_density_attribute(
        diagnostic      = density_diagnostic,
        density_method  = density_method,
        metadata        = metadata,
        density_control = density_control
      )
    }

    if ("ordinate" %in% outputs) {
      value <- suppressWarnings(as.numeric(values))[1]
      level <- metadata[["level"]]
      include_bf <- is.finite(value) &&
        abs(value) <= 1 &&
        !identical(level %in% omit_ordinate_levels, TRUE)

      ordinate_diagnostic <- base_diagnostic
      ordinate_diagnostic[["iwmde"]] <- list(
        x         = value,
        y         = if (isTRUE(include_bf)) .4 else NA_real_,
        estimator = estimator
      )
      ordinate_diagnostic[["diagnostics"]] <- .mock_iwmde_good_diagnostics(
        estimator  = estimator,
        rows       = rows,
        value      = value,
        include_bf = include_bf
      )
      out[["diagnostics"]][["ordinate"]] <- ordinate_diagnostic
      out[["posterior_ordinate"]] <- .iwmde_posterior_ordinate_attribute(
        diagnostic      = ordinate_diagnostic,
        density_method  = density_method,
        metadata        = metadata,
        density_control = density_control
      )
    }

    class(out) <- c("iwmde_estimate", "list")
    return(out)
  }
}


.mock_marginal_means_ordinate_request_provenance <- function(
    object, parameter, null_hypothesis, density_method) {

  if (is.null(object)) {
    return(list())
  }

  specs <- tryCatch(
    .iwmde_marginal_means_specs(
      marginal_means_object = object,
      parameter             = parameter,
      type                  = "conditional",
      levels                = NULL
    ),
    error = function(e) list()
  )
  if (length(specs) == 0L) {
    return(list())
  }

  density_method <- .density_method_normalize_precomputed(density_method)
  density_control <- .marginal_means_iwmde_settings_control(
    object         = object,
    density_method = density_method,
    display_grid   = "ordinate"
  )
  method <- .density_method_iwmde_estimator(density_method)
  out    <- list()

  for (name in names(specs)) {
    spec     <- specs[[name]]
    metadata <- .marginal_means_iwmde_metadata(
      marginal_means_object = object,
      type                  = "conditional",
      spec                  = spec
    )
    out[[spec[["level"]]]] <- .iwmde_provenance_request(
      density_method  = density_method,
      method          = method,
      metadata        = metadata,
      density_control = density_control,
      value           = null_hypothesis,
      attribute       = "ordinate"
    )
  }

  return(out)
}


.local_mock_iwmde_estimate_success <- function(omit_ordinate_levels = character(),
                                               .env = parent.frame()) {

  testthat::local_mocked_bindings(
    .iwmde_estimate = .mock_iwmde_estimate_success(omit_ordinate_levels),
    .package        = "RoBMA",
    .env            = .env
  )
}


.local_mock_marginal_means_iwmde_success <- function(
    omit_ordinate_levels = character(), .env = parent.frame()) {

  testthat::local_mocked_bindings(
    .iwmde_estimate = .mock_iwmde_estimate_success(omit_ordinate_levels),
    .marginal_means_ordinate_request_provenance =
      .mock_marginal_means_ordinate_request_provenance,
    .package        = "RoBMA",
    .env            = .env
  )
}


.expect_iwmde_batch_equals_scalar <- function(fit_name, parameter,
                                             parameter_spec = NULL,
                                             n_rows = 3L,
                                             n_values = 3L) {

  .skip_if_missing_raw_fits(fit_name)

  context <- .iwmde_context(load_fit(fit_name, validate = FALSE))
  spec    <- .iwmde_parameter_spec(context, parameter, parameter_spec)
  values  <- .iwmde_parameter_values(context, parameter, spec)
  finite  <- is.finite(values)
  component <- .iwmde_parameter_components(context, parameter, spec)
  rows      <- which(component[["active"]] & finite)
  rows      <- head(rows, n_rows)

  expect_gt(length(rows), 0L)

  row_states <- .iwmde_row_states(context, rows, parameter, spec)
  keep       <- vapply(row_states, function(state) {
    is.finite(state[["baseline_log_q"]])
  }, logical(1))
  row_states <- row_states[keep]
  rows       <- rows[keep]

  expect_gt(length(rows), 0L)

  active_values <- values[component[["active"]] & finite]
  grid_values   <- as.numeric(stats::quantile(
    active_values,
    probs = seq(.25, .75, length.out = n_values),
    names = FALSE,
    type  = 8
  ))
  grid_values <- unique(grid_values[is.finite(grid_values)])

  expect_gt(length(grid_values), 0L)

  replacement <- .iwmde_replacement_spec(context, parameter, spec)
  batch       <- .iwmde_log_q_grid(
    context     = context,
    parameter   = parameter,
    values      = grid_values,
    row_states  = row_states,
    replacement = replacement
  )
  scalar      <- .iwmde_log_q_grid_scalar(
    context     = context,
    parameter   = parameter,
    values      = grid_values,
    row_states  = row_states,
    replacement = replacement
  )

  expect_equal(dim(batch), dim(scalar))
  expect_equal(is.finite(batch), is.finite(scalar))
  finite_terms <- is.finite(batch) & is.finite(scalar)
  expect_true(any(finite_terms))
  if (any(finite_terms)) {
    expect_equal(batch[finite_terms], scalar[finite_terms], tolerance = 1e-8)
  }
}


.expect_iwmde_predictor_fast_equals_scalar <- function(fit_name, parameter,
                                                       parameter_spec = NULL,
                                                       n_rows = 3L,
                                                       n_values = 3L) {

  .skip_if_missing_raw_fits(fit_name)

  context <- .iwmde_context(load_fit(fit_name, validate = FALSE))
  spec    <- .iwmde_parameter_spec(context, parameter, parameter_spec)
  values  <- .iwmde_parameter_values(context, parameter, spec)
  finite  <- is.finite(values)
  component <- .iwmde_parameter_components(context, parameter, spec)
  rows      <- which(component[["active"]] & finite)
  rows      <- head(rows, n_rows)

  expect_gt(length(rows), 0L)

  row_states <- .iwmde_row_states(context, rows, parameter, spec)
  keep       <- vapply(row_states, function(state) {
    is.finite(state[["baseline_log_q"]])
  }, logical(1))
  row_states <- row_states[keep]

  expect_gt(length(row_states), 0L)

  active_values <- values[component[["active"]] & finite]
  grid_values   <- as.numeric(stats::quantile(
    active_values,
    probs = seq(.25, .75, length.out = n_values),
    names = FALSE,
    type  = 8
  ))
  grid_values <- unique(grid_values[is.finite(grid_values)])

  expect_gt(length(grid_values), 0L)

  replacement <- .iwmde_replacement_spec(context, parameter, spec)
  fast        <- .iwmde_log_q_grid_predictor_batch(
    context     = context,
    parameter   = parameter,
    values      = grid_values,
    row_states  = row_states,
    replacement = replacement
  )
  scalar      <- .iwmde_log_q_grid_scalar(
    context     = context,
    parameter   = parameter,
    values      = grid_values,
    row_states  = row_states,
    replacement = replacement
  )

  expect_true(is.matrix(fast))
  expect_equal(dim(fast), dim(scalar))
  expect_equal(is.finite(fast), is.finite(scalar))
  finite_terms <- is.finite(fast) & is.finite(scalar)
  expect_true(any(finite_terms))
  if (any(finite_terms)) {
    expect_equal(fast[finite_terms], scalar[finite_terms], tolerance = 1e-8)
  }
}


.expect_iwmde_mass_contract <- function(diagnostic, estimator,
                                        mass_tolerance = NULL) {

  active_mass            <- diagnostic[["active_mass"]]
  normalization_integral <- if (identical(estimator, "q_grid_cmde")) {
    diagnostic[["diagnostics"]][["final_normalization_integral"]]
  } else {
    diagnostic[["diagnostics"]][["support_grid_normalization_integral"]]
  }
  normalization_relative_error <-
    diagnostic[["diagnostics"]][["normalization_relative_error"]]
  normalization_mass_ratio <- diagnostic[["diagnostics"]][["normalization_mass_ratio"]]

  expect_true(is.finite(active_mass))
  expect_gt(active_mass, 0)
  expect_true(is.finite(normalization_integral))
  expect_gt(normalization_integral, 0)
  expect_true(is.finite(normalization_mass_ratio))
  expect_gt(normalization_mass_ratio, 0)
  expect_true(is.finite(normalization_relative_error))
  expect_gte(normalization_relative_error, 0)
  if (identical(estimator, "q_grid_cmde")) {
    normalizer_change <- diagnostic[["diagnostics"]][["max_normalizer_relative_change"]]
    ordinate_change   <- diagnostic[["diagnostics"]][["max_ordinate_relative_change"]]
    expect_equal(normalization_mass_ratio, 1, tolerance = 1e-6)
    expect_true(is.finite(normalizer_change))
    expect_gte(normalizer_change, 0)
    expect_true(is.finite(ordinate_change))
    expect_gte(ordinate_change, 0)
  } else if (identical(estimator, "iwmde") && !is.null(mass_tolerance)) {
    expect_equal(
      normalization_mass_ratio,
      active_mass / normalization_integral,
      tolerance = 1e-6
    )
    expect_equal(normalization_integral, active_mass, tolerance = mass_tolerance)
  }
}


.expect_iwmde_normal_quadratic_equals_scalar <- function(fit_name, parameter,
                                                         parameter_spec = NULL,
                                                         n_rows = 3L,
                                                         n_values = 3L) {

  .skip_if_missing_raw_fits(fit_name)

  context <- .iwmde_context(load_fit(fit_name, validate = FALSE))
  spec    <- .iwmde_parameter_spec(context, parameter, parameter_spec)
  values  <- .iwmde_parameter_values(context, parameter, spec)
  finite  <- is.finite(values)
  component <- .iwmde_parameter_components(context, parameter, spec)
  rows      <- which(component[["active"]] & finite)
  rows      <- head(rows, n_rows)

  expect_gt(length(rows), 0L)

  row_states <- .iwmde_row_states(context, rows, parameter, spec)
  keep       <- vapply(row_states, function(state) {
    is.finite(state[["baseline_log_q"]])
  }, logical(1))
  row_states <- row_states[keep]

  expect_gt(length(row_states), 0L)

  active_values <- values[component[["active"]] & finite]
  grid_values   <- as.numeric(stats::quantile(
    active_values,
    probs = seq(.25, .75, length.out = n_values),
    names = FALSE,
    type  = 8
  ))
  grid_values <- unique(grid_values[is.finite(grid_values)])

  expect_gt(length(grid_values), 0L)

  replacement  <- .iwmde_replacement_spec(context, parameter, spec)
  active_setup <- row_states[[1L]][["active_setup"]]
  unit         <- if (.is_data_multilevel(context[["data"]])) {
    "cluster"
  } else {
    "estimate"
  }
  setup <- .iwmde_predictor_setup(
    context      = context,
    row_states   = row_states,
    active_setup = active_setup,
    unit         = unit
  )
  basis <- .iwmde_predictor_update_basis(
    context     = context,
    parameter   = parameter,
    row_states  = row_states,
    replacement = replacement,
    setup       = setup
  )
  fast <- .iwmde_log_q_grid_normal_location_group(
    context     = context,
    parameter   = parameter,
    values      = grid_values,
    row_states  = row_states,
    replacement = replacement,
    setup       = setup,
    basis       = basis
  )
  scalar <- .iwmde_log_q_grid_scalar(
    context     = context,
    parameter   = parameter,
    values      = grid_values,
    row_states  = row_states,
    replacement = replacement
  )

  expect_true(is.matrix(fast))
  expect_equal(dim(fast), dim(scalar))
  expect_equal(is.finite(fast), is.finite(scalar))
  finite_terms <- is.finite(fast) & is.finite(scalar)
  expect_true(any(finite_terms))
  if (any(finite_terms)) {
    expect_equal(fast[finite_terms], scalar[finite_terms], tolerance = 1e-8)
  }
}
