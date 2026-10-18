# ============================================================================ #
# test-02-iwmde-normal-predictor-grid.R
# ============================================================================ #

context("IWMDE normal predictor grid")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-iwmde.R"))
source(testthat::test_path("helper-selection-kernel.R"))


# Cached fixtures that could reach the native normal candidate grid: ordinary
# normal models without a selection model. The guard inside
# .iwmde_predictor_normal_grid_log_lik() makes the final decision.
.normal_grid_fit_names <- function() {

  fit_names <- unique(c(
    list_fits(class = c("brma", "BMA", "RoBMA", "NoBMA")),
    list_fits(family = "norm")
  ))
  available <- vapply(fit_names, function(fit_name) {
    object <- load_fit(fit_name, validate = FALSE)
    .iwmde_capability(object = object)[["available"]]
  }, logical(1))
  return(fit_names[available])
}


test_that("normal-grid discovery preserves the declared component-scale refusal", {

  fit_name <- "hypothesis_plan_two_scale"
  skip_if_missing_fits(fit_name)
  object <- load_fit(fit_name, validate = FALSE)
  capability <- .iwmde_capability(object = object)

  expect_false(capability[["available"]])
  expect_identical(capability[["class"]], c(
    "RoBMA_density_method_scale_components",
    "RoBMA_density_method_unavailable"
  ))
  expect_match(capability[["reason"]], "several scale formulas", fixed = TRUE)
  expect_error(.iwmde_context(object), class = "RoBMA_density_method_scale_components")
  expect_false(fit_name %in% .normal_grid_fit_names())
})


# One location target, one moderator slope, the heterogeneity scale, and a
# scale-regression coefficient where the fixture carries one.
.normal_grid_parameters <- function(context) {

  columns    <- colnames(context[["posterior_samples"]])
  location   <- intersect(c("mu", "mu_intercept"), columns)
  moderators <- setdiff(grep("^mu_", columns, value = TRUE), "mu_intercept")
  scale      <- intersect(c("tau", "log_tau_intercept"), columns)
  scale_mods <- setdiff(grep("^log_tau_", columns, value = TRUE), "log_tau_intercept")

  return(unique(c(
    utils::head(location, 1L),
    utils::head(moderators, 1L),
    scale,
    utils::head(scale_mods, 1L)
  )))
}


.normal_grid_inputs <- function(context, parameter, n_rows = 6L, n_values = 5L) {

  # These structural exclusions match the native normal-grid guard. Their
  # separate likelihoods do not reach this kernel, so do not construct row
  # states for them merely to discover the same ineligibility afterward.
  data <- context[["data"]]
  if (!identical(.data_outcome_type(data), "norm") ||
      .is_data_joint_selection(data) || .is_data_known_v(data) ||
      .is_data_random(data) || .is_data_multilevel(data)) {
    return(NULL)
  }
  spec <- tryCatch(.iwmde_parameter_spec(context, parameter, NULL),
    BayesTools_ordered_expression_unavailable = function(e) NULL)
  if (is.null(spec) || !identical(spec[["status"]], "ok")) {
    return(NULL)
  }
  if (spec[["type"]] %in% c("ordered", "ordered_source")) {
    return(NULL)
  }
  values <- .iwmde_parameter_values(context, parameter, spec)
  component <- .iwmde_parameter_components(context, parameter, spec)
  if (is.null(values) || is.null(component)) {
    return(NULL)
  }
  active <- component[["active"]] & is.finite(values)
  if (sum(active) < 2L) {
    return(NULL)
  }
  row_states <- tryCatch(.iwmde_row_states(
    context, utils::head(which(active), n_rows), parameter, spec
  ), BayesTools_ordered_expression_unavailable = function(e) NULL)
  if (is.null(row_states)) {
    return(NULL)
  }
  row_states <- row_states[vapply(row_states, function(state) {
    is.finite(state[["baseline_log_q"]])
  }, logical(1))]
  if (length(row_states) == 0L) {
    return(NULL)
  }
  grid_values <- as.numeric(stats::quantile(
    values[active],
    probs = seq(.1, .9, length.out = n_values),
    names = FALSE,
    type  = 8
  ))

  return(list(
    row_states  = row_states,
    values      = grid_values[is.finite(grid_values)],
    replacement = .iwmde_replacement_spec(context, parameter, spec)
  ))
}


test_that("normal-grid fixtures apply structural eligibility before row preparation", {

  testthat::local_mocked_bindings(
    .iwmde_parameter_spec = function(...) stop("eligible row preparation reached"),
    .package = "RoBMA"
  )
  normal <- structure(list(), outcome_type = "norm")
  for (flag in c("known_V", "random", "cluster")) {
    data <- normal
    attr(data, flag) <- TRUE
    expect_null(.normal_grid_inputs(list(data = data), "mu"))
  }
  selection <- normal
  attr(selection, "selection_model") <- structure(list(schema_version = 3L),
    class = "RoBMA_selection_model")
  expect_null(.normal_grid_inputs(list(data = selection), "mu"))
  expect_null(.normal_grid_inputs(list(data = structure(list(), outcome_type = "bin")), "mu"))
  expect_error(.normal_grid_inputs(list(data = normal), "mu"),
    "eligible row preparation reached", fixed = TRUE)
})


# The two routes on one batch: the native grid where its guard admits the
# batch, and the candidate route with the native route mocked away. 'served'
# counts the state groups the native route answered and 'groups' how many were
# asked, so a batch whose branches split between the two routes is visible.
.normal_grid_both_routes <- function(context, parameter, inputs) {

  served <- 0L
  groups <- 0L
  native <- local({
    original <- .iwmde_predictor_normal_grid_log_lik
    testthat::local_mocked_bindings(
      .iwmde_predictor_normal_grid_log_lik = function(...) {
        value <- original(...)
        groups <<- groups + 1L
        served <<- served + as.integer(!is.null(value))
        value
      },
      .package = "RoBMA"
    )
    .iwmde_log_q_grid_predictor_batch(
      context     = context,
      parameter   = parameter,
      values      = inputs[["values"]],
      row_states  = inputs[["row_states"]],
      replacement = inputs[["replacement"]]
    )
  })
  candidates <- local({
    testthat::local_mocked_bindings(
      .iwmde_predictor_normal_grid_log_lik = function(...) NULL,
      .package = "RoBMA"
    )
    .iwmde_log_q_grid_predictor_batch(
      context     = context,
      parameter   = parameter,
      values      = inputs[["values"]],
      row_states  = inputs[["row_states"]],
      replacement = inputs[["replacement"]]
    )
  })

  return(list(native = native, candidates = candidates,
              served = served, groups = groups))
}


.expect_normal_grid_agreement <- function(routes, label) {

  expect_true(is.matrix(routes[["candidates"]]), info = label)
  expect_identical(dim(routes[["native"]]), dim(routes[["candidates"]]), info = label)
  expect_identical(routes[["native"]] == -Inf, routes[["candidates"]] == -Inf,
                   info = paste0(label, ": -Inf pattern"))
  expect_identical(is.finite(routes[["native"]]), is.finite(routes[["candidates"]]),
                   info = label)
  finite <- is.finite(routes[["candidates"]])
  expect_true(any(finite), info = paste0(label, ": no finite log density"))
  deviation <- max(abs(routes[["native"]][finite] - routes[["candidates"]][finite]) /
                     pmax(abs(routes[["candidates"]][finite]), 1))
  expect_lt(deviation, 1e-12, label = paste0(label, ": max scaled deviation"))

  return(invisible(NULL))
}


test_that("dependent Student-t focal coordinates retain the joint prior on both row routes", {

  fit_name <- "hypothesis_plan_meandif_mt"
  skip_if_missing_fits(fit_name)
  object <- load_fit(fit_name)
  context <- .iwmde_context(object)
  coordinates <- c("mu_g1[1]", "mu_g1[2]")
  prior <- attr(object[["fit"]], "prior_list")[["mu_g1"]]
  expect_identical(prior[["distribution"]], "mt")
  expect_true(is.na(prior[["parameters"]][["K"]]))
  expect_identical(ncol(attr(prior, "factor_design")), 2L)

  for (parameter in coordinates) {
    spec <- .iwmde_parameter_spec(context, parameter, NULL)
    nuisance <- setdiff(coordinates, parameter)
    nuisance_values <- context[["posterior_samples"]][, nuisance]
    rows <- unique(c(which.min(abs(nuisance_values)), which.max(abs(nuisance_values))))
    expect_length(rows, 2L)
    scalar <- .iwmde_row_states(context, rows, parameter, spec)
    grouped <- .iwmde_row_states_grouped_marginal(context, rows, parameter, spec)
    expect_length(grouped, 2L)
    for (states in list(scalar, grouped)) {
      for (state in states) {
        expect_true(is.finite(state[["baseline_log_q"]]))
        expect_identical(state[["focal_prior"]], prior)
        expect_false(state[["use_focal_prior_delta"]])
        expect_identical(state[["baseline_focal_log_prior"]], NA_real_)
      }
    }
    expect_equal(vapply(scalar, `[[`, numeric(1), "baseline_log_q"),
      vapply(grouped, `[[`, numeric(1), "baseline_log_q"), tolerance = 1e-12)
    values <- c(-.6, 0, .4, 1)
    replacement <- .iwmde_replacement_spec(context, parameter, spec)
    routes <- lapply(list(scalar, grouped), function(states) {
      inputs <- list(values = values, row_states = states, replacement = replacement)
      production <- .iwmde_log_q_grid_predictor_batch(context, parameter, values, states, replacement)
      # The ordinary location optimization answers before the native-grid
      # dispatcher. Exercise its two existing alternatives on the same
      # declared affine basis without changing either guard.
      both <- local({
        testthat::local_mocked_bindings(
          .iwmde_log_q_grid_normal_location_group = function(...) NULL,
          .package = "RoBMA"
        )
        .normal_grid_both_routes(context, parameter, inputs)
      })
      expect_gt(both[["served"]], 0L)
      .expect_normal_grid_agreement(both, parameter)
      direct <- vapply(states, function(state) {
        vapply(values, function(value) .iwmde_log_q_replacement(
          context, state, parameter, value, replacement), numeric(1))
      }, numeric(length(values)))
      expect_identical(is.finite(both[["native"]]), is.finite(direct))
      expect_identical(both[["native"]] == -Inf, direct == -Inf)
      finite <- is.finite(direct)
      expect_lt(max(abs(both[["native"]][finite] - direct[finite])), 1e-10)
      expect_identical(is.finite(production), is.finite(direct))
      expect_lt(max(abs(production[finite] - direct[finite])), 1e-10)
      expect_identical(is.finite(production), is.finite(both[["candidates"]]))
      expect_lt(max(abs(production[finite] - both[["candidates"]][finite])), 1e-10)
      both[["native"]]
    })
    expect_identical(is.finite(routes[[1L]]), is.finite(routes[[2L]]))
    finite <- is.finite(routes[[1L]])
    expect_lt(max(abs(routes[[1L]][finite] - routes[[2L]][finite])), 1e-10)

    # A vector t has dependent coordinates even with diagonal scale. The
    # conditional ratio at each actual nuisance row must use the joint law.
    joint_ratios <- marginal_ratios <- numeric(length(scalar))
    for (i in seq_along(scalar)) {
      state <- scalar[[i]]
      before <- matrix(state[["row"]], nrow = 1L,
        dimnames = list(NULL, names(state[["row"]])))
      after <- before
      after[, parameter] <- 1
      direct_joint <- function(samples) mvtnorm::dmvt(samples[, coordinates, drop = FALSE],
        delta = rep(0, 2L), sigma = diag(.5^2, 2L), df = 3, type = "shifted", log = TRUE)
      joint_ratios[[i]] <- direct_joint(after) - direct_joint(before)
      observed <- state[["prior_evaluator"]](after) - state[["prior_evaluator"]](before)
      expect_lt(abs(observed - joint_ratios[[i]]), 1e-12)
      marginal_ratios[[i]] <- stats::dt(1 / .5, df = 3, log = TRUE) -
        stats::dt(before[, parameter] / .5, df = 3, log = TRUE)
    }
    expect_gt(max(abs(joint_ratios - marginal_ratios)), 1e-5)
  }

  # This default native-positive scale batch carries the unchanged full vector
  # prior and exercises a finite value outside the nonnegative SD support.
  tau_spec <- .iwmde_parameter_spec(context, "tau", NULL)
  tau_states <- .iwmde_row_states(context, 1:2, "tau", tau_spec)
  tau_values <- c(-1, 0, .05, .2, .6)
  tau_replacement <- .iwmde_replacement_spec(context, "tau", tau_spec)
  tau_routes <- .normal_grid_both_routes(context, "tau", list(values = tau_values,
    row_states = tau_states, replacement = tau_replacement))
  expect_gt(tau_routes[["served"]], 0L)
  .expect_normal_grid_agreement(tau_routes, "dependent Student-t nuisance / tau")
  valid_tau <- tau_values >= 0
  tau_direct <- vapply(tau_states, function(state) vapply(tau_values[valid_tau],
    function(value) .iwmde_log_q_replacement(context, state, "tau", value, tau_replacement),
    numeric(1)), numeric(sum(valid_tau)))
  native_valid_tau <- tau_routes[["native"]][valid_tau, , drop = FALSE]
  expect_identical(native_valid_tau == -Inf, tau_direct == -Inf)
  expect_true(all(tau_routes[["native"]][1L, ] == -Inf))
  finite <- is.finite(tau_direct)
  expect_lt(max(abs(native_valid_tau[finite] - tau_direct[finite])), 1e-10)
  expect_error(.iwmde_log_q_replacement(context, tau_states[[1L]], "tau", -1, tau_replacement),
    "must be a finite, non-negative numeric matrix.", fixed = TRUE)

  spec <- .iwmde_parameter_spec(context, "mu_intercept", NULL)
  for (states in list(.iwmde_row_states(context, 1:2, "mu_intercept", spec),
                     .iwmde_row_states_grouped_marginal(context, 1:2, "mu_intercept", spec))) {
    expect_true(all(vapply(states, `[[`, logical(1), "use_focal_prior_delta")))
    expect_true(all(is.finite(vapply(states, `[[`, numeric(1), "baseline_focal_log_prior"))))
  }
})


test_that("ordered nuisance priors use the same numerical kernel on scalar and grid rows", {

  fit_name <- "lane_carrier_gated"
  skip_if_missing_fits(fit_name)
  object <- load_fit(fit_name)
  context <- .iwmde_context(object)
  indicator <- "mu_g1_ordered_total_indicator"
  keys <- .iwmde_active_keys(context)
  expect_identical(.iwmde_active_keys_matrix(context, context[["posterior_samples"]]), keys)
  rows <- vapply(c(0, 1), function(value) {
    which(context[["posterior_samples"]][, indicator] == value)[[1L]]
  }, integer(1))
  for (row in rows) {
    expect_identical(.iwmde_active_key(context, context[["posterior_samples"]][row, ]), keys[[row]])
    expect_identical(.iwmde_active_keys_matrix(context, context[["posterior_samples"]][row, , drop = FALSE]), keys[row])
  }
  for (value in c(-1, .5, NA_real_)) {
    row <- context[["posterior_samples"]][rows[[1L]], ]
    row[[indicator]] <- value
    if (is.finite(value) && value == trunc(value)) {
      expect_identical(.iwmde_active_key(context, row), paste0(indicator, "=-1"))
      expect_identical(.iwmde_active_keys_matrix(context, matrix(row, nrow = 1L,
        dimnames = list(NULL, names(row)))), paste0(indicator, "=-1"))
    } else {
      expect_error(.iwmde_active_key(context, row))
      expect_error(.iwmde_active_keys_matrix(context, matrix(row, nrow = 1L,
        dimnames = list(NULL, names(row)))))
    }
    prior <- attr(object[["fit"]], "prior_list")[["mu_g1"]]
    kernel <- BayesTools::JAGS_ordered_density_kernel(list(mu_g1 = prior))
    expect_error(kernel(matrix(row, nrow = 1L, dimnames = list(NULL, names(row)))),
      class = "BayesTools_ordered_invalid_state")
  }
  mixture <- BayesTools::prior_mixture(list(BayesTools::prior("normal", list(0, 1)), BayesTools::prior("normal", list(1, 1))))
  for (value in c(-1, 0, 3)) {
    expect_error(.iwmde_select_prior_component(mixture, c(indicator = value), "indicator"),
      "outside the available prior components", fixed = TRUE)
  }
  parameter <- "tau"
  spec <- .iwmde_parameter_spec(context, parameter, NULL)
  states <- .iwmde_row_states(context, rows, parameter, spec)
  grouped <- .iwmde_row_states_grouped_marginal(context, rows, parameter, spec)
  expect_equal(vapply(states, `[[`, numeric(1), "baseline_log_q"),
    vapply(grouped, `[[`, numeric(1), "baseline_log_q"), tolerance = 1e-12)
  inputs <- list(values = c(.05, .2, .6, -1, 0), row_states = states,
    replacement = .iwmde_replacement_spec(context, parameter, spec))
  routes <- .normal_grid_both_routes(context, parameter, inputs)
  expect_gt(routes[["served"]], 0L)
  .expect_normal_grid_agreement(routes, "ordered nuisance / tau")
  for (state in states) {
    samples <- matrix(state[["row"]], nrow = 1L, dimnames = list(NULL, names(state[["row"]])))
    observed <- .iwmde_log_prior_rows(samples, state[["prior_list"]])
    expect_equal(observed, state[["baseline_log_prior"]], tolerance = 1e-12)
    expect_equal(.iwmde_log_prior_row(state[["row"]], state[["prior_list"]]), observed[[1L]], tolerance = 1e-12)
    expect_true(is.finite(state[["baseline_log_q"]]))
    prior <- state[["prior_list"]][["mu_g1"]]
    ordered <- BayesTools::JAGS_ordered_density_kernel(list(mu_g1 = prior))(samples)
    total <- state[["row"]][["mu_g1_ordered_total"]]
    expected_ordered <- if (state[["row"]][[indicator]] == 0) 0 else stats::dgamma(total, 2, rate = 2, log = TRUE)
    expect_lt(abs(ordered - expected_ordered), 1e-12)
    expect_lt(abs(observed - ordered - BayesTools::lpdf(state[["prior_list"]][["tau"]], samples[, "tau"])), 1e-12)
    mean <- 1 + c(0, .25, 1)[as.integer(context[["data"]][["mods"]][["g1"]])] * total
    yi <- .outcome_data_yi(object)
    sei <- .outcome_data_sei(object)
    for (value in inputs[["values"]][inputs[["values"]] > 0]) {
      expected <- sum(stats::dnorm(yi, mean, sqrt(sei^2 + value^2), log = TRUE)) +
        expected_ordered + BayesTools::lpdf(state[["prior_list"]][["tau"]], value)
      actual <- .iwmde_log_q_replacement(context, state, parameter, value, inputs[["replacement"]])
      expect_lt(abs(actual - expected), 1e-10)
    }
  }
  g_spec <- .iwmde_parameter_spec(context, "mu_g1[20]", list(type = "linear",
    weights = c(`mu_g1[1]` = 1, `mu_g1[2]` = 1)))
  expect_identical(g_spec[["type"]], "ordered")
  ordered_state <- .iwmde_ordered_row_state(context, g_spec, rows[[2L]])
  expect_true(is.finite(ordered_state[["baseline_log_q"]]))
  source <- ordered_state[["ordered_chart"]][["current"]]
  expect_true(all(is.finite(.iwmde_ordered_log_q_source(context, ordered_state,
    c(source / 2, source, source * 2)))))
  components <- .iwmde_parameter_components(context, "mu_g1[20]", g_spec)
  expect_identical(sum(components[["active"]]), 5L)
  curve <- graphics::plot(object, parameter = "tau", prior = FALSE,
    density_method = "qCMDE", plot_type = "ggplot")
  expect_s3_class(curve, "ggplot")
  curve_data <- ggplot2::ggplot_build(curve)$data[[1L]]
  expect_true(any(is.finite(curve_data$y) & curve_data$y > 0))
  expect_error(graphics::plot(object, parameter = "g1[20]", prior = FALSE,
    density_method = "qCMDE", plot_type = "ggplot"), class = "RoBMA_density_plot_error")
  expect_error(graphics::plot(object, parameter = "g1[20]", prior = FALSE,
    density_method = "qCMDE", plot_type = "ggplot", conditional = TRUE),
    "'conditional' plots are available only for RoBMA objects.", fixed = TRUE)
})


test_that("expression nuisance and selected cells retain their declared typed refusals", {

  fit_name <- "marginal_means_expression"
  skip_if_missing_fits(fit_name)
  object <- load_fit(fit_name)
  context <- .iwmde_context(object)
  for (parameter in c("mu_intercept", "tau")) {
    spec <- .iwmde_parameter_spec(context, parameter, NULL)
    expect_error(.iwmde_row_states(context, 1:2, parameter, spec),
      class = "BayesTools_ordered_expression_unavailable")
    expect_null(.normal_grid_inputs(context, parameter))
  }
  for (parameter in c("intercept", "tau")) {
    expect_error(graphics::plot(object, parameter = parameter, prior = FALSE,
      density_method = "qCMDE", plot_type = "ggplot"),
      class = "BayesTools_ordered_expression_unavailable")
  }
  expect_error(graphics::plot(object, parameter = "g[hi]", prior = FALSE,
    density_method = "qCMDE", plot_type = "ggplot"),
    class = "BayesTools_formula_prior_density_unavailable")
  whole <- tryCatch(graphics::plot(object, parameter = "g", prior = FALSE,
    density_method = "qCMDE", plot_type = "ggplot"), error = identity)
  expect_true(inherits(whole, "BayesTools_ordered_expression_unavailable") ||
    inherits(whole, "BayesTools_formula_measure_unavailable"))
})


test_that("the native normal candidate grid matches the candidate route", {

  fit_names <- .normal_grid_fit_names()
  skip_if(length(fit_names) == 0L, "No cached normal fixtures are available.")

  served   <- character()
  declined <- character()

  for (fit_name in fit_names) {
    object <- tryCatch(load_fit(fit_name, validate = FALSE), error = function(e) NULL)
    if (is.null(object)) {
      next
    }
    context <- .iwmde_context(object)

    for (parameter in .normal_grid_parameters(context)) {
      inputs <- .normal_grid_inputs(context, parameter)
      if (is.null(inputs)) {
        next
      }
      label <- paste0(fit_name, " / ", parameter)

      routes <- .normal_grid_both_routes(context, parameter, inputs)
      if (routes[["served"]] == 0L || !is.matrix(routes[["native"]])) {
        declined <- c(declined, label)
        next
      }
      .expect_normal_grid_agreement(routes, label)
      served <- c(served, label)

      # A value outside the parameter's support exercises the validity mask of
      # both routes. Only a batch the native route serves for every state group
      # can take it: the candidate route evaluates a scale formula on every
      # candidate row, valid or not, and a formula batch is not this route's.
      if (routes[["served"]] != routes[["groups"]]) {
        next
      }
      outside <- inputs
      outside[["values"]] <- c(inputs[["values"]], -1)
      masked <- .normal_grid_both_routes(context, parameter, outside)
      if (masked[["served"]] != masked[["groups"]]) {
        next
      }
      .expect_normal_grid_agreement(masked, paste0(label, " (outside support)"))
    }
  }

  skip_if(length(served) == 0L,
          "No cached normal fixture reaches the native candidate grid.")
  expect_gt(length(served), 0L)
  # Recorded so a guard change is visible in the test output.
  cat("\nnative normal candidate grid served:\n  ",
      paste(served, collapse = "\n  "), "\n", sep = "")
  if (length(declined)) {
    cat("native grid not observed (earlier or candidate route retained):\n  ",
        paste(declined, collapse = "\n  "), "\n", sep = "")
  }
})


# The `tau` of a scale-regression model is the logged formula intercept. It is
# affine in its own logarithm, so the batched route carries a log-tau basis and
# forms the update from log(value) - log(current); the generic route rebuilds
# the whole scale formula for every candidate row. The two must agree to
# floating-point rounding.
test_that("the log-tau intercept basis matches the generic formula evaluator", {

  fit_names <- .normal_grid_fit_names()
  skip_if(length(fit_names) == 0L, "No cached normal fixtures are available.")

  parameter <- "log_tau_intercept"
  served    <- character()

  for (fit_name in fit_names) {
    object <- tryCatch(load_fit(fit_name, validate = FALSE), error = function(e) NULL)
    if (is.null(object) || !.is_scale(object)) {
      next
    }
    context <- .iwmde_context(object)
    if (!parameter %in% colnames(context[["posterior_samples"]])) {
      next
    }
    inputs <- .normal_grid_inputs(context, parameter)
    if (is.null(inputs)) {
      next
    }
    label <- paste0(fit_name, " / ", parameter)

    log_bases <- 0L
    batched <- local({
      original <- .iwmde_predictor_materialize_formula_basis
      testthat::local_mocked_bindings(
        .iwmde_predictor_materialize_formula_basis = function(...) {
          basis <- original(...)
          if (identical(basis[["log_tau_basis_coordinate"]], "log")) {
            log_bases <<- log_bases + 1L
          }
          basis
        },
        .package = "RoBMA"
      )
      .iwmde_log_q_grid_predictor_batch(
        context     = context,
        parameter   = parameter,
        values      = inputs[["values"]],
        row_states  = inputs[["row_states"]],
        replacement = inputs[["replacement"]]
      )
    })
    # Restoring the previous non-affine verdict keeps the generic evaluator.
    generic <- local({
      original <- .iwmde_predictor_materialize_formula_basis
      testthat::local_mocked_bindings(
        .iwmde_predictor_materialize_formula_basis = function(...) {
          basis <- original(...)
          if (identical(basis[["log_tau_basis_coordinate"]], "log")) {
            basis[["log_tau_basis"]]            <- NULL
            basis[["log_tau_basis_coordinate"]] <- NULL
            basis[["formula_logtau"]]           <- TRUE
            basis[["formula_logtau_columns"]]   <- parameter
          }
          basis
        },
        .package = "RoBMA"
      )
      .iwmde_log_q_grid_predictor_batch(
        context     = context,
        parameter   = parameter,
        values      = inputs[["values"]],
        row_states  = inputs[["row_states"]],
        replacement = inputs[["replacement"]]
      )
    })
    if (log_bases == 0L || !is.matrix(batched) || !is.matrix(generic)) {
      next
    }

    expect_identical(dim(batched), dim(generic), info = label)
    expect_identical(is.finite(batched), is.finite(generic), info = label)
    finite <- is.finite(generic)
    expect_true(any(finite), info = paste0(label, ": no finite log density"))
    deviation <- max(abs(batched[finite] - generic[finite]) /
                       pmax(abs(generic[finite]), 1))
    expect_lt(deviation, 1e-10, label = paste0(label, ": max scaled deviation"))

    # A non-positive candidate has no logarithm on either route.
    outside <- inputs
    outside[["values"]] <- c(inputs[["values"]], -1, 0)
    masked <- .iwmde_log_q_grid_predictor_batch(
      context     = context,
      parameter   = parameter,
      values      = outside[["values"]],
      row_states  = outside[["row_states"]],
      replacement = outside[["replacement"]]
    )
    expect_true(is.matrix(masked), info = paste0(label, " (outside support)"))
    expect_true(all(masked[length(inputs[["values"]]) + 1:2, ] == -Inf),
                info = paste0(label, " (outside support)"))

    served <- c(served, label)
  }

  skip_if(length(served) == 0L,
          "No cached scale-regression fixture carries a log-tau intercept basis.")
  cat("\nlog-tau intercept basis served:\n  ",
      paste(served, collapse = "\n  "), "\n", sep = "")
})


# The batched route is answered by one estimate of one batch above. What a user
# reads is the whole `tau` density line of a scale model and its point
# ordinate, which the estimator assembles from many batches, its own
# normalization grid and the row weights. Drive both routes through the
# estimator and compare what it returns.
.log_tau_route_estimates <- function(context, parameter, value) {

  control <- list(n_points = 20L, samples = 200L)
  run <- function() {
    line <- .iwmde_estimate(
      context         = context,
      parameter       = parameter,
      density_method  = "qCMDE",
      density_control = control,
      outputs         = "density",
      parameter_spec  = list(
        type             = "primitive",
        conditional      = NULL,
        conditional_rule = "AND"
      ),
      metadata        = list(parameter = parameter),
      cache           = .iwmde_estimate_cache()
    )
    point <- .iwmde_estimate(
      context         = context,
      parameter       = parameter,
      density_method  = "qCMDE",
      density_control = c(control, list(display_grid = "ordinate")),
      outputs         = "ordinate",
      values          = value,
      parameter_spec  = list(
        type             = "primitive",
        conditional      = NULL,
        conditional_rule = "AND"
      ),
      metadata        = list(parameter = parameter),
      cache           = .iwmde_estimate_cache()
    )
    list(
      x        = as.numeric(line[["posterior_density"]][["x"]]),
      y        = as.numeric(line[["posterior_density"]][["y"]]),
      ordinate = as.numeric(point[["posterior_ordinate"]][["ordinate"]]),
      density  = line[["diagnostics"]][["density"]]
    )
  }

  log_bases <- 0L
  batched <- local({
    original <- .iwmde_predictor_materialize_formula_basis
    testthat::local_mocked_bindings(
      .iwmde_predictor_materialize_formula_basis = function(...) {
        basis <- original(...)
        if (identical(basis[["log_tau_basis_coordinate"]], "log")) {
          log_bases <<- log_bases + 1L
        }
        basis
      },
      .package = "RoBMA"
    )
    run()
  })
  # Restoring the previous non-affine verdict sends every candidate row through
  # .iwmde_predictor_evaluate_tau(), which rebuilds the whole scale formula.
  generic <- local({
    original <- .iwmde_predictor_materialize_formula_basis
    testthat::local_mocked_bindings(
      .iwmde_predictor_materialize_formula_basis = function(...) {
        basis <- original(...)
        if (identical(basis[["log_tau_basis_coordinate"]], "log")) {
          basis[["log_tau_basis"]]            <- NULL
          basis[["log_tau_basis_coordinate"]] <- NULL
          basis[["formula_logtau"]]           <- TRUE
          basis[["formula_logtau_columns"]]   <- parameter
        }
        basis
      },
      .package = "RoBMA"
    )
    run()
  })

  return(list(batched = batched, generic = generic, log_bases = log_bases))
}


test_that("the scale model's `tau` density line matches the generic evaluator", {

  skip_if_not_certification(
    "The generic-evaluator density-line comparison runs in the normal-models case."
  )
  fit_names <- .normal_grid_fit_names()
  skip_if(length(fit_names) == 0L, "No cached normal fixtures are available.")

  parameter <- "log_tau_intercept"
  served    <- character()
  rejected  <- character()

  for (fit_name in fit_names) {
    object <- tryCatch(load_fit(fit_name, validate = FALSE), error = function(e) NULL)
    if (is.null(object) || !.is_scale(object)) {
      next
    }
    context <- .iwmde_context(object)
    if (!parameter %in% colnames(context[["posterior_samples"]])) {
      next
    }
    draws <- context[["posterior_samples"]][, parameter]
    draws <- draws[is.finite(draws) & draws > 0]
    if (length(draws) < 2L) {
      next
    }
    label     <- paste0(fit_name, " / ", parameter)
    estimates <- .log_tau_route_estimates(
      context   = context,
      parameter = parameter,
      value     = as.numeric(stats::median(draws))
    )
    if (estimates[["log_bases"]] == 0L) {
      next
    }
    batched <- estimates[["batched"]]
    generic <- estimates[["generic"]]

    # A model-averaged scale model spreads the fixed 200 estimator rows over
    # several active states, below the density minimum for such lines; the
    # estimator then rejects the line on both routes and returns none, so there
    # is no line to compare. Anything else that empties a line fails here.
    if (length(batched[["y"]]) == 0L || length(generic[["y"]]) == 0L) {
      for (route in list(batched, generic)) {
        diagnostics <- route[["density"]][["diagnostics"]]
        expect_length(route[["y"]], 0L)
        expect_identical(route[["density"]][["status"]], "ok",
                         info = paste0(label, ": rejected, not unavailable"))
        expect_false(is.null(.iwmde_diagnostics_density_failure_reason(diagnostics)),
                     info = paste0(label, ": rejected by diagnostics"))
        expect_gt(diagnostics[["n_active_state_keys"]], 1)
        expect_lt(diagnostics[["n_estimator_rows"]],
                  .iwmde_density_min_estimator_rows())
      }
      rejected <- c(rejected, label)
      next
    }

    expect_identical(batched[["x"]], generic[["x"]],
                     info = paste0(label, ": display grid"))
    expect_gt(length(generic[["y"]]), 1L)
    expect_length(batched[["y"]], length(generic[["y"]]))
    expect_true(all(is.finite(generic[["y"]])),
                info = paste0(label, ": generic density line"))
    line_deviation <- max(abs(batched[["y"]] - generic[["y"]]) /
                            pmax(abs(generic[["y"]]), .Machine$double.xmin))
    expect_lt(line_deviation, 1e-10,
              label = paste0(label, ": max relative density-line difference"))

    expect_true(is.finite(generic[["ordinate"]]),
                info = paste0(label, ": generic ordinate"))
    expect_gt(generic[["ordinate"]], 0)
    point_deviation <- abs(batched[["ordinate"]] - generic[["ordinate"]]) /
      abs(generic[["ordinate"]])
    expect_lt(point_deviation, 1e-10,
              label = paste0(label, ": relative ordinate difference"))

    served <- c(served, sprintf(
      "%s (line %.3g, ordinate %.3g)", label, line_deviation, point_deviation
    ))
  }

  skip_if(length(served) == 0L,
          "No cached scale-regression fixture carries a log-tau intercept basis.")
  cat("\n`tau` density line, batched vs generic:\n  ",
      paste(served, collapse = "\n  "), "\n", sep = "")
  if (length(rejected)) {
    cat("rejected by the density diagnostics on both routes:\n  ",
        paste(rejected, collapse = "\n  "), "\n", sep = "")
  }
})


test_that("the native normal candidate grid is thread invariant", {

  .with_kernel_thread_worker({

    expect_true(is.loaded("RoBMA_norm_predictor_grid_loglik", PACKAGE = "RoBMA"))
    previous_threads <- RoBMA.get_option("native_threads")
    withr::defer(RoBMA.options(native_threads = previous_threads))

    withr::local_seed(20260918)
    S   <- 1000L
    K   <- 9L
    G   <- 7L
    yi  <- stats::rnorm(K, sd = .3)
    sei <- stats::runif(K, .05, .40)
    mu  <- matrix(stats::rnorm(S * K, sd = .2), nrow = S)
    mu_basis      <- matrix(stats::rnorm(S * K), nrow = S)
    tau           <- matrix(stats::runif(S * K, 0, .5), nrow = S)
    log_tau_basis <- matrix(stats::rnorm(S * K, sd = .3), nrow = S)
    current       <- stats::rnorm(S, sd = .2)
    values        <- c(-0.4, 0, 0.2, 0.7, 1.5, NaN, Inf)
    weights       <- stats::runif(K, .5, 2)

    call_grid <- function(scale_tau, use_mu_basis, use_log_tau, use_weights,
                          log_delta = FALSE) {
      .Call("RoBMA_norm_predictor_grid_loglik",
        yi, sei, if (use_weights) weights else NULL, mu,
        if (use_mu_basis) mu_basis else NULL,
        if (scale_tau) NULL else tau,
        if (use_log_tau) log_tau_basis else NULL,
        current, values, scale_tau, log_delta, PACKAGE = "RoBMA")
    }

    for (scale_tau in c(FALSE, TRUE)) {
      for (use_log_tau in c(FALSE, TRUE)) {
        reference <- NULL
        for (threads in c(1L, 2L, 8L)) {
          RoBMA.options(native_threads = threads)
          if (threads > 1L) {
            schedule <- .Call("RoBMA_selnorm_row_schedule", as.double(S),
              as.double(G * K), PACKAGE = "RoBMA")
            expect_gt(schedule[["threads"]], 1L)
          }
          value <- call_grid(scale_tau, TRUE, use_log_tau, TRUE)
          if (is.null(reference)) {
            expect_true(any(is.finite(value[["log_lik"]])))
            expect_true(any(value[["valid"]]))
            reference <- value
          } else {
            expect_identical(value, reference,
              info = paste("scale_tau", scale_tau, "log_tau", use_log_tau,
                           "threads", threads))
          }
        }
      }
    }

    # The log-coordinate delta of a logged formula intercept is equally thread
    # invariant, and it accompanies no location basis.
    positive_current <- abs(current) + .1
    log_reference <- NULL
    for (threads in c(1L, 2L, 8L)) {
      RoBMA.options(native_threads = threads)
      if (threads > 1L) {
        schedule <- .Call("RoBMA_selnorm_row_schedule", as.double(S),
          as.double(G * K), PACKAGE = "RoBMA")
        expect_gt(schedule[["threads"]], 1L)
      }
      value <- .Call("RoBMA_norm_predictor_grid_loglik",
        yi, sei, weights, mu, NULL, tau, log_tau_basis,
        positive_current, values, FALSE, TRUE, PACKAGE = "RoBMA")
      if (is.null(log_reference)) {
        expect_true(any(is.finite(value[["log_lik"]])))
        expect_true(any(value[["valid"]]))
        log_reference <- value
      } else {
        expect_identical(value, log_reference, info = paste("log delta, threads", threads))
      }
    }
    # Only a positive candidate and a positive current value have a logarithm.
    non_positive <- !is.na(values) & values <= 0
    expect_true(all(!log_reference[["valid"]][rep(non_positive, times = S)]))
    expect_true(any(log_reference[["valid"]]))
    expect_error(
      .Call("RoBMA_norm_predictor_grid_loglik", yi, sei, weights, mu, mu_basis,
        tau, log_tau_basis, positive_current, values, FALSE, TRUE, PACKAGE = "RoBMA"),
      "location basis"
    )
    expect_error(
      .Call("RoBMA_norm_predictor_grid_loglik", yi, sei, weights, mu, NULL,
        tau, NULL, positive_current, values, FALSE, TRUE, PACKAGE = "RoBMA"),
      "log-tau basis"
    )

    # A grid without a mu basis and without weights keeps the same shapes.
    plain <- .expect_thread_invariant("plain normal predictor", function() {
      call_grid(FALSE, FALSE, FALSE, FALSE)
    }, rows = S, work_per_row = G * K)
    expect_identical(length(plain[["log_lik"]]), length(values) * S)
    expect_identical(length(plain[["valid"]]), length(values) * S)
    expect_true(is.logical(plain[["valid"]]))
  })
})
