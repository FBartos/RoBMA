# ============================================================================ #
# test-02-iwmde-normal-predictor-grid.R
# ============================================================================ #

context("IWMDE normal predictor grid")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-iwmde.R"))


# Cached fixtures that could reach the native normal candidate grid: ordinary
# normal models without a selection model. The guard inside
# .iwmde_predictor_normal_grid_log_lik() makes the final decision.
.normal_grid_fit_names <- function() {

  unique(c(
    list_fits(class = c("brma", "BMA", "RoBMA", "NoBMA")),
    list_fits(family = "norm")
  ))
}


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
  spec <- .iwmde_parameter_spec(context, parameter, NULL)
  if (is.null(spec) || !identical(spec[["status"]], "ok")) {
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
  row_states <- .iwmde_row_states(
    context, utils::head(which(active), n_rows), parameter, spec
  )
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
    cat("declined (candidate route retained):\n  ",
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
      ordinate = as.numeric(point[["posterior_ordinate"]][["ordinate"]])
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
})


test_that("the native normal candidate grid is thread invariant", {

  skip_if_not(is.loaded("RoBMA_norm_predictor_grid_loglik", PACKAGE = "RoBMA"))
  previous_threads <- RoBMA.get_option("native_threads")
  withr::defer(RoBMA.options(native_threads = previous_threads))

  withr::local_seed(20260918)
  S   <- 96L
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
        value <- call_grid(scale_tau, TRUE, use_log_tau, TRUE)
        if (is.null(reference)) {
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
    value <- .Call("RoBMA_norm_predictor_grid_loglik",
      yi, sei, weights, mu, NULL, tau, log_tau_basis,
      positive_current, values, FALSE, TRUE, PACKAGE = "RoBMA")
    if (is.null(log_reference)) {
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
  plain <- call_grid(FALSE, FALSE, FALSE, FALSE)
  expect_identical(length(plain[["log_lik"]]), length(values) * S)
  expect_identical(length(plain[["valid"]]), length(values) * S)
  expect_true(is.logical(plain[["valid"]]))
})
