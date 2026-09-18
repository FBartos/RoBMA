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

  spec <- tryCatch(
    .iwmde_parameter_spec(context, parameter, NULL),
    error = function(e) NULL
  )
  if (is.null(spec) || !identical(spec[["status"]], "ok")) {
    return(NULL)
  }
  values <- tryCatch(
    .iwmde_parameter_values(context, parameter, spec),
    error = function(e) NULL
  )
  component <- tryCatch(
    .iwmde_parameter_components(context, parameter, spec),
    error = function(e) NULL
  )
  if (is.null(values) || is.null(component)) {
    return(NULL)
  }
  active <- component[["active"]] & is.finite(values)
  if (sum(active) < 2L) {
    return(NULL)
  }
  row_states <- tryCatch(
    .iwmde_row_states(context, utils::head(which(active), n_rows), parameter, spec),
    error = function(e) NULL
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


test_that("the native normal candidate grid is thread invariant", {

  skip_if_not(is.loaded("RoBMA_norm_predictor_grid_loglik", PACKAGE = "RoBMA"))

  set.seed(20260918)
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

  call_grid <- function(scale_tau, use_mu_basis, use_log_tau, use_weights) {
    .Call("RoBMA_norm_predictor_grid_loglik",
      yi, sei, if (use_weights) weights else NULL, mu,
      if (use_mu_basis) mu_basis else NULL,
      if (scale_tau) NULL else tau,
      if (use_log_tau) log_tau_basis else NULL,
      current, values, scale_tau, PACKAGE = "RoBMA")
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
  RoBMA.options(native_threads = 1L)

  # A grid without a mu basis and without weights keeps the same shapes.
  plain <- call_grid(FALSE, FALSE, FALSE, FALSE)
  expect_identical(length(plain[["log_lik"]]), length(values) * S)
  expect_identical(length(plain[["valid"]]), length(values) * S)
  expect_true(is.logical(plain[["valid"]]))
})
