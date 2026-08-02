source(testthat::test_path("common-functions.R"))


test_that("known-V tau q-grid factors covariance once per value and block", {

  fit_name <- "brma.mv_block_mvn"
  skip_if_missing_fits(fit_name)

  context   <- .iwmde_context(load_fit(fit_name))
  parameter <- "tau"
  spec      <- .iwmde_parameter_spec(context, parameter, NULL)
  samples   <- .iwmde_parameter_values(context, parameter, spec)
  component <- .iwmde_parameter_components(context, parameter, spec)
  rows      <- which(component[["active"]] & is.finite(samples))
  row_states <- .iwmde_row_states(context, rows, parameter, spec)
  keep       <- vapply(row_states, function(state) {
    is.finite(state[["baseline_log_q"]])
  }, logical(1))
  row_states <- row_states[keep]

  expect_gt(length(row_states), 1L)

  grid <- c(
    0,
    sqrt(.Machine$double.eps),
    as.numeric(stats::quantile(
      samples[component[["active"]] & is.finite(samples)],
      probs = c(.25, .50, .75),
      names = FALSE,
      type  = 8
    ))
  )
  grid        <- unique(grid)
  replacement <- .iwmde_replacement_spec(context, parameter, spec)

  original_factor <- .known_v_chol_covariance
  factor_calls     <- 0L
  testthat::local_mocked_bindings(
    .known_v_chol_covariance = function(covariance, context) {

      factor_calls <<- factor_calls + 1L
      original_factor(covariance = covariance, context = context)
    },
    .package = "RoBMA"
  )

  fast <- .iwmde_log_q_grid_predictor_batch(
    context     = context,
    parameter   = parameter,
    values      = grid,
    row_states  = row_states,
    replacement = replacement
  )
  fast_factor_calls <- factor_calls
  factor_calls      <- 0L
  scalar <- .iwmde_log_q_grid_scalar(
    context     = context,
    parameter   = parameter,
    values      = grid,
    row_states  = row_states,
    replacement = replacement
  )

  block_count <- length(.known_v_blocks(
    .data_known_v_data(context[["data"]])
  ))
  expect_true(is.matrix(fast))
  expect_equal(dim(fast), dim(scalar))
  expect_equal(is.finite(fast), is.finite(scalar))
  expect_equal(fast, scalar, tolerance = 1e-8)
  expect_equal(fast_factor_calls, length(grid) * block_count)
  expect_equal(factor_calls, length(grid) * length(row_states) * block_count)
})


test_that("ordinary tau q-grid preserves the exact zero boundary", {

  fit_name <- "bcg_meta-analysis"
  skip_if_missing_fits(fit_name)

  context   <- .iwmde_context(load_fit(fit_name))
  parameter <- "tau"
  spec      <- .iwmde_parameter_spec(context, parameter, NULL)
  samples   <- .iwmde_parameter_values(context, parameter, spec)
  component <- .iwmde_parameter_components(context, parameter, spec)
  rows      <- head(
    which(component[["active"]] & is.finite(samples)),
    3L
  )
  row_states <- .iwmde_row_states(context, rows, parameter, spec)
  keep       <- vapply(row_states, function(state) {
    is.finite(state[["baseline_log_q"]])
  }, logical(1))
  row_states <- row_states[keep]

  expect_gt(length(row_states), 0L)

  grid <- c(
    0,
    sqrt(.Machine$double.eps),
    stats::median(samples[component[["active"]] & is.finite(samples)])
  )
  replacement <- .iwmde_replacement_spec(context, parameter, spec)
  fast <- .iwmde_log_q_grid_predictor_batch(
    context     = context,
    parameter   = parameter,
    values      = grid,
    row_states  = row_states,
    replacement = replacement
  )
  scalar <- .iwmde_log_q_grid_scalar(
    context     = context,
    parameter   = parameter,
    values      = grid,
    row_states  = row_states,
    replacement = replacement
  )

  expect_true(is.matrix(fast))
  expect_equal(dim(fast), dim(scalar))
  expect_true(all(is.finite(fast[grid == 0, ])))
  expect_equal(fast, scalar, tolerance = 1e-8)
})


test_that("IWMDE preflight rejects an unregularized singular known-V tau null", {

  K       <- 3L
  known_V <- .known_v_prepare(
    V                         = matrix(1, nrow = K, ncol = K),
    keep_rows                 = rep(TRUE, K),
    known_v_parameterization  = "block_mvn",
    known_v_residual_fraction = NULL,
    warn_singular             = FALSE
  )
  data    <- structure(list(), class = "RoBMA_data")
  attr(data, "known_V")      <- TRUE
  attr(data, "known_V_data") <- known_V
  attr(data, "random")       <- FALSE

  context <- list(
    data   = data,
    object = list(data = data)
  )
  plan    <- list(
    target         = list(parameter = "tau"),
    parameter_spec = list(type = "primitive"),
    outputs        = list(requested_values = 0)
  )
  actual  <- .iwmde_plan_prepare_contract(context, plan)

  expect_equal(actual[["status"]], "unsupported")
  expect_match(actual[["reason"]], "tau = 0 likelihood is degenerate")
  expect_match(actual[["reason"]], "singular known-V dependency block")
})
