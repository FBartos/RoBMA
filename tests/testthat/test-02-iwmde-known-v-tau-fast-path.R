source(testthat::test_path("common-functions.R"))


test_that("known-V estimate SD grids reuse the declared covariance plan", {

  fit_name <- "brma.mv_block_mvn_random"
  skip_if_missing_fits(fit_name)

  object    <- load_fit(fit_name)
  info      <- load_info(fit_name)
  context   <- .iwmde_context(object)
  target    <- .brma_random_parameter_density_target(object, "tau")
  expect_null(target[["reason"]])
  parameter <- target[["parameter"]]
  spec      <- .iwmde_parameter_spec(context, parameter, target[["parameter_spec"]])
  samples   <- .iwmde_parameter_values(context, parameter, spec)
  component <- .iwmde_parameter_components(context, parameter, spec)
  rows      <- which(component[["active"]] & is.finite(samples))
  row_states <- .iwmde_row_states(context, rows, parameter, spec)
  keep       <- vapply(row_states, function(state) {
    is.finite(state[["baseline_log_q"]])
  }, logical(1))
  row_states <- row_states[keep]

  expect_gt(length(row_states), 1L)

  grid <- unique(c(
    0,
    sqrt(.Machine$double.eps),
    as.numeric(stats::quantile(
      samples[component[["active"]] & is.finite(samples)],
      probs = c(.25, .50, .75), names = FALSE, type = 8
    ))
  ))
  replacement <- .iwmde_replacement_spec(context, parameter, spec)
  scalar <- .iwmde_log_q_grid_scalar(
    context = context, parameter = parameter, values = grid,
    row_states = row_states, replacement = replacement
  )

  original_grid <- .marglik_covariance_plan_group_iid_variance_grid_loglik
  grid_calls <- 0L
  testthat::local_mocked_bindings(
    .marglik_covariance_plan_group_iid_variance_grid_loglik = function(...) {

      grid_calls <<- grid_calls + 1L
      original_grid(...)
    },
    .package = "RoBMA"
  )
  fast <- .iwmde_log_q_grid_predictor_batch(
    context = context, parameter = parameter, values = grid,
    row_states = row_states, replacement = replacement
  )

  # The explicit unique estimate factor adds sd^2 I. Independently materialize
  # that Gaussian law and evaluate its priors using BayesTools' scalar API.
  expected <- vapply(row_states, function(state) {
    vapply(grid, function(sd) {
      row <- state[["row"]]
      row[[parameter]] <- sd
      mvtnorm::dmvnorm(
        info[["data"]][["yi"]],
        mean = rep(row[["mu_intercept"]], nrow(info[["V"]])),
        sigma = info[["V"]] + diag(sd^2, nrow(info[["V"]])),
        log = TRUE
      ) + BayesTools::JAGS_marglik_priors(row, state[["prior_list"]])
    }, numeric(1L))
  }, numeric(length(grid)))

  expect_true(is.matrix(fast))
  expect_equal(dim(fast), dim(scalar))
  expect_equal(is.finite(fast), is.finite(scalar))
  expect_true(all(is.finite(fast[grid == 0, ])))
  expect_equal(fast, scalar, tolerance = 1e-8)
  expect_equal(fast, expected, tolerance = 1e-10)
  expect_identical(grid_calls, 1L)
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
