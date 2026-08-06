context("IWMDE conditional-weight fallbacks")

source(testthat::test_path("common-functions.R"))


test_that("IWMDE parameter setup propagates backend programming errors", {

  testthat::local_mocked_bindings(
    JAGS_marglik_parameters = function(...) {

      stop("sentinel marginal-parameter error", call. = FALSE)
    },
    .package = "BayesTools"
  )

  expect_error(
    .iwmde_row_parameters(
      context      = list(data = list(outcome = data.frame(yi = 0))),
      row          = list(mu = 0),
      active_setup = list(
        fit_priors = list(mu = BayesTools::prior(
          "normal",
          parameters = list(mean = 0, sd = 1)
        )),
        priors = list()
      )
    ),
    "sentinel marginal-parameter error",
    fixed = TRUE
  )
})


test_that("IWMDE prior reconstruction propagates programming errors", {

  data <- list(outcome = data.frame(yi = 0))
  attr(data, "outcome_type") <- "norm"

  testthat::local_mocked_bindings(
    .create_fit_priors = function(...) {

      stop("sentinel prior-construction error", call. = FALSE)
    },
    .package = "RoBMA"
  )

  expect_error(
    .iwmde_resolve_fixed_prior_sample_columns(
      context      = list(data = data),
      samples      = matrix(0, nrow = 1L, ncol = 1L),
      active_setup = list(fit_priors = NULL, priors = list(mu = list()))
    ),
    "sentinel prior-construction error",
    fixed = TRUE
  )
})


test_that("conditional-weight fallback preserves method and failure reason", {

  success <- .iwmde_chen_try_weight(
    expr          = list(log_weight = 0, method = "conditional"),
    fallback      = list(log_weight = -1, method = "marginal"),
    fallback_from = "conditional"
  )

  expect_equal(success[["method"]], "conditional")
  expect_true(is.na(success[["fallback_from"]]))
  expect_true(is.na(success[["fallback_reason"]]))

  fallback <- .iwmde_chen_try_weight(
    expr          = .iwmde_chen_conditional_stop(
      "singular conditional covariance"
    ),
    fallback      = list(log_weight = -1, method = "marginal"),
    fallback_from = "conditional"
  )

  expect_equal(fallback[["method"]], "marginal")
  expect_equal(fallback[["log_weight"]], -1)
  expect_equal(fallback[["fallback_from"]], "conditional")
  expect_equal(
    fallback[["fallback_reason"]],
    paste0(
      "IWMDE Chen conditional-normal weights are unavailable: ",
      "singular conditional covariance."
    )
  )

  expect_error(
    .iwmde_chen_try_weight(
      expr          = stop("programming defect", call. = FALSE),
      fallback      = list(log_weight = -1, method = "marginal"),
      fallback_from = "conditional"
    ),
    "programming defect"
  )
})


test_that("IWMDE parameter and prior evaluation errors propagate", {

  testthat::local_mocked_bindings(
    is.prior.none = function(...) FALSE,
    lpdf = function(...) {
      stop("prior density defect", call. = FALSE)
    },
    .package = "BayesTools"
  )
  expect_error(
    .iwmde_focal_log_prior_values(
      prior  = structure(list(), class = "mock_prior"),
      values = 0
    ),
    "prior density defect"
  )
})


test_that("Chen dispatcher aggregates fallback partitions and rows", {

  supports <- matrix(
    c(-Inf, Inf,
      -Inf, Inf,
      0, 1,
      0, 1),
    ncol  = 2L,
    byrow = TRUE
  )
  rownames(supports) <- as.character(seq_len(4L))

  testthat::local_mocked_bindings(
    .iwmde_parameter_row_supports = function(context, parameter, rows,
                                             parameter_spec) {

      context[["supports"]][as.character(rows), , drop = FALSE]
    },
    .iwmde_chen_conditional_normal_log_weight = function(...) {

      .iwmde_chen_conditional_stop("singular covariance")
    },
    .iwmde_chen_marginal_normal_log_weight = function(active_values,
                                                      weight_values) {

      list(
        log_weight = rep(log(.25), length(active_values)),
        method     = "chen_marginal_normal"
      )
    },
    .iwmde_chen_logit_conditional_normal_log_weight = function(...) {

      .iwmde_chen_conditional_stop("singular covariance")
    },
    .iwmde_chen_beta_log_weight = function(active_values, weight_values,
                                           support) {

      list(
        log_weight = rep(log(.5), length(active_values)),
        method     = "chen_beta"
      )
    },
    .package = "RoBMA"
  )

  weight <- .iwmde_chen_log_weight(
    context        = list(supports = supports),
    parameter      = "theta",
    parameter_spec = list(type = "primitive"),
    active_rows    = seq_len(4L),
    active_values  = c(-1, 1, .25, .75),
    weight_rows    = seq_len(4L),
    weight_values  = c(-1, 1, .25, .75),
    support        = c(-Inf, Inf)
  )

  expect_equal(weight[["fallback_count"]], 2L)
  expect_equal(weight[["fallback_rows"]], 4L)
  expect_equal(
    weight[["fallback_from"]],
    c("chen_conditional_normal", "chen_logit_conditional_normal")
  )
  fallback_reason <- paste0(
    "IWMDE Chen conditional-normal weights are unavailable: ",
    "singular covariance."
  )
  expect_equal(weight[["fallback_reason"]], fallback_reason)
  expect_equal(
    weight[["fallback_reasons"]],
    stats::setNames(2L, fallback_reason)
  )
  expect_equal(
    vapply(weight[["partitions"]], `[[`, character(1), "fallback_reason"),
    rep(fallback_reason, 2L)
  )
})


test_that("Chen conditioning drops exact boundary constants before transforming", {

  samples <- cbind(
    mu  = c(-.2, -.1, .1, .2),
    rho = rep(0, 4L)
  )
  context <- list(
    posterior_samples = samples,
    indicator_names   = character(),
    priors            = list(),
    flat_prior_list   = list(),
    formula_inputs    = list()
  )

  conditioning <- .iwmde_chen_conditioning_matrix(
    context        = context,
    parameter      = "mu",
    parameter_spec = list(type = "primitive", parameter = "mu"),
    active_rows    = 1:2,
    weight_rows    = 3:4
  )

  expect_identical(conditioning[["columns"]], character())
  expect_equal(dim(conditioning[["fit"]]), c(2L, 0L))
  expect_equal(dim(conditioning[["eval"]]), c(2L, 0L))
})


test_that("IWMDE density surfaces aggregate fallback diagnostics", {

  grid <- seq(0, 1, length.out = 21L)
  normalization_grid <- list(
    x            = grid,
    z            = grid,
    log_jacobian = rep(0, length(grid))
  )
  row_states <- list(
    list(baseline_log_q = 0),
    list(baseline_log_q = 0)
  )

  testthat::local_mocked_bindings(
    .iwmde_chen_log_weight = function(...) {

      list(
        log_weight       = rep(log(.5), 2L),
        method           = "chen_marginal_normal",
        fallback_from    = "chen_conditional_normal",
        fallback_reason  = "singular covariance",
        fallback_count   = 1L,
        fallback_rows    = 2L,
        fallback_reasons = c("singular covariance" = 1L),
        partitions       = list()
      )
    },
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {

      matrix(0, nrow = length(values), ncol = length(row_states))
    },
    .package = "RoBMA"
  )

  density <- .iwmde_density_iwmde(
    context            = list(),
    parameter          = "mu",
    parameter_spec     = list(type = "primitive"),
    display_grid       = grid,
    row_states         = row_states,
    active_rows        = 1:2,
    active_values      = c(.25, .75),
    weight_rows        = 1:2,
    weight_values      = c(.25, .75),
    support            = c(0, 1),
    active_mass        = 1,
    replacement        = list(type = "scalar"),
    normalization_grid = normalization_grid
  )

  expect_equal(density[["n_weight_fallbacks"]], 1L)
  expect_equal(density[["n_weight_fallback_rows"]], 2L)
  expect_equal(
    density[["weight_fallback_from"]],
    "chen_conditional_normal"
  )
  expect_equal(
    density[["weight_fallback_reasons"]],
    c("singular covariance" = 1L)
  )

  plan <- list(
    target = list(
      parameter  = "mu",
      target_key = "mu"
    ),
    rows = list(
      samples            = c(.25, .75),
      continuous_values  = c(.25, .75),
      active_mass        = 1,
      point_mass_total   = 0,
      point_masses       = numeric(),
      n_candidate_rows   = 2L,
      n_denominator_rows = 2L,
      n_estimator_rows   = 2L,
      n_total            = 2L
    ),
    support = list(
      support      = c(0, 1),
      xlim         = c(0, 1),
      density_xlim = c(0, 1)
    ),
    grids = list(
      requested_values = 0
    ),
    outputs = list(
      requested_values = 0
    ),
    control = list(
      display_grid = "ordinate"
    )
  )
  execution <- list(
    active_rows     = 1:2,
    row_states      = list(list(active_key = "all"), list(active_key = "all"))
  )
  diagnostic <- .iwmde_plan_diagnostic_result(
    plan      = plan,
    output    = "ordinate",
    execution = execution,
    density   = .iwmde_new_density_result(density, estimator = "iwmde")
  )

  expect_equal(diagnostic[["diagnostics"]][["n_weight_fallbacks"]], 1L)
  expect_equal(diagnostic[["diagnostics"]][["n_weight_fallback_rows"]], 2L)
  expect_equal(
    diagnostic[["diagnostics"]][["weight_fallback_from"]],
    "chen_conditional_normal"
  )
  expect_equal(
    diagnostic[["diagnostics"]][["weight_fallback_reasons"]],
    c("singular covariance" = 1L)
  )
  expect_true(any(grepl(
    "normalized fallback weights",
    .iwmde_diagnostics_density_warning(diagnostic[["diagnostics"]])
  )))
  expect_true(any(grepl(
    "normalized fallback weights",
    .iwmde_diagnostics_bf_warning(diagnostic[["diagnostics"]])
  )))
})
