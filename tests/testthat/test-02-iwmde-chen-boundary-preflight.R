# ============================================================================ #
# IWMDE Chen boundary preflight
# ============================================================================ #


test_that("Chen boundary fallback does not construct later columns", {

  samples <- cbind(
    mu         = c(-.4, -.1, .2, .5),
    PET        = c(-1, -.2, .4, 1.1),
    rho        = c(.2, .4, .6, .8),
    `omega[1]` = rep(1, 4L),
    `omega[2]` = c(.4, 1, .6, 1)
  )
  context <- list(
    posterior_samples = samples,
    indicator_names   = character(),
    selection_spec    = list(jags_omega = "omega"),
    priors            = list(),
    flat_prior_list   = list(
      mu   = TRUE,
      PET  = TRUE,
      rho  = TRUE,
      bias = BayesTools::prior_weightfunction(
        side    = "one-sided",
        steps   = .05,
        weights = BayesTools::wf_cumulative(c(1, 1))
      )
    ),
    formula_inputs    = list()
  )
  rows          <- seq_len(nrow(samples))
  active_values <- samples[, "mu"] + .05
  expected <- .iwmde_chen_marginal_normal_log_weight(
    active_values = active_values,
    weight_values = samples[, "mu"]
  )
  original_column_values <- .iwmde_parameter_column_values
  visited_columns        <- character()

  testthat::local_mocked_bindings(
    .iwmde_parameter_row_supports = function(context, parameter, rows,
                                             parameter_spec) {

      matrix(
        rep(c(-Inf, Inf), length(rows)),
        ncol  = 2L,
        byrow = TRUE
      )
    },
    .iwmde_parameter_column_values = function(context, samples, parameter) {

      visited_columns <<- c(visited_columns, parameter)
      if (identical(parameter, "PET")) {
        stop("The full nuisance matrix was constructed.", call. = FALSE)
      }

      original_column_values(
        context   = context,
        samples   = samples,
        parameter = parameter
      )
    },
    .package = "RoBMA"
  )

  weight <- .iwmde_chen_log_weight(
    context        = context,
    parameter      = "mu",
    parameter_spec = list(type = "primitive", parameter = "mu"),
    active_rows    = rows,
    active_values  = active_values,
    weight_rows    = rows,
    weight_values  = samples[, "mu"],
    support        = c(-Inf, Inf)
  )

  expect_equal(weight[["log_weight"]], expected[["log_weight"]])
  expect_identical(weight[["method"]], "chen_marginal_normal")
  expect_identical(
    weight[["fallback_reason"]],
    paste0(
      "IWMDE Chen conditional-normal weights are unavailable: ",
      "conditioning fit columns contain non-finite transformed values."
    )
  )
  expect_identical(
    visited_columns,
    c("omega[1]", "omega[1]", "omega[2]", "omega[2]")
  )
})


test_that("Chen boundary preflight retains constant-column semantics", {

  samples <- cbind(
    mu         = c(-.4, -.1, .2, .5),
    PET        = c(-1, -.2, .4, 1.1),
    rho        = rep(0, 4L),
    `omega[1]` = rep(1, 4L)
  )
  context <- list(
    posterior_samples = samples,
    indicator_names   = character(),
    selection_spec    = list(jags_omega = "omega"),
    priors            = list(),
    flat_prior_list   = list(mu = TRUE, PET = TRUE, rho = TRUE),
    formula_inputs    = list()
  )

  conditioning <- .iwmde_chen_conditioning_matrix(
    context        = context,
    parameter      = "mu",
    parameter_spec = list(type = "primitive", parameter = "mu"),
    active_rows    = 1:2,
    weight_rows    = 3:4
  )

  expect_identical(conditioning[["columns"]], "PET")
  expect_true(all(is.finite(conditioning[["fit"]])))
  expect_true(all(is.finite(conditioning[["eval"]])))
})


test_that("Chen boundary preflight does not reject evaluation-only boundaries", {

  samples <- cbind(
    mu         = c(-.4, -.1, .2, .5),
    PET        = c(-1, -.2, .4, 1.1),
    `omega[1]` = c(.2, .4, .6, 1)
  )
  context <- list(
    posterior_samples = samples,
    indicator_names   = character(),
    selection_spec    = list(jags_omega = "omega"),
    priors            = list(),
    flat_prior_list   = list(
      mu   = TRUE,
      PET  = TRUE,
      bias = BayesTools::prior_weightfunction(
        side    = "one-sided",
        steps   = .05,
        weights = BayesTools::wf_cumulative(c(1, 1))
      )
    ),
    formula_inputs    = list()
  )

  conditioning <- .iwmde_chen_conditioning_matrix(
    context        = context,
    parameter      = "mu",
    parameter_spec = list(type = "primitive", parameter = "mu"),
    active_rows    = 4L,
    weight_rows    = 1:3
  )

  expect_true(all(is.finite(conditioning[["fit"]])))
  expect_true(any(!is.finite(conditioning[["eval"]])))
})


test_that("Chen omega transforms use declared weight support", {

  bounded_context <- list(
    flat_prior_list = list(
      bias = BayesTools::prior_weightfunction(
        side    = "one-sided",
        steps   = .05,
        weights = BayesTools::wf_cumulative(c(1, 1))
      )
    ),
    selection_spec = list(jags_omega = "omega")
  )
  nonnegative_context <- list(
    flat_prior_list = list(
      bias = BayesTools::prior_weightfunction(
        side  = "one-sided",
        steps = .05,
        weights = BayesTools::wf_independent(
          BayesTools::prior(
            "gamma",
            parameters = list(shape = 2, rate = 1)
          )
        )
      )
    ),
    selection_spec = list(jags_omega = "omega")
  )

  bounded <- .iwmde_chen_transform_conditioning_column(
    context     = bounded_context,
    fit_values  = c(.25, .5),
    eval_values = 1,
    column      = "omega[2]"
  )
  nonnegative <- .iwmde_chen_transform_conditioning_column(
    context     = nonnegative_context,
    fit_values  = c(.25, .5),
    eval_values = 1,
    column      = "omega[2]"
  )

  expect_equal(bounded[["fit"]], stats::qlogis(c(.25, .5)))
  expect_identical(bounded[["eval"]], Inf)
  expect_equal(nonnegative[["fit"]], log1p(c(.25, .5)))
  expect_identical(nonnegative[["eval"]], log(2))
})
