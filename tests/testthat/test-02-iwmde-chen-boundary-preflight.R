# ============================================================================ #
# IWMDE Chen boundary preflight
# ============================================================================ #


test_that("Chen boundary preflight preserves fallback without full construction", {

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
    flat_prior_list   = list(mu = TRUE, PET = TRUE, rho = TRUE),
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
    c("rho", "rho", "omega[1]", "omega[1]", "omega[2]", "omega[2]")
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
    flat_prior_list   = list(mu = TRUE, PET = TRUE),
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
