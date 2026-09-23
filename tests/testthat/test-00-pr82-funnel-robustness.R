test_that("common funnel heterogeneity tolerates only double roundoff", {

  tau <- matrix(c(1, 2, 1 + .Machine$double.eps, 2), nrow = 2L)
  testthat::local_mocked_bindings(
    .get_posterior_samples = function(...) matrix(0, 2L, 1L),
    .funnel_row_heterogeneity_samples = function(...) tau,
    .funnel_selection_source_variances = function(...) NULL,
    .package = "RoBMA"
  )
  expect_true(.funnel_common_heterogeneity(list())$common)
  tau[1L, 2L] <- 1 + 1e-8
  expect_false(.funnel_common_heterogeneity(list())$common)
  tau[1L, 2L] <- NA_real_
  expect_error(.funnel_common_heterogeneity(list()),
               "finite row-marginal heterogeneity", fixed = TRUE)
})

test_that("conditioned funnel quadrature retains zero-width mixture branches", {

  expect_identical(.funnel_conditioned_order(c(1, 0), c(0.1, 1)), 401L)
  expect_lt(.funnel_conditioned_order(c(1, 0), c(0.1, 0)), 401L)
})

test_that("funnel clipping handles singleton edges and rejects invalid coordinates", {

  expect_equal(.clip_line_x(2, 1, c(0, 1)), data.frame(x = numeric(), y = numeric()))
  expect_equal(.clip_line_x(0.5, 1, c(0, 1)), data.frame(x = 0.5, y = 1))
  expect_equal(.clip_line_x(c(-1, 2), c(0, 3), c(0, 1)),
               data.frame(x = c(0, 1), y = c(1, 2)))
  expect_error(.clip_line_x(c(0, NA), c(0, 1), c(0, 1)),
               "Funnel contour coordinates must be finite.", fixed = TRUE)
})

test_that("funnel data rejects missing standard errors before constructing axes", {

  testthat::local_mocked_bindings(
    .effect_direction = function(...) "positive",
    .outcome_type = function(...) "norm",
    .outcome_data_yi = function(...) c(0, 1),
    .outcome_data_sei = function(...) c(0.1, NA_real_),
    .package = "RoBMA"
  )
  expect_error(
    .funnel_data_outcome(list(), TRUE, TRUE, 10, "plugin", list(), list()),
    "Funnel-plot standard errors must be finite and non-negative.", fixed = TRUE
  )
})
