context("Selection boundary sentinels")

source(testthat::test_path("common-functions.R"))

skip_on_cran()


test_that("finite selected-normal CDF queries retain their numeric value", {

  skip_if_not(.has_native_selnorm_kernel())

  q        <- 1e290
  location <- 1e300
  scale    <- 1e300

  out <- .Call(
    "RoBMA_selnorm_kernel_cdf_matrix",
    as.numeric(q),
    matrix(location, nrow = 1L, ncol = 1L),
    matrix(scale, nrow = 1L, ncol = 1L),
    as.numeric(1),
    matrix(c(1, 1), nrow = 1L),
    as.numeric(0),
    as.integer(0),
    as.integer(SELKERNEL_STEP),
    as.numeric(c(0, -Inf)),
    as.numeric(c(Inf, 0)),
    as.integer(1),
    as.integer(0),
    as.numeric(c(0, 0)),
    as.numeric(c(0, 0)),
    as.numeric(c(-Inf, 0, Inf)),
    as.integer(c(2, 1)),
    as.integer(c(0, 0)),
    TRUE,
    TRUE,
    PACKAGE = "RoBMA"
  )

  expected <- stats::pnorm(q, mean = location, sd = scale)

  expect_equal(as.numeric(out), expected, tolerance = 8 * .Machine$double.eps)
  expect_lt(out, .5)
})


test_that("marginal-likelihood reconstruction restores only JAGS sentinels", {

  values <- c(-1e300, -1e290, 0, 1e290, 1e300)

  expect_identical(
    .marglik_restore_jags_selection_bounds(values),
    c(-Inf, -1e290, 0, 1e290, Inf)
  )

  data <- list(
    sel_kernel_mode          = SELKERNEL_STEP,
    sel_z_lower              = c(0, -1e300),
    sel_z_upper              = c(1e300, 0),
    sel_obs_bin              = 1L,
    sel_sign                 = 1L,
    sel_segment_bounds       = c(-1e300, 0, 1e300),
    sel_segment_step_bin     = c(2L, 1L),
    sel_segment_phack_region = c(0L, 0L),
    sel_omega                = c(1, .5)
  )
  selection <- .marglik_selection_context(parameters = list(), data = data)

  expect_identical(selection[["z_lower"]], c(0, -Inf))
  expect_identical(selection[["z_upper"]], c(Inf, 0))
  expect_identical(selection[["segments"]][["bounds"]], c(-Inf, 0, Inf))
})
