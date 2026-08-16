test_that("zplot total SD matches replicated-SE evaluation", {

  tau_within <- matrix(
    c(0, .1, 1, 10, .25, .5, 2, 20, .75, 1.5, 3, 30),
    nrow = 4,
    ncol = 3
  )
  sei     <- c(.05, .2, 2)
  sei_mat <- matrix(sei, nrow = nrow(tau_within), ncol = ncol(tau_within),
                    byrow = TRUE)
  expected <- .root_sum_squares(tau_within, sei_mat)

  expect_identical(.zplot_total_sd(tau_within, sei), expected)
  expect_equal(
    .zplot_total_sd(tau_within, sei),
    sqrt(sweep(tau_within^2, 2, sei^2, "+")),
    tolerance = 1e-15
  )
})


test_that("zplot reuses invariant predictive components", {

  calls <- new.env(parent = emptyenv())
  calls$predictive <- 0L
  calls$density    <- 0L
  calls$paired     <- 0L
  calls$selection  <- FALSE

  predictive <- list(
    mu         = matrix(c(.1, .2, .3, .4), nrow = 2),
    tau_within = matrix(c(.05, .1, .15, .2), nrow = 2),
    sei        = c(.1, .2)
  )
  object <- list(fit = structure(list(), class = "BayesTools_fit"))

  testthat::local_mocked_bindings(
    .get_posterior_samples = function(...) matrix(0, nrow = 2, ncol = 1),
    .thin_sample_rows = function(...) NULL,
    .is_PET = function(...) FALSE,
    .is_PEESE = function(...) FALSE,
    .is_weightfunction = function(...) calls$selection,
    .effect_direction = function(...) "positive",
    .zplot_predictive_components = function(..., extrapolate) {
      calls$predictive <- calls$predictive + 1L
      if (extrapolate) {
        stop("Invariant predictive components were recomputed.")
      }
      predictive
    },
    .zplot_selection_context = function(...) {
      if (calls$selection) list(selection = TRUE) else NULL
    },
    .zplot_density_vectorized = function(...) {
      calls$density <- calls$density + 1L
      matrix(1, nrow = 1, ncol = 1)
    },
    .zplot_selnorm_density_pair = function(...) {
      calls$paired <- calls$paired + 1L
      list(
        fitted       = matrix(2, nrow = 1, ncol = 1),
        extrapolated = matrix(3, nrow = 1, ncol = 1)
      )
    },
    .package = "RoBMA"
  )

  ordinary <- .zplot_density_pair(object, z_sequence = 0, max_samples = 10)
  expect_identical(calls$predictive, 1L)
  expect_identical(calls$density, 1L)
  expect_identical(calls$paired, 0L)
  expect_identical(ordinary$fitted, ordinary$extrapolated)

  calls$predictive <- 0L
  calls$density    <- 0L
  calls$paired     <- 0L
  calls$selection  <- TRUE
  selected <- .zplot_density_pair(object, z_sequence = 0, max_samples = 10)
  expect_identical(calls$predictive, 1L)
  expect_identical(calls$density, 0L)
  expect_identical(calls$paired, 1L)
  expect_false(identical(selected$fitted, selected$extrapolated))
})


test_that("zplot retains separate PET and PEESE predictive components", {

  calls <- new.env(parent = emptyenv())
  calls$predictive <- 0L
  calls$density    <- 0L

  object <- list(fit = structure(list(), class = "BayesTools_fit"))
  testthat::local_mocked_bindings(
    .get_posterior_samples = function(...) matrix(0, nrow = 2, ncol = 1),
    .thin_sample_rows = function(...) NULL,
    .is_PET = function(...) TRUE,
    .is_PEESE = function(...) FALSE,
    .is_weightfunction = function(...) FALSE,
    .effect_direction = function(...) "positive",
    .zplot_predictive_components = function(..., extrapolate) {
      calls$predictive <- calls$predictive + 1L
      list(
        mu         = matrix(as.numeric(extrapolate), nrow = 2, ncol = 1),
        tau_within = matrix(.1, nrow = 2, ncol = 1),
        sei        = .2
      )
    },
    .zplot_selection_context = function(...) NULL,
    .zplot_density_vectorized = function(..., mu_samples) {
      calls$density <- calls$density + 1L
      mu_samples
    },
    .package = "RoBMA"
  )

  result <- .zplot_density_pair(object, z_sequence = 0, max_samples = 10)
  expect_identical(calls$predictive, 2L)
  expect_identical(calls$density, 2L)
  expect_equal(result$fitted, matrix(0, nrow = 2, ncol = 1))
  expect_equal(result$extrapolated, matrix(1, nrow = 2, ncol = 1))
})
