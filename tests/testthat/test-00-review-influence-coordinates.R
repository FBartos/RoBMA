context("Influence diagnostics use the fitted coordinates of factor levels")

test_that("COVRATIO of a dependent factor-level parameterization uses the fitted coordinates", {

  # The levels of a mean-difference factor are linearly dependent (they sum to
  # zero), so their covariance is singular; the fitted coordinates are not.
  n      <- 100L
  object <- review_factor_replay(c("A", "B"), "meandif", independent = TRUE, n = n)
  weights <- matrix(1 / n, nrow = n, ncol = 6L)

  levels <- .diagnostic_location_parameter_samples(object)
  expect_equal(qr(scale(levels, center = TRUE, scale = FALSE))[["rank"]], 2L)
  expect_equal(ncol(levels), 3L)

  observed <- suppressWarnings(covratio(object, .weights = weights))
  coordinates <- .diagnostic_location_parameter_samples(
    object, transform_factors = FALSE
  )
  expected <- .covratio_internal(coordinates, weights)

  expect_true(all(is.finite(unclass(observed))))
  expect_equal(as.numeric(observed), expected[["values"]], tolerance = 1e-12)
  # equal weights leave the posterior covariance unchanged
  expect_equal(as.numeric(observed), rep(1, 6L), tolerance = 1e-10)
  expect_null(attr(observed, "note"))
})
