test_that("invalid covariance blocks fail before zero-status classification", {

  for (value in c(NA_real_, NaN, Inf)) {
    expect_error(.block_covariance(list(1L), list(array(value, c(1L, 1L, 1L))), 1L, 1L),
                 "Block covariance values must be finite.", fixed = TRUE)
  }
})

test_that("parameter draw selection validates its container and extraction key", {

  expect_error(parameter_draws.brma(list(), 1),
    "'selection' must contain one resolved BayesTools parameter quantity.", fixed = TRUE)
  selection <- structure(list(quantities = data.frame(provider = "RoBMA")),
    class = "BayesTools_parameter_selection")
  expect_error(parameter_draws.brma(list(), selection),
    "RoBMA parameter quantity has no valid extraction key.", fixed = TRUE)
})

test_that("missing heterogeneity components do not erase all total draws", {

  expect_error(.total_brma_mv_heterogeneity_samples(list(matrix(1, 2L), NULL)),
    "Heterogeneity components must contain numeric samples.", fixed = TRUE)
  expect_equal(.total_brma_mv_heterogeneity_samples(list(matrix(3, 2L), matrix(4, 2L))),
                matrix(5, 2L))
})
