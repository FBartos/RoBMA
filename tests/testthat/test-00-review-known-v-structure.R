test_that("Gaussian construction avoids selection-only exact factor recovery", {

  K <- 80L
  loading <- cbind(seq(.1, .3, length.out = K), rep(c(.2, .4), K / 2))
  V <- tcrossprod(loading) + diag(.2, K)
  calls <- 0L
  original <- .covariance_minimum_rank_factor
  testthat::local_mocked_bindings(
    .covariance_minimum_rank_factor = function(...) {
      calls <<- calls + 1L
      original(...)
    }, .package = "RoBMA"
  )
  ordinary <- brma.mv(yi = rep(0, K), V = V, random = NULL,
    prior_effect = prior("point", list(0)), measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE)
  expect_identical(calls, 0L)
  expect_identical(.known_v_covariance_matrix(.data_known_v_data(ordinary$data)), V)

  selected <- bselmodel.mv(yi = rep(0, K), V = V, random = NULL,
    prior_effect = prior("point", list(0)), measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE)
  expect_gt(calls, 0L)
  expect_true(.known_v_has_certified_factor(.data_known_v_data(selected$data)))
  expect_true(all(.data_selection_execution_plan(selected$data)$block_methods %in%
    c("factor", "rank_one")))
})

test_that("correlation levels preserve adjacent-gap membership within the cap", {

  expect_null(.covariance_correlation_levels(numeric(), .1, max_levels = 3L))
  expect_equal(.covariance_correlation_levels(c(.2, .2), .1, max_levels = 1L), .2)
  expect_equal(.covariance_correlation_levels(c(.1, .19, .28, .5), .1,
    max_levels = 2L), c(.19, .5))
  expect_equal(.covariance_correlation_levels(c(0, .25, .5), .25,
    max_levels = 1L), .25)
  expect_null(.covariance_correlation_levels(c(.2, .4, .7), .1, max_levels = 2L))
})
