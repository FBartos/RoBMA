test_that("loo_weights rejects additional models with an actionable error", {

  fit_1 <- structure(list(), class = "brma")
  fit_2 <- structure(list(), class = "brma")

  expect_error(
    loo_weights(fit_1, fit_2),
    "loo_model_weights",
    fixed = TRUE
  )
})


.make_loo_model_weights_fixture <- function(shift, data_hash = "same-data") {

  draws <- seq(-1.5, 1.5, length.out = 200)
  y     <- c(-0.5, 0, 0.75, 1)
  log_lik <- vapply(y, function(observation) {
    stats::dnorm(observation, mean = draws + shift, sd = 1, log = TRUE)
  }, numeric(length(draws)))
  loo_result <- suppressWarnings(loo::loo(log_lik))
  attr(loo_result, "RoBMA_target") <- list(
    unit               = "estimate",
    conditioning_depth = "estimate",
    target             = "factorized_estimate",
    data_hash          = data_hash
  )

  return(loo_result)
}


test_that("loo_model_weights delegates stored brma LOO results", {

  loo_1 <- .make_loo_model_weights_fixture(0)
  loo_2 <- .make_loo_model_weights_fixture(0.3)
  fit_1 <- structure(list(stored_loo = loo_1), class = "brma")
  fit_2 <- structure(list(stored_loo = loo_2), class = "brma")

  testthat::local_mocked_bindings(
    loo.brma = function(x, unit = "estimate", ...) {
      expect_identical(unit, "estimate")
      x[["stored_loo"]]
    },
    .package = "RoBMA"
  )

  expected_pseudobma <- loo::loo_model_weights(
    list(fit_1 = loo_1, fit_2 = loo_2),
    method = "pseudobma",
    BB     = FALSE
  )
  expect_equal(
    loo_model_weights(fit_1, fit_2, method = "pseudobma", BB = FALSE),
    expected_pseudobma
  )

  expected_stacking <- loo::loo_model_weights(
    list(fit_1 = loo_1, fit_2 = loo_2),
    method = "stacking"
  )
  expect_equal(
    loo::loo_model_weights(fit_1, fit_2, method = "stacking"),
    expected_stacking
  )
})


test_that("loo_model_weights rejects incompatible stored targets", {

  fit_1 <- structure(
    list(stored_loo = .make_loo_model_weights_fixture(0)),
    class = "brma"
  )
  fit_2 <- structure(
    list(stored_loo = .make_loo_model_weights_fixture(0.3, "other-data")),
    class = "brma"
  )

  testthat::local_mocked_bindings(
    loo.brma = function(x, unit = "estimate", ...) x[["stored_loo"]],
    .package = "RoBMA"
  )

  expect_error(
    loo_model_weights(fit_1, fit_2),
    "different data",
    fixed = TRUE
  )
})
