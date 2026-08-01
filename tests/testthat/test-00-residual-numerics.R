test_that("LOO-PIT transformation preserves extreme normal quantiles", {

  z <- c(-100, -40, 0, 40, 100)
  log_lower <- matrix(stats::pnorm(z, log.p = TRUE), nrow = 1L)
  log_upper <- matrix(stats::pnorm(
    z,
    lower.tail = FALSE,
    log.p      = TRUE
  ), nrow = 1L)

  expect_equal(
    .loo_pit_z_from_log_tails(
      log_lower    = log_lower,
      log_upper    = log_upper,
      psis_weights = matrix(1, nrow = 1L, ncol = length(z))
    ),
    z,
    tolerance = 1e-12
  )
})


test_that("normal LOO-PIT tails remain finite after probability underflow", {

  setup <- list(
    outcome_type = "norm",
    data         = list(),
    yi           = c(-100, 100),
    sei          = c(1, 1),
    S            = 2L,
    K            = 2L,
    mu           = matrix(0, nrow = 2L, ncol = 2L),
    tau_within   = matrix(0, nrow = 2L, ncol = 2L)
  )
  tails <- .loo_predictive_log_tails_estimate(setup)

  expect_true(all(is.finite(tails[["log_lower"]])))
  expect_true(all(is.finite(tails[["log_upper"]])))
  expect_equal(
    .loo_pit_z_from_log_tails(
      tails[["log_lower"]],
      tails[["log_upper"]],
      matrix(0.5, nrow = 2L, ncol = 2L)
    ),
    c(-100, 100),
    tolerance = 1e-12
  )
})


test_that("GLMM LOO-PIT is rejected before PSIS work", {

  expect_error(
    .check_residual_type_availability("LOO-PIT", "bin", FALSE),
    "discrete PIT convention"
  )
  expect_error(
    .loo_predictive_log_tails_estimate(list(outcome_type = "pois")),
    "discrete PIT convention"
  )
})
