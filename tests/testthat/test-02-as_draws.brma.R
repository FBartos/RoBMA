context("as_draws methods for brma objects")

# Load common test helpers
source(testthat::test_path("common-functions.R"))

# load prefitted model
skip_if_no_fits()
skip_if_not_installed("posterior")
fit <- load_fit("bcg_meta-analysis")

test_that("as_draws methods work for fit", {

  # test as_draws and consistency between methods
  draws <- RoBMA::as_draws(fit)
  expect_s3_class(draws, "draws")
  expect_true(posterior::ndraws(draws)      == 15000)
  expect_true(posterior::nchains(draws)     == 3)
  expect_true(posterior::niterations(draws) == 5000)

  draws <- posterior::as_draws(fit)
  expect_s3_class(draws, "draws")
  expect_true(posterior::ndraws(draws)      == 15000)
  expect_true(posterior::nchains(draws)     == 3)
  expect_true(posterior::niterations(draws) == 5000)

  # test as_draws_array
  draws <- RoBMA::as_draws_array(fit)
  expect_s3_class(draws, "draws_array")
  expect_true(posterior::ndraws(draws)      == 15000)
  expect_true(posterior::nchains(draws)     == 3)
  expect_true(posterior::niterations(draws) == 5000)

  draws <- posterior::as_draws_array(fit)
  expect_s3_class(draws, "draws_array")
  expect_true(posterior::ndraws(draws)      == 15000)
  expect_true(posterior::nchains(draws)     == 3)
  expect_true(posterior::niterations(draws) == 5000)

  # test as_draws_df
  draws <- RoBMA::as_draws_df(fit)
  expect_s3_class(draws, "draws_df")
  expect_true(posterior::ndraws(draws)      == 15000)
  expect_true(posterior::nchains(draws)     == 3)
  expect_true(posterior::niterations(draws) == 5000)

  draws <- posterior::as_draws_df(fit)
  expect_s3_class(draws, "draws_df")
  expect_true(posterior::ndraws(draws)      == 15000)
  expect_true(posterior::nchains(draws)     == 3)
  expect_true(posterior::niterations(draws) == 5000)

  # test as_draws_list
  draws <- RoBMA::as_draws_list(fit)
  expect_s3_class(draws, "draws_list")
  expect_true(posterior::ndraws(draws)      == 15000)
  expect_true(posterior::nchains(draws)     == 3)
  expect_true(posterior::niterations(draws) == 5000)

  draws <- posterior::as_draws_list(fit)
  expect_s3_class(draws, "draws_list")
  expect_true(posterior::ndraws(draws)      == 15000)
  expect_true(posterior::nchains(draws)     == 3)
  expect_true(posterior::niterations(draws) == 5000)

  # test as_draws_rvars
  draws <- RoBMA::as_draws_rvars(fit)
  expect_s3_class(draws, "draws_rvars")
  expect_true(posterior::ndraws(draws)      == 15000)
  expect_true(posterior::nchains(draws)     == 3)
  expect_true(posterior::niterations(draws) == 5000)
  expect_true("mu" %in% posterior::variables(draws))

  draws <- posterior::as_draws_rvars(fit)
  expect_s3_class(draws, "draws_rvars")
  expect_true(posterior::ndraws(draws)      == 15000)
  expect_true(posterior::nchains(draws)     == 3)
  expect_true(posterior::niterations(draws) == 5000)
  expect_true("mu" %in% posterior::variables(draws))
})
