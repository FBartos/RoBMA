context("as_draws methods for brma_samples objects")

# Load common test helpers
source(testthat::test_path("common-functions.R"))

# load prefitted model
skip_if_no_fits()
skip_if_not_installed("posterior")

fit <- load_fit("bcg_meta-analysis")

test_that("as_draws methods return posterior draws for brma_samples", {

  # obtain brma_samples
  brma_samples <- predict(fit)
  n_draws       <- nrow(brma_samples)
  n_chains      <- attr(brma_samples, "nchains")
  n_iter        <- attr(brma_samples, "niter")

  # test as_draws and consistency between methods
  draws <- RoBMA::as_draws(brma_samples)
  expect_s3_class(draws, "draws")
  expect_true(posterior::ndraws(draws)      == n_draws)
  expect_true(posterior::nchains(draws)     == n_chains)
  expect_true(posterior::niterations(draws) == n_iter)

  draws <- posterior::as_draws(brma_samples)
  expect_s3_class(draws, "draws")
  expect_true(posterior::ndraws(draws)      == n_draws)
  expect_true(posterior::nchains(draws)     == n_chains)
  expect_true(posterior::niterations(draws) == n_iter)

  # test as_draws_array
  draws <- RoBMA::as_draws_array(brma_samples)
  expect_s3_class(draws, "draws_array")
  expect_true(posterior::ndraws(draws)      == n_draws)
  expect_true(posterior::nchains(draws)     == n_chains)
  expect_true(posterior::niterations(draws) == n_iter)

  draws <- posterior::as_draws_array(brma_samples)
  expect_s3_class(draws, "draws_array")
  expect_true(posterior::ndraws(draws)      == n_draws)
  expect_true(posterior::nchains(draws)     == n_chains)
  expect_true(posterior::niterations(draws) == n_iter)

  # test as_draws_df
  draws <- RoBMA::as_draws_df(brma_samples)
  expect_s3_class(draws, "draws_df")
  expect_true(posterior::ndraws(draws)      == n_draws)
  expect_true(posterior::nchains(draws)     == n_chains)
  expect_true(posterior::niterations(draws) == n_iter)

  draws <- posterior::as_draws_df(brma_samples)
  expect_s3_class(draws, "draws_df")
  expect_true(posterior::ndraws(draws)      == n_draws)
  expect_true(posterior::nchains(draws)     == n_chains)
  expect_true(posterior::niterations(draws) == n_iter)

  # test as_draws_list
  draws <- RoBMA::as_draws_list(brma_samples)
  expect_s3_class(draws, "draws_list")
  expect_true(posterior::ndraws(draws)      == n_draws)
  expect_true(posterior::nchains(draws)     == n_chains)
  expect_true(posterior::niterations(draws) == n_iter)

  draws <- posterior::as_draws_list(brma_samples)
  expect_s3_class(draws, "draws_list")
  expect_true(posterior::ndraws(draws)      == n_draws)
  expect_true(posterior::nchains(draws)     == n_chains)
  expect_true(posterior::niterations(draws) == n_iter)

  # test as_draws_rvars
  draws <- RoBMA::as_draws_rvars(brma_samples)
  expect_s3_class(draws, "draws_rvars")
  expect_true(posterior::ndraws(draws)      == n_draws)
  expect_true(posterior::nchains(draws)     == n_chains)
  expect_true(posterior::niterations(draws) == n_iter)
  expect_true("mu" == posterior::variables(draws))

  draws <- posterior::as_draws_rvars(brma_samples)
  expect_s3_class(draws, "draws_rvars")
  expect_true(posterior::ndraws(draws)      == n_draws)
  expect_true(posterior::nchains(draws)     == n_chains)
  expect_true(posterior::niterations(draws) == n_iter)
  expect_true("mu" == posterior::variables(draws))
})

test_that("as_draws methods preserve brma_samples values and chain order", {

  brma_samples    <- predict(fit)
  expected_matrix <- unclass(brma_samples)
  draws_matrix    <- as.matrix(RoBMA::as_draws_matrix(brma_samples))
  selected_rows   <- unique(round(seq(1, nrow(expected_matrix), length.out = 7)))
  selected_vars   <- colnames(expected_matrix)

  expect_equal(
    as.vector(draws_matrix[selected_rows, selected_vars, drop = FALSE]),
    as.vector(expected_matrix[selected_rows, selected_vars, drop = FALSE]),
    tolerance = 0
  )

  draws_array <- RoBMA::as_draws_array(brma_samples)
  first_var   <- selected_vars[[1]]
  expect_equal(
    as.numeric(draws_array[1, 1, first_var]),
    as.numeric(expected_matrix[1, first_var]),
    tolerance = 0,
    info      = "first brma_samples chain and iteration are preserved"
  )
})
