context("as_draws methods for brma objects")

# Load common test helpers
source(testthat::test_path("common-functions.R"))

# load prefitted model
skip_if_no_fits()
fit <- load_fit("bcg_meta-analysis")


test_that("as_draws methods work for brma objects", {

  skip_if_not_installed("posterior")

  # test as_draws
  draws <- as_draws(fit)
  expect_s3_class(draws, "draws")
  expect_true(posterior::ndraws(draws) > 0)

  # test as_draws_array
  draws_array <- as_draws_array(fit)
  expect_s3_class(draws_array, "draws_array")
  expect_equal(length(dim(draws_array)), 3)  # iteration x chain x variable

  # test as_draws_df
  draws_df <- as_draws_df(fit)
  expect_s3_class(draws_df, "draws_df")
  expect_true(is.data.frame(draws_df))
  expect_true(".iteration" %in% colnames(draws_df))
  expect_true(".chain" %in% colnames(draws_df))

  # test as_draws_list
  draws_list <- as_draws_list(fit)
  expect_s3_class(draws_list, "draws_list")
  expect_true(is.list(draws_list))

  # test as_draws_matrix
  draws_matrix <- as_draws_matrix(fit)
  expect_s3_class(draws_matrix, "draws_matrix")
  expect_equal(length(dim(draws_matrix)), 2)  # draw x variable

  # test as_draws_rvars
  draws_rvars <- as_draws_rvars(fit)
  expect_s3_class(draws_rvars, "draws_rvars")

  # verify that key parameters are present in the draws
  var_names <- posterior::variables(draws)
  expect_true("mu" %in% var_names)
  expect_true("tau" %in% var_names)

})
