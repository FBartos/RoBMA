context("as_draws methods for brma objects")

# Load common test helpers
source(testthat::test_path("common-functions.R"))

# load prefitted model
skip_if_no_fits()
skip_if_not_installed("posterior")
fit <- load_fit("bcg_meta-analysis")

test_that("as_draws methods return posterior draws for fits", {

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

test_that("as_draws methods preserve product-space BMA and RoBMA indicators", {

  product_names <- c(
    "dat.lehmann2018_BMA.norm",
    "bcg_BMA.glmm",
    "dat.lehmann2018_RoBMA"
  )
  skip_if_missing_fits(product_names)

  for (name in product_names) {

    fit_product <- load_fit(name, validate = FALSE)
    n_chains    <- length(fit_product[["fit"]][["mcmc"]])
    n_iter      <- fit_product[["fit"]][["sample"]]
    n_draws     <- n_chains * n_iter

    draws <- RoBMA::as_draws(fit_product)
    expect_s3_class(draws, "draws")
    expect_equal(posterior::ndraws(draws), n_draws)
    expect_equal(posterior::nchains(draws), n_chains)
    expect_equal(posterior::niterations(draws), n_iter)

    draws_array <- RoBMA::as_draws_array(fit_product)
    expect_s3_class(draws_array, "draws_array")
    expect_equal(posterior::ndraws(draws_array), n_draws)
    expect_equal(posterior::nchains(draws_array), n_chains)
    expect_equal(posterior::niterations(draws_array), n_iter)

    draws_df <- RoBMA::as_draws_df(fit_product)
    expect_s3_class(draws_df, "draws_df")
    expect_equal(posterior::ndraws(draws_df), n_draws)
    expect_equal(posterior::nchains(draws_df), n_chains)
    expect_equal(posterior::niterations(draws_df), n_iter)

    draws_matrix <- RoBMA::as_draws_matrix(fit_product)
    expect_s3_class(draws_matrix, "draws_matrix")
    expect_equal(posterior::ndraws(draws_matrix), n_draws)

    draws_rvars <- RoBMA::as_draws_rvars(fit_product)
    expect_s3_class(draws_rvars, "draws_rvars")
    expect_equal(posterior::ndraws(draws_rvars), n_draws)

    variables <- posterior::variables(draws_df)
    expect_true(any(grepl("_indicator$", variables)),
                info = paste(name, "exposes product-space indicators"))

    if (!is.null(fit_product[["priors"]][["outcome"]][["mu"]]) &&
        BayesTools::is.prior.mixture(fit_product[["priors"]][["outcome"]][["mu"]])) {
      expect_true("mu_indicator" %in% variables, info = name)
    }
    if (!is.null(fit_product[["priors"]][["outcome"]][["tau"]]) &&
        BayesTools::is.prior.mixture(fit_product[["priors"]][["outcome"]][["tau"]])) {
      expect_true("tau_indicator" %in% variables, info = name)
    }
    if (!is.null(fit_product[["priors"]][["outcome"]][["bias"]])) {
      expect_true("bias_indicator" %in% variables, info = name)
    }
  }
})
