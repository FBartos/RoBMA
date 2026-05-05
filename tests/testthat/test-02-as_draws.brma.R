context("as_draws methods for brma objects")

# Load common test helpers
source(testthat::test_path("common-functions.R"))

skip_if_not_installed("posterior")

test_that("RoBMA as_draws defaults forward to posterior", {

  x <- matrix(seq_len(6), ncol = 2)
  colnames(x) <- c("alpha", "beta")

  expect_s3_class(RoBMA::as_draws(x), "draws")
  expect_s3_class(RoBMA::as_draws_array(x), "draws_array")
  expect_s3_class(RoBMA::as_draws_df(x), "draws_df")
  expect_s3_class(RoBMA::as_draws_list(x), "draws_list")
  expect_s3_class(RoBMA::as_draws_matrix(x), "draws_matrix")
  expect_s3_class(RoBMA::as_draws_rvars(x), "draws_rvars")
})

# load prefitted model
skip_if_no_fits()
fit <- load_fit("bcg_meta-analysis")

.fit_draw_dimensions <- function(fit) {

  mcmc_list <- RoBMA:::.brma_to_mcmc.list(fit)
  n_iter    <- vapply(mcmc_list, function(x) nrow(as.matrix(x)), integer(1))

  return(list(
    n_draws  = sum(n_iter),
    n_chains = length(mcmc_list),
    n_iter   = unique(n_iter)
  ))
}

fit_dims <- .fit_draw_dimensions(fit)

test_that("as_draws methods return posterior draws for fits", {

  # test as_draws and consistency between methods
  draws <- RoBMA::as_draws(fit)
  expect_s3_class(draws, "draws")
  expect_true(posterior::ndraws(draws)      == fit_dims[["n_draws"]])
  expect_true(posterior::nchains(draws)     == fit_dims[["n_chains"]])
  expect_true(posterior::niterations(draws) == fit_dims[["n_iter"]])

  draws <- posterior::as_draws(fit)
  expect_s3_class(draws, "draws")
  expect_true(posterior::ndraws(draws)      == fit_dims[["n_draws"]])
  expect_true(posterior::nchains(draws)     == fit_dims[["n_chains"]])
  expect_true(posterior::niterations(draws) == fit_dims[["n_iter"]])

  # test as_draws_array
  draws <- RoBMA::as_draws_array(fit)
  expect_s3_class(draws, "draws_array")
  expect_true(posterior::ndraws(draws)      == fit_dims[["n_draws"]])
  expect_true(posterior::nchains(draws)     == fit_dims[["n_chains"]])
  expect_true(posterior::niterations(draws) == fit_dims[["n_iter"]])

  draws <- posterior::as_draws_array(fit)
  expect_s3_class(draws, "draws_array")
  expect_true(posterior::ndraws(draws)      == fit_dims[["n_draws"]])
  expect_true(posterior::nchains(draws)     == fit_dims[["n_chains"]])
  expect_true(posterior::niterations(draws) == fit_dims[["n_iter"]])

  # test as_draws_df
  draws <- RoBMA::as_draws_df(fit)
  expect_s3_class(draws, "draws_df")
  expect_true(posterior::ndraws(draws)      == fit_dims[["n_draws"]])
  expect_true(posterior::nchains(draws)     == fit_dims[["n_chains"]])
  expect_true(posterior::niterations(draws) == fit_dims[["n_iter"]])

  draws <- posterior::as_draws_df(fit)
  expect_s3_class(draws, "draws_df")
  expect_true(posterior::ndraws(draws)      == fit_dims[["n_draws"]])
  expect_true(posterior::nchains(draws)     == fit_dims[["n_chains"]])
  expect_true(posterior::niterations(draws) == fit_dims[["n_iter"]])

  # test as_draws_list
  draws <- RoBMA::as_draws_list(fit)
  expect_s3_class(draws, "draws_list")
  expect_true(posterior::ndraws(draws)      == fit_dims[["n_draws"]])
  expect_true(posterior::nchains(draws)     == fit_dims[["n_chains"]])
  expect_true(posterior::niterations(draws) == fit_dims[["n_iter"]])

  draws <- posterior::as_draws_list(fit)
  expect_s3_class(draws, "draws_list")
  expect_true(posterior::ndraws(draws)      == fit_dims[["n_draws"]])
  expect_true(posterior::nchains(draws)     == fit_dims[["n_chains"]])
  expect_true(posterior::niterations(draws) == fit_dims[["n_iter"]])

  # test as_draws_rvars
  draws <- RoBMA::as_draws_rvars(fit)
  expect_s3_class(draws, "draws_rvars")
  expect_true(posterior::ndraws(draws)      == fit_dims[["n_draws"]])
  expect_true(posterior::nchains(draws)     == fit_dims[["n_chains"]])
  expect_true(posterior::niterations(draws) == fit_dims[["n_iter"]])
  expect_true("mu" %in% posterior::variables(draws))

  draws <- posterior::as_draws_rvars(fit)
  expect_s3_class(draws, "draws_rvars")
  expect_true(posterior::ndraws(draws)      == fit_dims[["n_draws"]])
  expect_true(posterior::nchains(draws)     == fit_dims[["n_chains"]])
  expect_true(posterior::niterations(draws) == fit_dims[["n_iter"]])
  expect_true("mu" %in% posterior::variables(draws))
})

test_that("as_draws methods preserve raw MCMC values and order", {

  mcmc_list       <- RoBMA:::.brma_to_mcmc.list(fit)
  expected_matrix <- do.call(rbind, lapply(mcmc_list, as.matrix))
  draws_matrix    <- as.matrix(RoBMA::as_draws_matrix(fit))
  variables       <- intersect(colnames(expected_matrix), colnames(draws_matrix))
  selected_rows   <- unique(round(seq(1, nrow(expected_matrix), length.out = 7)))
  selected_vars   <- head(variables, 8)

  expect_true(length(selected_vars) > 0L)
  expect_equal(
    as.vector(draws_matrix[selected_rows, selected_vars, drop = FALSE]),
    as.vector(expected_matrix[selected_rows, selected_vars, drop = FALSE]),
    tolerance = 0
  )

  draws_array <- RoBMA::as_draws_array(fit)
  first_var   <- selected_vars[[1]]
  expect_equal(
    as.numeric(draws_array[1, 1, first_var]),
    as.numeric(as.matrix(mcmc_list[[1]])[1, first_var]),
    tolerance = 0,
    info      = "first chain and iteration are preserved"
  )
})

test_that("as_draws methods preserve product-space BMA and RoBMA indicators", {

  product_names <- c(
    "dat.lehmann2018_BMA.norm",
    "bcg_BMA.glmm",
    "nielweise2008_BMA.glmm",
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
