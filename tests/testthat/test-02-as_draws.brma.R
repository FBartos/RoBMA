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

  mcmc_list <- RoBMA:::.brma_to_mcmc.list(
    fit,
    include_auxiliary = TRUE
  )
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

  mcmc_list       <- RoBMA:::.brma_to_mcmc.list(
    fit,
    include_auxiliary = TRUE
  )
  expected_matrix <- do.call(rbind, lapply(mcmc_list, as.matrix))
  draws_matrix    <- as.matrix(RoBMA::as_draws_matrix(
    fit,
    include_auxiliary = TRUE
  ))
  variables       <- intersect(colnames(expected_matrix), colnames(draws_matrix))
  selected_rows   <- unique(round(seq(1, nrow(expected_matrix), length.out = 7)))
  selected_vars   <- head(variables, 8)

  expect_true(length(selected_vars) > 0L)
  expect_equal(
    as.vector(draws_matrix[selected_rows, selected_vars, drop = FALSE]),
    as.vector(expected_matrix[selected_rows, selected_vars, drop = FALSE]),
    tolerance = 0
  )

  draws_array <- RoBMA::as_draws_array(
    fit,
    include_auxiliary = TRUE
  )
  first_var   <- selected_vars[[1]]
  expect_equal(
    as.numeric(draws_array[1, 1, first_var]),
    as.numeric(as.matrix(mcmc_list[[1]])[1, first_var]),
    tolerance = 0,
    info      = "first chain and iteration are preserved"
  )
})

test_that("auxiliary catalog separates backend state from model parameters", {

  variables <- c(
    "mu",
    "tau_indicator",
    "theta[1]",
    "gamma[1]",
    "sampling_z[1]",
    "mu__xREx__study_xRE_Zx[1,1]",
    "mu__xREx__study_xRE_CORx_L[1,1]",
    "mu__xREx__study_xRE_CORx_lkj_u[1]",
    "mu__xREx__study_xRE_CORx_lkj_cpc[1]",
    "mu__xREx__study_xRE_CORx_R[1,1]",
    "prior_par_eta_mu__xRE_ALLOCx_total__weight[1]",
    "inv_tau",
    "eta[1]",
    "log_omega[1]",
    "alpha",
    "pi_null",
    "mu__xRE_ALLOCx_total__weight[1]",
    "mu_variable",
    "mu_beta",
    "mu_beta_indicator",
    "mu_beta_variable",
    "mu_beta_inclusion"
  )
  fit_metadata <- list()
  attr(fit_metadata, "prior_list") <- list(
    mu_beta = BayesTools::prior_spike_and_slab(
      BayesTools::prior(
        distribution = "normal",
        parameters   = list(mean = 0, sd = 1)
      )
    )
  )
  object <- structure(list(fit = fit_metadata), class = "brma")
  catalog <- RoBMA:::.brma_auxiliary_variable_catalog(object, variables)

  expect_equal(
    catalog[["variable"]],
    variables[c(3:9, 11:16, 21:22)]
  )
  expect_false("tau_indicator" %in% catalog[["variable"]])
  expect_false("mu__xREx__study_xRE_CORx_R[1,1]" %in% catalog[["variable"]])
  expect_false("mu__xRE_ALLOCx_total__weight[1]" %in% catalog[["variable"]])
  expect_false("mu_variable" %in% catalog[["variable"]])
  expect_false("mu_beta_indicator" %in% catalog[["variable"]])
  expect_true(all(nzchar(catalog[["category"]])))
})

test_that("default draw extraction filters auxiliary variables only", {

  values <- matrix(
    seq_len(32),
    nrow     = 4,
    dimnames = list(NULL, c(
      "mu",
      "tau_indicator",
      "sampling_z[1]",
      "theta[1]",
      "mu__xREx__study_xRE_Zx[1,1]",
      "mu__xREx__study_xRE_CORx_lkj_u[1]",
      "mu__xREx__study_xRE_CORx_lkj_cpc[1]",
      "mu__xREx__study_xRE_CORx_R[1,1]"
    ))
  )
  chain <- coda::mcmc(values, start = 5, thin = 2)
  materialized <- coda::mcmc.list(chain)
  fit <- structure(list(sentinel = TRUE), class = "BayesTools_fit")
  attr(fit, "prior_list") <- list()
  object <- structure(list(fit = fit), class = "brma")
  testthat::local_mocked_bindings(
    JAGS_validate_fit_contract = function(...) invisible(TRUE),
    JAGS_materialize_draws = function(...) materialized,
    .package = "BayesTools"
  )

  filtered <- RoBMA:::.brma_to_mcmc.list(object)
  raw      <- RoBMA:::.brma_to_mcmc.list(
    object,
    include_auxiliary = TRUE
  )

  expect_equal(
    colnames(as.matrix(filtered[[1]])),
    c("mu", "tau_indicator", "mu__xREx__study_xRE_CORx_R[1,1]")
  )
  expect_equal(coda::mcpar(filtered[[1]]), coda::mcpar(chain))
  expect_equal(raw, materialized, tolerance = 0)
})

test_that("all brma draw formats share auxiliary filtering and validation", {

  converters <- list(
    RoBMA::as_draws,
    RoBMA::as_draws_array,
    RoBMA::as_draws_df,
    RoBMA::as_draws_list,
    RoBMA::as_draws_matrix,
    RoBMA::as_draws_rvars
  )
  expected_default <- posterior::variables(RoBMA::as_draws_matrix(fit))
  expected_raw <- colnames(as.matrix(
    RoBMA:::.brma_to_mcmc.list(fit, include_auxiliary = TRUE)[[1]]
  ))

  for (converter in converters) {
    expect_equal(
      posterior::variables(converter(fit)),
      expected_default
    )
    expect_equal(
      posterior::variables(converter(fit, include_auxiliary = TRUE)),
      expected_raw
    )
    expect_error(
      converter(fit, include_auxiliary = 1),
      "include_auxiliary"
    )
  }
})

test_that("all raw draw formats preserve auxiliary values", {

  name <- "brma.mv_latent"
  skip_if_missing_fits(name)

  object     <- load_fit(name, validate = FALSE)
  converters <- list(
    RoBMA::as_draws,
    RoBMA::as_draws_array,
    RoBMA::as_draws_df,
    RoBMA::as_draws_list,
    RoBMA::as_draws_matrix,
    RoBMA::as_draws_rvars
  )
  raw_mcmc <- RoBMA:::.brma_to_mcmc.list(
    object,
    include_auxiliary = TRUE
  )
  expected <- do.call(rbind, lapply(raw_mcmc, as.matrix))

  expect_true(any(startsWith(colnames(expected), "sampling_z[")))
  for (converter in converters) {
    actual <- as.matrix(posterior::as_draws_matrix(converter(
      object,
      include_auxiliary = TRUE
    )))

    expect_equal(colnames(actual), colnames(expected))
    expect_equal(as.vector(actual), as.vector(expected), tolerance = 0)
  }
})

test_that("draw schemas are stable across known-V parameterizations", {

  names <- c(
    "brma.mv_latent",
    "brma.mv_whitened",
    "brma.mv_block_mvn",
    "brma.mv_block_mvn_random",
    "brma.mv_block_mvn_random_sampled",
    "brma.mv_block_mvn_known_R"
  )
  skip_if_missing_fits(names)

  objects <- lapply(names, load_fit, validate = FALSE)
  names(objects) <- names
  schemas <- lapply(objects, function(object) {
    posterior::variables(RoBMA::as_draws_matrix(object))
  })
  raw_schemas <- lapply(objects, function(object) {
    posterior::variables(RoBMA::as_draws_matrix(
      object,
      include_auxiliary = TRUE
    ))
  })

  expect_equal(schemas[["brma.mv_latent"]], c("mu", "tau"))
  expect_equal(schemas[["brma.mv_whitened"]], c("mu", "tau"))
  expect_equal(schemas[["brma.mv_block_mvn"]], c("mu", "tau"))
  expect_equal(
    schemas[["brma.mv_block_mvn_random"]],
    schemas[["brma.mv_block_mvn_random_sampled"]]
  )
  expect_true(any(startsWith(
    raw_schemas[["brma.mv_latent"]],
    "sampling_z["
  )))
  expect_false(any(startsWith(
    schemas[["brma.mv_latent"]],
    "sampling_z["
  )))
  expect_true(any(grepl(
    "_xRE_Zx[",
    raw_schemas[["brma.mv_block_mvn_random_sampled"]],
    fixed = TRUE
  )))
  expect_false(any(grepl(
    "_xRE_Zx[",
    schemas[["brma.mv_block_mvn_random_sampled"]],
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "_intercept$",
    schemas[["brma.mv_block_mvn_known_R"]]
  )))
})

test_that("draw schemas filter model-family-specific auxiliary state", {

  names <- c(
    "konstantopoulos2011_3lvl",
    "bcg_glmm_reg",
    "dat.lehmann2018-3PSM",
    "brma.mv_block_mvn_random_mods_scale",
    "brma.mv_v14_konstantopoulos2011_cs",
    "brma.mv_v14_ishak2007_har"
  )
  skip_if_missing_fits(names)

  schemas <- lapply(names, function(name) {
    posterior::variables(RoBMA::as_draws_matrix(
      load_fit(name, validate = FALSE)
    ))
  })
  names(schemas) <- names

  expect_false(any(startsWith(schemas[["konstantopoulos2011_3lvl"]], "gamma[")))
  expect_true("rho" %in% schemas[["konstantopoulos2011_3lvl"]])
  expect_false(any(startsWith(schemas[["bcg_glmm_reg"]], "theta[")))
  expect_true(any(startsWith(schemas[["bcg_glmm_reg"]], "pi[")))
  expect_false(any(startsWith(schemas[["dat.lehmann2018-3PSM"]], "eta[")))
  expect_true(any(startsWith(schemas[["dat.lehmann2018-3PSM"]], "omega[")))

  for (name in names[4:6]) {
    expect_false(any(grepl("_xRE_Zx[", schemas[[name]], fixed = TRUE)))
    expect_false(any(grepl("_xRE_CORx_L[", schemas[[name]], fixed = TRUE)))
    expect_false(any(startsWith(schemas[[name]], "prior_par_eta_")))
  }
  expect_true(any(grepl(
    "_xRE_CORx_R[",
    schemas[["brma.mv_v14_konstantopoulos2011_cs"]],
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "_xRE_CORx_R[",
    schemas[["brma.mv_v14_ishak2007_har"]],
    fixed = TRUE
  )))
})

test_that("known-V backend posterior draws agree for common parameters", {

  names <- c(
    "brma.mv_latent",
    "brma.mv_whitened",
    "brma.mv_block_mvn"
  )
  skip_if_missing_fits(names)

  draws <- lapply(names, function(name) {
    as.matrix(RoBMA::as_draws_matrix(
      load_fit(name, validate = FALSE)
    ))
  })
  estimates <- vapply(draws, function(x) {
    c(mu = mean(x[, "mu"]), tau = mean(x[, "tau"]))
  }, numeric(2))

  expect_lt(diff(range(estimates["mu", ])), 0.08)
  expect_lt(diff(range(estimates["tau", ])), 0.08)
})

test_that("equivalent brma and brma.mv fits expose the same model schema", {

  names <- c("vif_parity_brma", "vif_parity_brma_mv")
  skip_if_missing_fits(names)

  schemas <- lapply(names, function(name) {
    posterior::variables(RoBMA::as_draws_matrix(
      load_fit(name, validate = FALSE)
    ))
  })

  expect_equal(schemas[[1]], schemas[[2]])
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
