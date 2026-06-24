# ============================================================================ #
#
# Test bridge sampling marginal likelihood computation for brma objects.
# This tests the logml, bf, bayes_factor, and post_prob S3 methods.
# The bridge_sampler function is tested directly alongside the fit functions
#
# ============================================================================ #
context("bridgesampling methods for brma objects")

source(testthat::test_path("common-functions.R"))
skip_on_cran()
skip_if_no_fits()
skip_if_not_installed("metafor")

# list cached marginal-likelihood fits lazily
marglik_names <- list_fits(has_marglik = TRUE)
fits          <- lazy_fits(marglik_names, validate = FALSE)
info          <- lazy_infos(marglik_names, validate = FALSE)

# ---------------------------------------------------------------------------- #
# bridge_sampler function works with all model types
# ---------------------------------------------------------------------------- #

test_that("bridge_sampler extracts bridge sampling object", {

  for (name in marglik_names) {
    # the marginal likelihood inherits the correct class
    expect_s3_class(bridge_sampler(fits[[name]]), "bridge")
  }
})

test_that("add_marglik computes bridge sampling for brma.mv known-V fits", {

  mv_names <- c(
    "brma.mv_latent",
    "brma.mv_whitened",
    "brma.mv_block_mvn",
    "brma.mv_block_mvn_fixed_random_null",
    "brma.mv_block_mvn_random",
    "brma.mv_block_mvn_known_R",
    "brma.mv_latent_estimate_scale",
    "brma.mv_block_mvn_estimate_scale"
  )
  skip_if_missing_fits(mv_names)

  same_model_names <- c(
    "brma.mv_latent",
    "brma.mv_whitened",
    "brma.mv_block_mvn"
  )
  estimate_scale_names <- c(
    "brma.mv_latent_estimate_scale",
    "brma.mv_block_mvn_estimate_scale"
  )
  backend_logml        <- numeric()
  estimate_scale_logml <- numeric()
  same_models          <- list()
  for (name in mv_names) {
    set.seed(100)
    fit_brma <- add_marglik(load_fit(name, validate = FALSE))
    target   <- attr(fit_brma[["marglik"]], "RoBMA_target", exact = TRUE)
    expect_true(inherits(fit_brma[["marglik"]], "bridge"), info = name)
    expect_true(is.finite(logml(fit_brma)), info = name)
    expect_equal(target[["reported_target"]], "full joint fitted likelihood", info = name)
    expect_equal(
      target[["known_v_parameterization"]],
      .data_known_v_data(fit_brma[["data"]])[["parameterization"]],
      info = name
    )
    expect_equal(
      target[["known_v_parameterization_requested"]],
      .data_known_v_data(fit_brma[["data"]])[["parameterization_requested"]],
      info = name
    )
    if (name %in% same_model_names) {
      backend_logml[[name]] <- logml(fit_brma)
      same_models[[name]]   <- fit_brma
    }
    if (name %in% estimate_scale_names) {
      estimate_scale_logml[[name]] <- logml(fit_brma)
    }
  }
  expect_lt(diff(range(backend_logml)), 0.50)
  expect_lt(diff(range(estimate_scale_logml)), 0.75)

  bf_result <- bf(
    same_models[["brma.mv_latent"]],
    same_models[["brma.mv_whitened"]]
  )
  bf_alias <- bayes_factor(
    same_models[["brma.mv_latent"]],
    same_models[["brma.mv_block_mvn"]]
  )
  pp_result <- post_prob(
    same_models[["brma.mv_latent"]],
    same_models[["brma.mv_whitened"]],
    same_models[["brma.mv_block_mvn"]],
    model_names = names(same_models)
  )

  expect_s3_class(bf_result, "bf_default")
  expect_s3_class(bf_alias, "bf_default")
  expect_true(is.finite(bf_result[["bf"]]))
  expect_true(is.finite(bf_alias[["bf"]]))
  expect_equal(length(pp_result), length(same_models))
  expect_true(all(is.finite(pp_result)))
  expect_equal(sum(pp_result), 1, tolerance = sqrt(.Machine$double.eps))
})

test_that("add_marglik computes bridge sampling for marginalized random allocation", {

  name <- "brma.mv_block_mvn_random_scale"
  skip_if_missing_fits(name)

  set.seed(100)
  fit_brma <- add_marglik(load_fit(name, validate = FALSE))
  target   <- attr(fit_brma[["marglik"]], "RoBMA_target", exact = TRUE)
  expect_true(inherits(fit_brma[["marglik"]], "bridge"))
  expect_true(is.finite(logml(fit_brma)))
  expect_equal(target[["reported_target"]], "full joint fitted likelihood")
})

test_that("v14 brma.mv metafor fixtures cache usable marginal likelihoods", {

  mv_names <- c(
    "brma.mv_v14_konstantopoulos2011_cs",
    "brma.mv_v14_assink2016_nested",
    "brma.mv_v14_ishak2007_har",
    "brma.mv_v14_begg1989_study_treatment"
  )
  skip_if_missing_fits(mv_names)

  for (name in mv_names) {
    fit_brma <- fits[[name]]
    bridge   <- bridge_sampler(fit_brma)
    target   <- attr(bridge, "RoBMA_target", exact = TRUE)

    expect_true(inherits(bridge, "bridge"), info = name)
    expect_true(is.finite(logml(fit_brma)), info = name)
    expect_true(is.finite(bridge[["mcse_logml"]]), info = name)
    expect_true(bridge[["mcse_logml"]] < 1, info = name)
    expect_equal(target[["reported_target"]], "full joint fitted likelihood",
                 info = name)
    expect_true(isTRUE(target[["known_v"]]), info = name)
    expect_equal(
      target[["known_v_parameterization"]],
      .data_known_v_data(fit_brma[["data"]])[["parameterization"]],
      info = name
    )
    expect_equal(
      target[["known_v_parameterization_requested"]],
      .data_known_v_data(fit_brma[["data"]])[["parameterization_requested"]],
      info = name
    )
  }
})

test_that("bridge sampling marginal likelihood is close to BIC for metafor-reference fits", {

  metafor_names <- list_fits(has_marglik = TRUE, has_metafor = TRUE)

  for (name in metafor_names) {

    # the marginal likelihood and BIC-based marglik are close
    # (the analyses don't use unit information priors but the
    #  difference should not be large)
    if (name %in% c("bcg_glmm", "nielweise2008_glmm", "bcg_glmm_reg")) {
      next
    } # skip glmm because of marginalization differences

    if (name %in% c("bcg_meta-regression", "bcg_meta-regression2", "bcg_meta-regression3", "bcg_meta-regression3b", "bcg_meta-regression4", "bcg_meta-regression4b")) {
      next
    } # skip because of scaled priors differences
    if (inherits(fits[[name]], "brma.mv")) {
      next
    } # skip because brma.mv marginal likelihood targets the full known-V joint likelihood

    BIC_metafor <- BIC(info[[name]][["metafor"]])
    marglik     <- logml(fits[[name]])

    expect_equal(-BIC_metafor / 2, marglik, tolerance = 0.15)
  }
})


test_that("effect-direction reversals preserve marginal likelihood", {

  fit1 <- fits[["dat.lehmann2018-PET"]]
  fit2 <- fits[["dat.lehmann2018-PET_neg"]]

  expect_equal(bridgesampling::bf(bridge_sampler(fit2), bridge_sampler(fit1))$bf, 1, tolerance = 0.01)

  fit1 <- fits[["dat.lehmann2018-3PSM"]]
  fit2 <- fits[["dat.lehmann2018-3PSM_neg"]]

  expect_equal(bridgesampling::bf(bridge_sampler(fit2), bridge_sampler(fit1))$bf, 1, tolerance = 0.01)
})

test_that("logml returns scalar log marginal likelihood, can be applied to both bridge and brma", {

  skip_if_missing_fits("bcg_meta-analysis")

  name     <- "bcg_meta-analysis"
  fit_brma <- fits[[name]]

  logml_result1 <- bridgesampling::logml(fit_brma)
  logml_result2 <- RoBMA::logml(fit_brma)

  # consistency check
  expect_equal(logml_result1, logml_result2, tolerance = 0.01)
})


test_that("bf computes Bayes factor between two models, can be applied to both bridge and brma", {

  skip_if_missing_fits(c("bcg_meta-analysis", "bcg_meta-regression"))

  fit_brma1 <- fits[["bcg_meta-analysis"]]
  fit_brma2 <- fits[["bcg_meta-regression"]]

  bf_result1 <- RoBMA::bf(fit_brma1, fit_brma2)
  bf_result2 <- bridgesampling::bf(fit_brma1, fit_brma2)
  bf_result3 <- RoBMA::bayes_factor(fit_brma1, fit_brma2)
  bf_result4 <- bridgesampling::bayes_factor(fit_brma1, fit_brma2)

  # consistency check
  expect_equal(bf_result1, bf_result2, tolerance = 0.01)
  expect_equal(bf_result1, bf_result3, tolerance = 0.01)
  expect_equal(bf_result1, bf_result4, tolerance = 0.01)
})


test_that("post_prob computes posterior model probabilities", {

  skip_if_missing_fits(c("bcg_meta-analysis", "bcg_meta-regression"))

  fit_brma1 <- fits[["bcg_meta-analysis"]]
  fit_brma2 <- fits[["bcg_meta-regression"]]

  pp_result1 <- RoBMA::post_prob(fit_brma1, fit_brma2)
  pp_result2 <- bridgesampling::post_prob(fit_brma1, fit_brma2)

  # consistency check
  expect_equal(pp_result1, pp_result2, tolerance = 0.01)
})


test_that("post_prob respects prior model probabilities", {

  skip_if_missing_fits(c("bcg_meta-analysis", "bcg_meta-regression"))

  fit_brma  <- fits[["bcg_meta-analysis"]]
  fit_brma2 <- fits[["bcg_meta-regression"]]

  # equal priors (default)
  pp_equal <- RoBMA::post_prob(fit_brma, fit_brma2)

  # unequal priors
  pp_unequal <- RoBMA::post_prob(fit_brma, fit_brma2, prior_prob = c(0.9, 0.1))

  # with stronger prior on first model, its posterior should increase
  expect_lt(unname(pp_equal[1]), unname(pp_unequal[1]))
  expect_equal(sum(pp_unequal), 1)
})


# ---------------------------------------------------------------------------- #
# error handling tests
# ---------------------------------------------------------------------------- #

test_that("bf rejects non-brma object", {

  skip_if_missing_fits("bcg_meta-analysis")

  fit_brma <- fits[["bcg_meta-analysis"]]
  expect_error(
    bf(fit_brma, "not a brma"),
    "needs to be of class 'brma'"
  )
})


test_that("post_prob rejects single model", {

  skip_if_missing_fits("bcg_meta-analysis")

  fit_brma <- fits[["bcg_meta-analysis"]]
  expect_error(
    post_prob(fit_brma),
    "Only one object"
  )
})


test_that("bridge_sampler rejects fits without marginal likelihood", {

  model_names <- c("bcg_meta-analysis", "brma.mv_block_mvn_fixed_random_null")
  skip_if_missing_fits(model_names)

  for (name in model_names) {
    fit_brma <- load_fit(name, validate = FALSE)
    fit_brma[["marglik"]] <- NULL
    expect_error(
      bridge_sampler(fit_brma),
      "add_marglik",
      info = name
    )
  }
})

test_that("add_marglik rejects product-space model-averaging objects", {

  product_names <- c(
    "dat.lehmann2018_BMA.norm",
    "bcg_BMA.glmm",
    "dat.lehmann2018_RoBMA"
  )
  skip_if_missing_fits(product_names)

  for (name in product_names) {
    expect_error(
      add_marglik(load_fit(name, validate = FALSE)),
      "Marginal likelihood is not available for product-space"
    )
    expect_error(
      bridge_sampler(load_fit(name, validate = FALSE)),
      "Marginal likelihood is not available for product-space"
    )
    expect_error(
      logml(load_fit(name, validate = FALSE)),
      "Marginal likelihood is not available for product-space"
    )
  }
})
