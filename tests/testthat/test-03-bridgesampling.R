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

test_that("bridge sampling marginal likelihood is close to BIC for metafor-reference fits", {

  metafor_names <- list_fits(has_marglik = TRUE, has_metafor = TRUE)

  for (name in metafor_names) {

    # the marginal likelihood and BIC-based marglik are close
    # (the analyses don't use unit information priors but the
    #  difference should not be large)
    BIC_metafor <- BIC(info[[name]][["metafor"]])
    marglik     <- logml(fits[[name]])

    if (name %in% c("bcg_glmm", "nielweise2008_glmm", "bcg_glmm_reg")) {
      next
    } # skip glmm because of marginalization differences

    if (name %in% c("bcg_meta-regression", "bcg_meta-regression2", "bcg_meta-regression3", "bcg_meta-regression3b", "bcg_meta-regression4", "bcg_meta-regression4b")) {
      next
    } # skip because of scaled priors differences

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
  expect_true(pp_equal[1] <  pp_unequal[1])
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

  skip_if_missing_fits("bcg_meta-analysis")

  # create a fit without marglik
  fit_brma <- fits[["bcg_meta-analysis"]]
  fit_brma[["marglik"]] <- NULL
  expect_error(
    bridge_sampler(fit_brma),
    "add_marglik"
  )
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
  }
})
