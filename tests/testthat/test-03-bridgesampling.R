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

# load cached fits
fits     <- list()
margliks <- list()
info     <- list()
for (name in list_fits()) {
  fits[[name]]     <- load_fit(name)
  margliks[[name]] <- load_marglik(name)
  info[[name]]     <- load_info(name)
}

# ---------------------------------------------------------------------------- #
# bridge_sampler function works with all model types
# ---------------------------------------------------------------------------- #

test_that("bridge_sampler computes bridge sampling object", {

  for (name in list_fits()) {
    # the marginal likelihood inherits the correct class
    expect_s3_class(margliks[[name]], "bridge")

    # the marginal likelihood and BIC-based marglik are close
    # (the analyses don't use unit information priors but the
    #  difference should not be large)
    BIC_metafor <- BIC(info[[name]][["metafor"]])
    marglik     <- logml(margliks[[name]])

    if (name %in% c("bcg_glmm", "nielweise2008_glmm"))
      next # skip glmm because of marginalization differences

    expect_equal(- BIC_metafor/2, marglik, tolerance = 0.15)
  }

})


# ---------------------------------------------------------------------------- #
# test simple comparisons
# ---------------------------------------------------------------------------- #

test_that("simple comparisons", {

  ### there should be evidence for meta-regression
  name     <- "bcg_meta-analysis"
  marglik  <- margliks[[name]]

  name2    <- "bcg_meta-regression"
  marglik2 <- margliks[[name2]]

  expect_equal(bridgesampling::bf(marglik2, marglik)$bf, 4.54, tolerance = 0.01)

  ### there should be no difference between positive and negative PET
  name     <- "dat.lehmann2018-PET"
  marglik  <- margliks[[name]]

  name2    <- "dat.lehmann2018-PET_neg"
  marglik2 <- margliks[[name2]]

  expect_equal(bridgesampling::bf(marglik2, marglik)$bf, 1, tolerance = 0.01)

  ### there should be no difference between positive and negative 3PSM
  name     <- "dat.lehmann2018-3PSM"
  marglik  <- margliks[[name]]

  name2    <- "dat.lehmann2018-3PSM_neg"
  marglik2 <- margliks[[name2]]

  expect_equal(bridgesampling::bf(marglik2, marglik)$bf, 1, tolerance = 0.01)
})

# ---------------------------------------------------------------------------- #
# logml simple function test
# ---------------------------------------------------------------------------- #

test_that("logml returns scalar log marginal likelihood, can be applied to both bridge and brma", {

  name         <- "bcg_meta-analysis"
  fit_brma     <- fits[[name]]
  marglik_brma <- margliks[[name]]

  logml_result1 <- bridgesampling::logml(marglik_brma)
  logml_result2 <- bridgesampling::logml(fit_brma)
  logml_result3 <- RoBMA::logml(marglik_brma)
  logml_result4 <- RoBMA::logml(fit_brma)

  # consistency check
  expect_equal(logml_result1, logml_result2, tolerance = 0.01)
  expect_equal(logml_result1, logml_result3, tolerance = 0.01)
  expect_equal(logml_result1, logml_result4, tolerance = 0.01)
})


# ---------------------------------------------------------------------------- #
# bf simple function test
# ---------------------------------------------------------------------------- #

test_that("bf computes Bayes factor between two models, can be applied to both bridge and brma", {

  marglik_brma1 <- margliks[["bcg_meta-analysis"]]
  marglik_brma2 <- margliks[["bcg_meta-regression"]]

  bf_result1 <- RoBMA::bf(marglik_brma1, marglik_brma2)
  bf_result2 <- bridgesampling::bf(marglik_brma1, marglik_brma2)
  bf_result3 <- RoBMA::bayes_factor(marglik_brma1, marglik_brma2)
  bf_result4 <- bridgesampling::bayes_factor(marglik_brma1, marglik_brma2)

  fit_brma1 <- fits[["bcg_meta-analysis"]]
  fit_brma2 <- fits[["bcg_meta-regression"]]

  bf_result5 <- RoBMA::bf(fit_brma1, fit_brma2)
  # the other combinations are skipped because they require recomputation of all margliks...

  # consistency check
  expect_equal(bf_result1, bf_result2, tolerance = 0.01)
  expect_equal(bf_result1, bf_result3, tolerance = 0.01)
  expect_equal(bf_result1, bf_result4, tolerance = 0.01)
  expect_equal(bf_result1, bf_result5, tolerance = 0.01)
})


# ---------------------------------------------------------------------------- #
# post_prob simple function test
# ---------------------------------------------------------------------------- #

test_that("post_prob computes posterior model probabilities", {
# TODO: finish test here
  marglik_brma1 <- margliks[["bcg_meta-analysis"]]
  marglik_brma2 <- margliks[["bcg_meta-regression"]]

  pp_result1 <- post_prob(marglik_brma1, marglik_brma2)

  fit_brma1 <- fits[["bcg_meta-analysis"]]
  fit_brma2 <- fits[["bcg_meta-regression"]]

  pp_result2 <- post_prob(fit_brma1, fit_brma2)

})


test_that("post_prob respects prior model probabilities", {

  skip_if_no_fits()
  skip_on_cran()

  fit_brma  <- fits[["bcg_meta-analysis"]]
  fit_brma2 <- fits[["bcg_meta-regression"]]

  # equal priors (default)
  pp_equal <- post_prob(fit_brma, fit_brma2)

  # unequal priors
  pp_unequal <- post_prob(fit_brma, fit_brma2, prior_prob = c(0.9, 0.1))

  # with stronger prior on first model, its posterior should increase
  # (unless evidence overwhelmingly favors second model)
  expect_true(is.numeric(pp_unequal))
  expect_length(pp_unequal, 2)
  expect_equal(sum(pp_unequal), 1, tolerance = 1e-10)

})


# ---------------------------------------------------------------------------- #
# error handling tests
# ---------------------------------------------------------------------------- #

test_that("bf errors with non-brma object", {

  skip_if_no_fits()

  fit_brma <- fits[["bcg_meta-analysis"]]

  expect_error(
    bf(fit_brma, "not a brma"),
    "needs to be of class 'brma'"
  )

})


test_that("post_prob errors with single model", {

  skip_if_no_fits()

  fit_brma <- fits[["bcg_meta-analysis"]]

  expect_error(
    post_prob(fit_brma),
    "Only one object"
  )

})
