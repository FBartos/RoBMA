context("Nested models for IWMDE numerical oracles")

source(testthat::test_path("common-functions.R"))

skip_on_cran()
skip_if_not_installed("metadat")
skip_refit_if_cached("iwmde-oracle-nested")


test_that("nested models support IWMDE bridge oracles", {

  skip_if_not_certification("These nested bridge fits certify numerical accuracy.")
  spike <- BayesTools::prior(
    "spike",
    parameters = list(location = 0)
  )

  data(dat.nielweise2008, package = "metadat")
  fit_glmm_null <- brma.glmm(
    x1i          = x1i,
    t1i          = t1i,
    x2i          = x2i,
    t2i          = t2i,
    data         = dat.nielweise2008,
    measure      = "IRR",
    prior_effect = spike,
    seed         = 1,
    silent       = TRUE
  )
  fit_glmm_null <- add_marglik(fit_glmm_null)
  save_fit("nielweise2008_glmm_effect_null", fit_glmm_null)

  data(dat.lehmann2018, package = "metadat")
  fit_selection_null <- bselmodel(
    yi           = yi,
    vi           = vi,
    data         = dat.lehmann2018,
    measure      = "SMD",
    prior_effect = spike,
    seed         = 1,
    silent       = TRUE
  )
  fit_selection_null <- add_marglik(fit_selection_null)
  save_fit("dat.lehmann2018-3PSM_effect_null", fit_selection_null)

  skip_if_missing_fits(c(
    "brma.mv_block_mvn",
    "brma.mv_block_mvn_fixed_random_null"
  ), active_only = FALSE)
  fit_known_v_full <- add_marglik(load_fit("brma.mv_block_mvn"))
  fit_known_v_null <- add_marglik(load_fit(
    "brma.mv_block_mvn_fixed_random_null"
  ))
  save_fit("iwmde_known_v_tau_full", fit_known_v_full)
  save_fit("iwmde_known_v_tau_null", fit_known_v_null)

  expect_s3_class(fit_glmm_null, "brma.glmm")
  expect_s3_class(fit_selection_null, "bselmodel")
  expect_true(is.finite(logml(fit_glmm_null)))
  expect_true(is.finite(logml(fit_selection_null)))
  expect_true(is.finite(logml(fit_known_v_full)))
  expect_true(is.finite(logml(fit_known_v_null)))
})
