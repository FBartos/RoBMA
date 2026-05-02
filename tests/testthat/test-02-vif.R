context("VIF")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-contracts.R"))
source(testthat::test_path("helper-test-matrix.R"))
source(testthat::test_path("helper-metafor.R"))

skip_if_no_fits()
skip_if_not_installed("metafor")

fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)

for_each_case(vif_cases(), function(case) {
  test_that_case("VIF matches metafor", case, {
    expect_vif_matches_metafor(case)
  })
})

test_that("VIF returns optional posterior correlation diagnostics", {

  name <- "bcg_meta-regression"
  skip_if_missing_fits(name)

  fit_brma      <- fits[[name]]
  vif_with_post <- vif(fit_brma)
  post_cor      <- vif_with_post[["posterior_correlation"]]

  expect_s3_class(vif_with_post, "vif.brma")
  expect_vif_table(vif_with_post[["vif"]], 2)
  expect_type(post_cor, "double")
  expect_equal(nrow(post_cor), ncol(post_cor))
  expect_equal(unname(diag(post_cor)), rep(1, nrow(post_cor)), tolerance = 1e-12)
  expect_equal(post_cor, t(post_cor), tolerance = 1e-12)
  expect_true(all(is.finite(post_cor)))

  vif_without_post <- vif(fit_brma, posterior_correlation = FALSE)
  expect_null(vif_without_post[["posterior_correlation"]])
})

test_that("VIF rejects models without moderators", {

  model_names <- c(
    "bcg_meta-analysis",
    "konstantopoulos2011_3lvl",
    "nielweise2008_glmm",
    "dat.lehmann2018-PET",
    "dat.lehmann2018-3PSM"
  )
  skip_if_missing_fits(model_names)

  for (name in model_names) {
    expect_error(
      vif(fits[[name]]),
      "only meaningful for models with moderators",
      info = paste("vif rejects model without moderators for", name)
    )
  }
})
