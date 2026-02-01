context("Model fitting for bselmodel")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
skip_on_cran()
skip_if_not_installed("metadat")
skip_refit_if_cached("RoBMA")
skip("TODO")

### Uses examples from the metafor package
test_that("Test RoBMA", {

  data(dat.lehmann2018, package = "metadat")
  fit.RoBMA <- RoBMA(yi, vi, data = dat.lehmann2018, measure = "SMD", seed = 1, silent = TRUE)
  fit.RoBMA <- suppressWarnings(add_loo(fit.RoBMA))
  save_fit("dat.lehmann2018-RoBMA", fit.RoBMA, info = list())

})

test_that("Test RoBMA (with negative effect sizes)", {

  data(dat.lehmann2018, package = "metadat")
  dat.lehmann2018$yi <- -dat.lehmann2018$yi
  fit.RoBMA <- RoBMA(yi, vi, data = dat.lehmann2018, measure = "SMD", seed = 1, silent = TRUE)
  fit.RoBMA <- suppressWarnings(add_loo(fit.RoBMA))

  save_fit("dat.lehmann2018-RoBMA-neg", fit.RoBMA, info = list())
})
