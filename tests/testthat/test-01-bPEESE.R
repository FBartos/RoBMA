context("Model fitting for bPEESE")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
skip_on_cran()
skip_if_not_installed("metadat")
skip_if_not_installed("metafor")
skip_refit_if_cached("bPEESE")

### Uses examples from the metafor package
test_that("bPEESE fits PEESE metafor-reference model", {
  ### fit PEESE model
  data(dat.lehmann2018, package = "metadat")
  fit_PEESE.metafor <- metafor::rma(yi, vi, mods = ~vi, data = dat.lehmann2018)

  # using RoBMA package
  fit.bPEESE <- bPEESE(yi, vi, data = dat.lehmann2018, measure = "SMD", seed = 1, silent = TRUE)
  fit.bPEESE <- add_marglik(fit.bPEESE)
  fit.bPEESE <- suppressWarnings(add_loo(fit.bPEESE))
  save_fit("dat.lehmann2018-PEESE", fit.bPEESE, info = list(metafor = fit_PEESE.metafor))

  expect_s3_class(fit.bPEESE, "bPEESE")
})
