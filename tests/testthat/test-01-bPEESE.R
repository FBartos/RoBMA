context("Model fitting for bPEESE")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
skip_on_cran()
skip_if_not_installed("metadat")
skip_if_not_installed("metafor")
skip_refit_if_cached("bPEESE")

### Uses examples from the metafor package
test_that("Test against metafor::rma.uni with mods = ~ sei^2 ", {

  ### fit simple meta-analytic model to difference in two proportions
  data(dat.lehmann2018, package = "metadat")
  fit_PEESE.metafor <- metafor::rma(yi, vi, mods = ~ vi, data = dat.lehmann2018)

  # using RoBMA package
  fit.bPEESE  <- bPEESE(yi, vi, data = dat.lehmann2018, measure = "SMD", seed = 1, silent = TRUE)
  fit.marglik <- bridge_sampler(fit.bPEESE)
  save_fit("dat.lehmann2018-PEESE", fit.bPEESE, fit.marglik, info = list(metafor = fit_PEESE.metafor))

  expect_equal(fit_PEESE.metafor$beta[[1]],   fit.bPEESE$summary["mu","Mean"],    tolerance = 0.05)
  expect_equal(sqrt(fit_PEESE.metafor$tau2),  fit.bPEESE$summary["tau","Mean"],   tolerance = 0.05)
  expect_equal(fit_PEESE.metafor$beta[[2]],   fit.bPEESE$summary["PEESE","Mean"], tolerance = 0.05)
})


