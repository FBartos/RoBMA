context("Model fitting for bselmodel")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
skip_refit_if_cached("bselmodel")

### Uses examples from the metafor package
test_that("Test against metafor::selmodel", {

  skip_on_cran()
  skip_if_not_installed("metadat")
  skip_if_not_installed("metafor")

  ### fit simple meta-analytic model to difference in two proportions
  data(dat.lehmann2018, package = "metadat")
  fit_rma.metafor      <- metafor::rma(yi, vi, data = dat.lehmann2018, method = "ML")
  fit_selmodel.metafor <- metafor::selmodel(fit_rma.metafor, type = "stepfun", alternative = "greater", steps = .025)

  # using RoBMA package
  fit.bselmodel <- bselmodel(yi, vi, data = dat.lehmann2018, measure = "SMD", seed = 1, silent = TRUE)
  save_fit("dat.lehmann2018-3PSM", fit.bselmodel, info = list(metafor = fit_selmodel.metafor))

  expect_equal(fit_selmodel.metafor$beta[[1]],  fit.bselmodel$summary["mu","Mean"],  tolerance = 0.01)
  expect_equal(sqrt(fit_selmodel.metafor$tau2), fit.bselmodel$summary["tau","Mean"], tolerance = 0.01)
  expect_equal(fit_selmodel.metafor$delta[2],   fit.bselmodel$summary["omega[0.025,1]","Mean"], tolerance = 0.10)
})


