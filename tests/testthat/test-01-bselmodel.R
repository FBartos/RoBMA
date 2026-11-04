context("Model fitting for bselmodel")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
skip_on_cran()
skip_if_not_installed("metadat")
skip_if_not_installed("metafor")
skip_refit_if_cached("bselmodel")

### Uses examples from the metafor package
test_that("bselmodel fits one-step metafor-reference selection model", {
  ### fit selection model
  data(dat.lehmann2018, package = "metadat")
  fit_rma.metafor <- metafor::rma(yi, vi, data = dat.lehmann2018, method = "ML")
  fit_selmodel.metafor <- metafor::selmodel(fit_rma.metafor, type = "stepfun", alternative = "greater", steps = .025)

  # using RoBMA package
  fit.bselmodel <- bselmodel(yi, vi, data = dat.lehmann2018, measure = "SMD", seed = 1, silent = TRUE)
  fit.bselmodel <- add_marglik(fit.bselmodel)
  fit.bselmodel <- suppressWarnings(add_loo(fit.bselmodel))
  save_fit("dat.lehmann2018-3PSM", fit.bselmodel, info = list(metafor = fit_selmodel.metafor))

  expect_s3_class(fit.bselmodel, "bselmodel")
})

test_that("bselmodel fits two-step metafor-reference selection model", {

  skip_if_fit_not_active("dat.lehmann2018-4PSM")
  ### fit selection model
  data(dat.lehmann2018, package = "metadat")
  fit_rma.metafor <- metafor::rma(yi, vi, data = dat.lehmann2018, method = "ML")
  fit_selmodel.metafor <- metafor::selmodel(fit_rma.metafor, type = "stepfun", alternative = "greater", steps = c(.025, 0.50))

  # using RoBMA package
  fit.bselmodel <- bselmodel(yi, vi, data = dat.lehmann2018, measure = "SMD", seed = 1, silent = TRUE, steps = c(0.025, 0.50))
  fit.bselmodel <- add_marglik(fit.bselmodel)
  fit.bselmodel <- suppressWarnings(add_loo(fit.bselmodel))
  save_fit("dat.lehmann2018-4PSM", fit.bselmodel, info = list(metafor = fit_selmodel.metafor))

  expect_s3_class(fit.bselmodel, "bselmodel")
})

test_that("bselmodel fits negative-direction metafor-reference selection model", {

  skip_if_fit_not_active("dat.lehmann2018-3PSM_neg")
  ### fit selection model
  data(dat.lehmann2018, package = "metadat")
  dat.lehmann2018$yi <- -dat.lehmann2018$yi
  fit_rma.metafor <- metafor::rma(yi, vi, data = dat.lehmann2018, method = "ML")
  fit_selmodel.metafor <- metafor::selmodel(fit_rma.metafor, type = "stepfun", alternative = "less", steps = .025)

  # using RoBMA package
  fit.bselmodel <- bselmodel(yi, vi, data = dat.lehmann2018, measure = "SMD", seed = 1, silent = TRUE)
  fit.bselmodel <- add_marglik(fit.bselmodel)
  fit.bselmodel <- suppressWarnings(add_loo(fit.bselmodel))
  save_fit("dat.lehmann2018-3PSM_neg", fit.bselmodel, info = list(metafor = fit_selmodel.metafor))

  expect_s3_class(fit.bselmodel, "bselmodel")
})

test_that("bselmodel fits selection meta-regression model", {

  skip_if_fit_not_active("dat.lehmann2018-3PSMreg")
  ### fit selection model
  data(dat.lehmann2018, package = "metadat")
  fit_rma.metafor <- metafor::rma(yi, vi, mods = ~ Preregistered, data = dat.lehmann2018, method = "ML")
  fit_selmodel.metafor <- metafor::selmodel(fit_rma.metafor, type = "stepfun", steps = .025)

  # using RoBMA package
  fit.bselmodel <- bselmodel(yi, vi, mods = ~ Preregistered, data = dat.lehmann2018, measure = "SMD", seed = 1, silent = TRUE, steps = .025)
  fit.bselmodel <- add_marglik(fit.bselmodel)
  fit.bselmodel <- suppressWarnings(add_loo(fit.bselmodel))
  save_fit("dat.lehmann2018-3PSMreg", fit.bselmodel, info = list(metafor = fit_selmodel.metafor, mods = "Preregistered"))

  expect_s3_class(fit.bselmodel, "bselmodel")
})
