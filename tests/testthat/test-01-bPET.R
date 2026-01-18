context("Model fitting for bPET")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
skip_on_cran()
skip_if_not_installed("metadat")
skip_if_not_installed("metafor")
skip_refit_if_cached("bPET")

### Uses examples from the metafor package
test_that("Test against metafor::rma.uni with mods = ~ sei ", {
  ### fit simple meta-analytic model to difference in two proportions
  data(dat.lehmann2018, package = "metadat")
  fit_PET.metafor <- metafor::rma(yi, vi, mods = ~ sqrt(vi), data = dat.lehmann2018)

  # using RoBMA package
  fit.bPET <- bPET(yi, vi, data = dat.lehmann2018, measure = "SMD", seed = 1, silent = TRUE)
  fit.bPET <- add_marglik(fit.bPET)
  fit.bPET <- suppressWarnings(add_loo(fit.bPET))
  save_fit("dat.lehmann2018-PET", fit.bPET, info = list(metafor = fit_PET.metafor))

  expect_equal(fit_PET.metafor$beta[[1]], fit.bPET$summary["mu", "Mean"], tolerance = 0.05)
  expect_equal(sqrt(fit_PET.metafor$tau2), fit.bPET$summary["tau", "Mean"], tolerance = 0.05)
  expect_equal(fit_PET.metafor$beta[[2]], fit.bPET$summary["PET", "Mean"], tolerance = 0.20) # the PET prior regularizes the estimate a bit
})


test_that("Test against metafor::rma.uni with mods = ~ sei (with negative effect sizes)", {
  ### fit simple meta-analytic model to difference in two proportions
  data(dat.lehmann2018, package = "metadat")
  dat.lehmann2018$yi <- -dat.lehmann2018$yi
  fit_PET.metafor <- metafor::rma(yi, vi, mods = ~ sqrt(vi), data = dat.lehmann2018)

  # using RoBMA package
  fit.bPET <- bPET(yi, vi, data = dat.lehmann2018, measure = "SMD", seed = 1, silent = TRUE)
  fit.bPET <- add_marglik(fit.bPET)
  fit.bPET <- suppressWarnings(add_loo(fit.bPET))
  save_fit("dat.lehmann2018-PET_neg", fit.bPET, info = list(metafor = fit_PET.metafor))

  expect_equal(fit_PET.metafor$beta[[1]], fit.bPET$summary["mu", "Mean"], tolerance = 0.05)
  expect_equal(sqrt(fit_PET.metafor$tau2), fit.bPET$summary["tau", "Mean"], tolerance = 0.05)
  expect_equal(fit_PET.metafor$beta[[2]], -fit.bPET$summary["PET", "Mean"], tolerance = 0.20) # the PET prior regularizes the estimate a bit (and bPET keeps it positive)
})
