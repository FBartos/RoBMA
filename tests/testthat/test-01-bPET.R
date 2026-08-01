context("Model fitting for bPET")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
skip_on_cran()
skip_if_not_installed("metadat")
skip_if_not_installed("metafor")
skip_refit_if_cached("bPET")

### Uses examples from the metafor package
test_that("bPET fits PET metafor-reference model", {
  ### fit PET model
  data(dat.lehmann2018, package = "metadat")
  fit_PET.metafor <- metafor::rma(yi, vi, mods = ~ sqrt(vi), data = dat.lehmann2018)

  # using RoBMA package
  fit.bPET <- bPET(yi, vi, data = dat.lehmann2018, measure = "SMD", seed = 1, silent = TRUE)
  fit.bPET <- add_marglik(fit.bPET)
  fit.bPET <- suppressWarnings(add_loo(fit.bPET))
  save_fit("dat.lehmann2018-PET", fit.bPET, info = list(metafor = fit_PET.metafor))

  expect_s3_class(fit.bPET, "bPET")
})


test_that("bPET fits negative-direction PET metafor-reference model", {

  skip_if_fit_not_active("dat.lehmann2018-PET_neg")
  ### fit PET model
  data(dat.lehmann2018, package = "metadat")
  dat.lehmann2018$yi <- -dat.lehmann2018$yi
  fit_PET.metafor <- metafor::rma(yi, vi, mods = ~ sqrt(vi), data = dat.lehmann2018)

  # using RoBMA package
  fit.bPET <- bPET(yi, vi, data = dat.lehmann2018, measure = "SMD", seed = 1, silent = TRUE)
  fit.bPET <- add_marglik(fit.bPET)
  fit.bPET <- suppressWarnings(add_loo(fit.bPET))
  save_fit("dat.lehmann2018-PET_neg", fit.bPET, info = list(metafor = fit_PET.metafor))

  expect_s3_class(fit.bPET, "bPET")
})

test_that("bPET fits PET meta-regression model", {

  skip_if_fit_not_active("dat.lehmann2018-PETreg")
  ### fit PET model
  data(dat.lehmann2018, package = "metadat")
  fit_PET.metafor <- metafor::rma(yi, vi, mods = ~ sqrt(vi) + Preregistered, data = dat.lehmann2018)

  # using RoBMA package
  fit.bPET <- bPET(yi, vi, mods = ~ Preregistered, data = dat.lehmann2018, measure = "SMD", seed = 1, silent = TRUE)
  fit.bPET <- add_marglik(fit.bPET)
  fit.bPET <- suppressWarnings(add_loo(fit.bPET))
  save_fit("dat.lehmann2018-PETreg", fit.bPET, info = list(metafor = fit_PET.metafor, mods = "Preregistered"))

  expect_s3_class(fit.bPET, "bPET")
})
