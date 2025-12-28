context("Model fitting for brma.uni")

### Uses examples from the metafor package
# - compares the results


test_that("Test against metafor::rma.uni", {

  skip_on_cran()
  skip_if_not_installed("metadat")
  skip_if_not_installed("metafor")

  ### fit simple meta-analytic model
  data(dat.bcg, package = "metadat")
  dat                 <- metafor::escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
  fit_simple.metafor  <- metafor::rma(yi = yi, vi = vi, data = dat, method = "REML")

  # using RoBMA package
  fit_simple.brma     <- brma.uni(yi = yi, vi = vi, data = dat, measure = "RR", seed = 1)

  expect_equal(fit_simple.metafor$beta[[1]],  fit_simple.brma$summary["mu","Mean"],  tolerance = 0.05)
  expect_equal(sqrt(fit_simple.metafor$tau2), fit_simple.brma$summary["tau","Mean"], tolerance = 0.05)


  ### fit meta-regression
  fit_mods.metafor <- metafor::rma(yi, vi, mods = ~ ablat + year, data = dat)

  # using RoBMA package
  fit_mods.brma    <- brma.uni(yi, vi, mods = ~ ablat + year, data = dat, measure = "RR", seed = 1)

})


test_that("Test against metafor::rma.ls", {

  skip_on_cran()
  skip_if_not_installed("metadat")
  skip_if_not_installed("metafor")

  ### fit location-scale meta-analytic model
  data(dat.bangertdrowns2004, package = "metadat")
  dat       <- dat.bangertdrowns2004
  dat$ni100 <- dat$ni/100
  fit_scale.metafor <- metafor::rma(yi, vi, mods = ~ ni100 + meta, scale = ~ ni100 + imag, data = dat)

  # using RoBMA package
  fit_scale.brma     <- brma.uni(yi, vi, mods = ~ ni100 + meta, scale = ~ ni100 + imag, data = dat, measure = "SMD", seed = 1)

  expect_equal(fit_simple.metafor$beta[[1]],  fit_simple.brma$summary["mu","Mean"],  tolerance = 0.05)
  expect_equal(sqrt(fit_simple.metafor$tau2), fit_simple.brma$summary["tau","Mean"], tolerance = 0.05)


  ### fit meta-regression
  fit_mods.metafor <- metafor::rma(yi, vi, mods = ~ ablat + year, data = dat)

  # using RoBMA package
  fit_mods.brma    <- brma.uni(yi, vi, mods = ~ ablat + year, data = dat, measure = "RR", seed = 1)

})
