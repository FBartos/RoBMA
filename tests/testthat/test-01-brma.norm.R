context("Model fitting for brma.norm")


# Load common test helpers
source(testthat::test_path("common-functions.R"))

### Uses examples from the metafor package
test_that("Test against metafor::rma.uni", {

  skip_on_cran()
  skip_if_not_installed("metadat")
  skip_if_not_installed("metafor")

  ### fit simple meta-analytic model
  data(dat.bcg, package = "metadat")
  dat                 <- metafor::escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
  fit_simple.metafor  <- metafor::rma(yi = yi, vi = vi, data = dat, method = "REML")

  # using RoBMA package
  fit_simple.brma <- brma(yi = yi, vi = vi, data = dat, measure = "RR", seed = 1)
  save_fit(fit_simple.brma, "bcg_meta-analysis")
  expect_equal(fit_simple.metafor$beta[[1]],  fit_simple.brma$summary["mu","Mean"],  tolerance = 0.05)
  expect_equal(sqrt(fit_simple.metafor$tau2), fit_simple.brma$summary["tau","Mean"], tolerance = 0.05)


  ### fit meta-regression
  fit_mods.metafor <- metafor::rma(yi, vi, mods = ~ ablat + year, data = dat)

  # using RoBMA package
  fit_mods.brma <- brma(yi, vi, mods = ~ ablat + year, data = dat, measure = "RR", seed = 1)
  save_fit(fit_mods.brma, "bcg_meta-regression")

  expect_equal(fit_mods.metafor$beta[[1]],  fit_mods.brma$summary["(mu) intercept","Mean"],    tolerance = 1.5)
  expect_equal(fit_mods.metafor$beta[[2]],  fit_mods.brma$summary["(mu) ablat","Mean"],        tolerance = 0.01)
  expect_equal(fit_mods.metafor$beta[[3]],  fit_mods.brma$summary["(mu) year","Mean"],         tolerance = 0.001)
  # the large difference in the intercept is a consequence in the large impact of small changes in the year coefficient
  # (they get multiplied by 2000)

  # verify that using scaled predictors removes the issue
  fit_mods.metafor.scaled <- metafor::rma(yi, vi, mods = ~ scale(ablat) + scale(year), data = dat)
  fit_mods.brma.scaled    <- BayesTools::JAGS_estimates_table(fit_mods.brma$fit)
  expect_equal(fit_mods.metafor.scaled$beta[[1]],  fit_mods.brma.scaled["(mu) intercept","Mean"],    tolerance = 0.05)
  expect_equal(fit_mods.metafor.scaled$beta[[2]],  fit_mods.brma.scaled["(mu) ablat","Mean"],        tolerance = 0.05)
  expect_equal(fit_mods.metafor.scaled$beta[[3]],  fit_mods.brma.scaled["(mu) year","Mean"],         tolerance = 0.05)

})

test_that("Test against metafor::rma.ls", {

  skip_on_cran()
  skip_if_not_installed("metadat")
  skip_if_not_installed("metafor")

  ### fit location-scale model
  data(dat.bangertdrowns2004, package = "metadat")
  dat       <- dat.bangertdrowns2004
  dat$ni100 <- dat$ni/100
  dat$meta  <- as.factor(dat$meta)
  dat$imag  <- as.factor(dat$imag)
  fit_scale.metafor <- suppressWarnings(metafor::rma(yi, vi, mods = ~ meta + ni100, scale = ~ ni100, data = dat, method = "REML"))

  # using RoBMA package (using wide priors for scale to remove shrinkage -- the likelihood is very wide)
  fit_scale.brma <- suppressWarnings(brma(yi, vi, mods = ~ meta + ni100, scale = ~ ni100, data = dat, measure = "SMD", seed = 1))
  save_fit(fit_scale.brma, "bangertdrowns2004_location-scale")

  expect_equal(fit_scale.metafor$beta[[1]],  fit_scale.brma$summary["(mu) intercept","Mean"],  tolerance = 0.05)
  expect_equal(fit_scale.metafor$beta[[2]],  fit_scale.brma$summary["(mu) meta[1]","Mean"],    tolerance = 0.05)
  expect_equal(fit_scale.metafor$beta[[3]],  fit_scale.brma$summary["(mu) ni100","Mean"],      tolerance = 0.05)

  expect_equal(0.5 * fit_scale.metafor$alpha[[1]], log(fit_scale.brma$summary["(log_tau) exp(intercept)","Mean"]), tolerance = 0.15)
  expect_equal(0.5 * fit_scale.metafor$alpha[[2]], fit_scale.brma$summary["(log_tau) ni100","Mean"],               tolerance = 0.15)

})

test_that("Test against metafor::rma.mv (3lvl)", {

  skip_on_cran()
  skip_if_not_installed("metadat")
  skip_if_not_installed("metafor")

  ### fit 3lvl model
  data(dat.konstantopoulos2011, package = "metadat")

  fit_3lvl.metafor <- metafor::rma.mv(yi, vi, random = ~ school | district, data = dat.konstantopoulos2011)

  # using RoBMA package
  fit_3lvl.brma <- brma(yi, vi, study_ids = district, data = dat.konstantopoulos2011, measure = "SMD", seed = 1)
  save_fit(fit_3lvl.brma, "konstantopoulos2011_3lvl")

  expect_equal(fit_3lvl.metafor$beta[[1]],  fit_3lvl.brma$summary["mu","Mean"],   tolerance = 0.01)
  expect_equal(sqrt(fit_3lvl.metafor$tau2), fit_3lvl.brma$summary["tau","Mean"],  tolerance = 0.01)
  expect_equal(fit_3lvl.metafor$rho,        fit_3lvl.brma$summary["rho","Mean"],  tolerance = 0.05)
})
