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
  fit_simple.brma <- brma(yi = yi, vi = vi, data = dat, measure = "RR", seed = 1)

  expect_equal(fit_simple.metafor$beta[[1]],  fit_simple.brma$summary["mu","Mean"],  tolerance = 0.05)
  expect_equal(sqrt(fit_simple.metafor$tau2), fit_simple.brma$summary["tau","Mean"], tolerance = 0.05)


  ### fit meta-regression
  fit_mods.metafor <- metafor::rma(yi, vi, mods = ~ ablat + year, data = dat)

  # using RoBMA package
  fit_mods.brma <- brma(yi, vi, mods = ~ ablat + year, data = dat, measure = "RR", seed = 1)

  expect_equal(fit_mods.metafor$beta[[1]],  fit_mods.brma$summary["(mu) intercept","Mean"],    tolerance = 1.5)
  expect_equal(fit_mods.metafor$beta[[2]],  fit_mods.brma$summary["(mu) ablat","Mean"],        tolerance = 0.01)
  expect_equal(fit_mods.metafor$beta[[3]],  fit_mods.brma$summary["(mu) year","Mean"],         tolerance = 0.001)
  fit_mods.brma
})


test_that("Test against metafor::rma.ls", {

  skip_on_cran()
  skip_if_not_installed("metadat")
  skip_if_not_installed("metafor")

  ### fit location-scale meta-analytic model
  data(dat.bangertdrowns2004, package = "metadat")
  dat       <- dat.bangertdrowns2004
  dat$ni100 <- dat$ni/100
  dat$meta  <- as.factor(dat$meta)
  dat$imag  <- as.factor(dat$imag)
  fit_scale.metafor <- metafor::rma(yi, vi, mods = ~ ni100 + meta, scale = ~ ni100, data = dat)

  # using RoBMA package (using wide priors for scale to remove shrinkage)
  fit_scale.brma <- brma(yi, vi, mods = ~ ni100 + meta, scale = ~ ni100, data = dat, measure = "SMD", seed = 1,
                         prior_scale = list(ni100 = prior("normal", list(0, 1))))

  fit_scale.brma
  expect_equal(fit_scale.metafor$beta[[1]],  fit_scale.brma$summary["(mu) intercept","Mean"],  tolerance = 0.05)
  expect_equal(fit_scale.metafor$beta[[2]],  fit_scale.brma$summary["(mu) ni100","Mean"],      tolerance = 0.05)
  expect_equal(fit_scale.metafor$beta[[3]],  fit_scale.brma$summary["(mu) meta[1]","Mean"],    tolerance = 0.05)

  # TODO: the intercept is not adjusted for the descaling --- how to deal with this from BayesTools???
  expect_equal(exp(0.5 * fit_scale.metafor$alpha[[1]]), fit_scale.brma$summary["(tau) intercept","Mean"],  tolerance = 0.05)
  expect_equal(0.5 * fit_scale.metafor$alpha[[2]],  fit_scale.brma$summary["(tau_exp) ni100","Mean"],    tolerance = 0.05)
})
