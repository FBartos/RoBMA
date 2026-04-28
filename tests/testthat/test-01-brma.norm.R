context("Model fitting for brma.norm")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
skip_on_cran()
skip_if_not_installed("metadat")
skip_if_not_installed("metafor")
skip_refit_if_cached("brma.norm")

### Uses examples from the metafor package
test_that("Test against metafor::rma.uni", {
  ### fit simple meta-analytic model
  data(dat.bcg, package = "metadat")
  dat <- metafor::escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg, data = dat.bcg)
  fit_simple.metafor <- metafor::rma(yi = yi, vi = vi, data = dat, method = "REML")

  # using RoBMA package
  fit_simple.brma <- brma(yi = yi, vi = vi, data = dat, measure = "RR", seed = 1, silent = TRUE)
  fit_simple.brma <- add_marglik(fit_simple.brma)
  fit_simple.brma <- add_loo(fit_simple.brma)
  save_fit("bcg_meta-analysis", fit_simple.brma, info = list(metafor = fit_simple.metafor))
  expect_equal(fit_simple.metafor$beta[[1]], fit_simple.brma$summary["mu", "Mean"], tolerance = 0.05)
  expect_equal(sqrt(fit_simple.metafor$tau2), fit_simple.brma$summary["tau", "Mean"], tolerance = 0.05)


  ### fit meta-regression (continuous predictors)
  fit_mods.metafor <- metafor::rma(yi, vi, mods = ~ ablat + year, data = dat)

  # using RoBMA package (rescale priors to reduce the shrinkage)
  fit_mods.brma <- brma(yi, vi, mods = ~ ablat + year, data = dat, measure = "RR", seed = 1, silent = TRUE)
  fit_mods.brma <- add_marglik(fit_mods.brma)
  fit_mods.brma <- suppressWarnings(add_loo(fit_mods.brma))
  save_fit("bcg_meta-regression", fit_mods.brma, info = list(mods = c("ablat", "year"), metafor = fit_mods.metafor))

  # expect_equal(fit_mods.metafor$beta[[1]], fit_mods.brma$summary["(mu) intercept", "Mean"], tolerance = 0.10)
  expect_equal(fit_mods.metafor$beta[[2]], fit_mods.brma$summary["(mu) ablat", "Mean"], tolerance = 0.10)
  expect_equal(fit_mods.metafor$beta[[3]], fit_mods.brma$summary["(mu) year", "Mean"], tolerance = 0.10)
  # the large difference in the intercept is a consequence in the large impact of small changes in the year coefficient
  # (they get multiplied by 2000)

  # verify that using scaled predictors removes the issue
  fit_mods.metafor.scaled <- metafor::rma(yi, vi, mods = ~ scale(ablat) + scale(year), data = dat)
  fit_mods.brma.scaled <- BayesTools::JAGS_estimates_table(fit_mods.brma$fit)
  expect_equal(fit_mods.metafor.scaled$beta[[1]], fit_mods.brma.scaled["(mu) intercept", "Mean"], tolerance = 0.10)
  expect_equal(fit_mods.metafor.scaled$beta[[2]], fit_mods.brma.scaled["(mu) ablat", "Mean"], tolerance = 0.10)
  expect_equal(fit_mods.metafor.scaled$beta[[3]], fit_mods.brma.scaled["(mu) year", "Mean"], tolerance = 0.10)


  ### fit meta-regression (factor predictor)
  fit_mods2.metafor <- metafor::rma(yi, vi, mods = ~ alloc, data = dat)

  # using RoBMA package (dummy coding)
  fit_mods2.brma <- brma(yi, vi, mods = ~ alloc, data = dat, measure = "RR", seed = 1, silent = TRUE)
  fit_mods2.brma <- add_marglik(fit_mods2.brma)
  fit_mods2.brma <- suppressWarnings(add_loo(fit_mods2.brma))
  save_fit("bcg_meta-regression2", fit_mods2.brma, info = list(mods = c("alloc"), metafor = fit_mods2.metafor))

  expect_equal(fit_mods2.metafor$beta[[1]], fit_mods2.brma$summary["(mu) intercept", "Mean"], tolerance = 0.05)
  expect_equal(fit_mods2.metafor$beta[[2]], fit_mods2.brma$summary["(mu) alloc[random]", "Mean"], tolerance = 0.15)
  expect_equal(fit_mods2.metafor$beta[[3]], fit_mods2.brma$summary["(mu) alloc[systematic]", "Mean"], tolerance = 0.05)

  # using RoBMA package (meandif coding)
  fit_mods2b.brma <- brma(yi, vi, mods = ~ alloc, data = dat, measure = "RR", seed = 1, silent = TRUE, rescale_priors = 2, set_contrast_factor_predictors = "meandif")
  fit_mods2b.brma <- add_marglik(fit_mods2b.brma)
  fit_mods2b.brma <- suppressWarnings(add_loo(fit_mods2b.brma))
  save_fit("bcg_meta-regression2b", fit_mods2b.brma, info = list(mods = c("alloc"), metafor = NULL))


  ### fit meta-regression with interaction (fact * cont)
  fit_mods3.metafor <- metafor::rma(yi, vi, mods = ~ alloc * year, data = dat)

  # using RoBMA package (increase scale because of shrinkage, dummy coding)
  fit_mods3.brma <- brma(yi, vi, mods = ~ alloc * year, data = dat, measure = "RR", seed = 1, silent = TRUE)
  fit_mods3.brma <- add_marglik(fit_mods3.brma)
  fit_mods3.brma <- suppressWarnings(add_loo(fit_mods3.brma))
  save_fit("bcg_meta-regression3", fit_mods3.brma, info = list(mods = c("alloc", "year", "alloc:year"), metafor = fit_mods3.metafor))

  # expect_equal(fit_mods3.metafor$beta[[1]], fit_mods3.brma$summary["(mu) intercept", "Mean"], tolerance = 0.15)
  # expect_equal(fit_mods3.metafor$beta[[2]], fit_mods3.brma$summary["(mu) alloc[random]", "Mean"], tolerance = 0.15)
  # expect_equal(fit_mods3.metafor$beta[[3]], fit_mods3.brma$summary["(mu) alloc[systematic]", "Mean"], tolerance = 0.15)
  expect_equal(fit_mods3.metafor$beta[[4]], fit_mods3.brma$summary["(mu) year", "Mean"], tolerance = 0.15)
  expect_equal(fit_mods3.metafor$beta[[5]], fit_mods3.brma$summary["(mu) alloc[systematic]:year", "Mean"], tolerance = 0.15)
  expect_equal(fit_mods3.metafor$beta[[6]], fit_mods3.brma$summary["(mu) alloc[random]:year", "Mean"], tolerance = 0.15)
  # lower levels are extremely shrunk due to sparse data, testing with wider priors shows correspondece

  # using RoBMA package (meandif coding)
  fit_mods3b.brma <- brma(yi, vi, mods = ~ alloc * year, data = dat, measure = "RR", seed = 1, silent = TRUE, set_contrast_factor_predictors = "meandif")
  fit_mods3b.brma <- add_marglik(fit_mods3b.brma)
  fit_mods3b.brma <- suppressWarnings(add_loo(fit_mods3b.brma))
  save_fit("bcg_meta-regression3b", fit_mods3b.brma, info = list(mods = c("alloc", "year", "alloc:year"), metafor = NULL))


  ### fit meta-regression with interaction (fact * fact)
  dat$year_before1969 <- factor(dat$year < 1969)
  fit_mods4.metafor <- metafor::rma(yi, vi, mods = ~ alloc * year_before1969, data = dat)

  # using RoBMA package (dummy coding)
  fit_mods4.brma <- brma(yi, vi, mods = ~ alloc * year_before1969, data = dat, measure = "RR", seed = 1, silent = TRUE)
  fit_mods4.brma <- add_marglik(fit_mods4.brma)
  fit_mods4.brma <- suppressWarnings(add_loo(fit_mods4.brma))
  save_fit("bcg_meta-regression4", fit_mods4.brma, info = list(mods = c("alloc", "year_before1969", "alloc:year_before1969"), metafor = fit_mods4.metafor))

  # expect_equal(fit_mods4.metafor$beta[[1]], fit_mods4.brma$summary["(mu) intercept", "Mean"], tolerance = 0.15)
  # expect_equal(fit_mods4.metafor$beta[[2]], fit_mods4.brma$summary["(mu) alloc[random]", "Mean"], tolerance = 0.15)
  expect_equal(fit_mods4.metafor$beta[[3]], fit_mods4.brma$summary["(mu) alloc[systematic]", "Mean"], tolerance = 0.15)
  # expect_equal(fit_mods4.metafor$beta[[4]], fit_mods4.brma$summary["(mu) year", "Mean"], tolerance = 0.15)
  # expect_equal(fit_mods4.metafor$beta[[5]], fit_mods4.brma$summary["(mu) alloc[random]:year_before1969[TRUE]", "Mean"], tolerance = 0.15)
  # expect_equal(fit_mods4.metafor$beta[[6]], fit_mods4.brma$summary["(mu) alloc[systematic]:year_before1969[TRUE]", "Mean"], tolerance = 0.15)
  # lower levels are extremely shrunk due to sparse data, testing with wider priors shows correspondece

  # using RoBMA package (increase scale because of shrinkage, meandif coding)
  fit_mods4b.brma <- brma(yi, vi, mods = ~ alloc * year_before1969, data = dat, measure = "RR", seed = 1, silent = TRUE, set_contrast_factor_predictors = "meandif")
  fit_mods4b.brma <- add_marglik(fit_mods4b.brma)
  fit_mods4b.brma <- suppressWarnings(add_loo(fit_mods4b.brma))
  save_fit("bcg_meta-regression4b", fit_mods4b.brma, info = list(mods = c("alloc", "year_before1969", "alloc:year_before1969"), metafor = NULL))
})

test_that("Test against metafor::rma.ls", {
  ### fit location-scale model
  data(dat.bangertdrowns2004, package = "metadat")
  dat <- dat.bangertdrowns2004
  dat$ni100 <- dat$ni / 100
  dat$meta <- as.factor(dat$meta)
  dat$imag <- as.factor(dat$imag)
  fit_scale.metafor <- suppressWarnings(metafor::rma(yi, vi, mods = ~ meta + ni100, scale = ~ni100, data = dat, method = "REML"))

  # using RoBMA package (using wide priors for scale to remove shrinkage -- the likelihood is very wide)
  fit_scale.brma <- suppressWarnings(brma(yi, vi, mods = ~ meta + ni100, scale = ~ni100, data = dat, measure = "SMD", seed = 1, silent = TRUE))
  fit_scale.brma <- add_marglik(fit_scale.brma)
  fit_scale.brma <- add_loo(fit_scale.brma)
  save_fit("bangertdrowns2004_location-scale", fit_scale.brma, info = list(mods = c("meta", "ni100"), scale = "ni100", metafor = fit_scale.metafor))

  expect_equal(fit_scale.metafor$beta[[1]], fit_scale.brma$summary["(mu) intercept", "Mean"], tolerance = 0.05)
  expect_equal(fit_scale.metafor$beta[[2]], fit_scale.brma$summary["(mu) meta[1]", "Mean"], tolerance = 0.05)
  expect_equal(fit_scale.metafor$beta[[3]], fit_scale.brma$summary["(mu) ni100", "Mean"], tolerance = 0.05)

  expect_equal(0.5 * fit_scale.metafor$alpha[[1]], log(fit_scale.brma$summary["(log_tau) exp(intercept)", "Mean"]), tolerance = 0.15)
  expect_equal(0.5 * fit_scale.metafor$alpha[[2]], fit_scale.brma$summary["(log_tau) ni100", "Mean"], tolerance = 0.15)
})

test_that("Test against metafor::rma.mv (3lvl)", {
  ### fit 3lvl model
  data(dat.konstantopoulos2011, package = "metadat")
  fit_3lvl.metafor <- metafor::rma.mv(yi, vi, random = ~ school | district, data = dat.konstantopoulos2011)

  # using RoBMA package
  fit_3lvl.brma <- brma(yi, vi, cluster = district, data = dat.konstantopoulos2011, measure = "SMD", seed = 1, silent = TRUE)
  fit_3lvl.brma <- add_marglik(fit_3lvl.brma)
  fit_3lvl.brma <- suppressWarnings(add_loo(fit_3lvl.brma))
  save_fit("konstantopoulos2011_3lvl", fit_3lvl.brma, info = list(metafor = fit_3lvl.metafor))

  expect_equal(fit_3lvl.metafor$beta[[1]], fit_3lvl.brma$summary["mu", "Mean"], tolerance = 0.01)
  expect_equal(sqrt(fit_3lvl.metafor$tau2), fit_3lvl.brma$summary["tau", "Mean"], tolerance = 0.01)
  expect_equal(fit_3lvl.metafor$rho, fit_3lvl.brma$summary["rho", "Mean"], tolerance = 0.05)

  ### fit 3lvl model with a predictor
  data(dat.konstantopoulos2011, package = "metadat")
  fit_3lvl2.metafor <- metafor::rma.mv(yi, vi, mods = ~ vi, random = ~ school | district, data = dat.konstantopoulos2011)

  # using RoBMA package
  fit_3lvl2.brma <- brma(yi, vi, mods = ~ vi, cluster = district, data = dat.konstantopoulos2011, measure = "SMD", seed = 1, silent = TRUE)
  fit_3lvl2.brma <- add_marglik(fit_3lvl2.brma)
  fit_3lvl2.brma <- suppressWarnings(add_loo(fit_3lvl2.brma))
  save_fit("konstantopoulos2011_3lvl2", fit_3lvl2.brma, info = list(metafor = fit_3lvl2.metafor))

  expect_equal(fit_3lvl2.metafor$beta[[1]], fit_3lvl2.brma$summary["(mu) intercept", "Mean"], tolerance = 0.10)
  expect_equal(fit_3lvl2.metafor$beta[[2]], fit_3lvl2.brma$summary["(mu) vi", "Mean"], tolerance = 0.10)
  expect_equal(sqrt(fit_3lvl2.metafor$tau2), fit_3lvl2.brma$summary["tau", "Mean"], tolerance = 0.10)
  expect_equal(fit_3lvl2.metafor$rho, fit_3lvl2.brma$summary["rho", "Mean"], tolerance = 0.10)
})
