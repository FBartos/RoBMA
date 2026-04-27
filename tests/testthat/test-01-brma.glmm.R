context("Model fitting for brma.glmm")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
skip_on_cran()
skip_if_not_installed("metadat")
skip_if_not_installed("metafor")
skip_refit_if_cached("brma.glmm")

### Uses examples from the metafor package
test_that("Test against metafor::rma.glmm", {
  ### fit generalized meta-analytic model to difference in two proportions
  data(dat.bcg, package = "metadat")
  fit_simple.metafor <- metafor::rma.glmm(measure = "OR", ai = tpos, bi = tneg, ci = cpos, di = cneg, data = dat.bcg, model = "UM.FS")

  # using RoBMA package
  fit_simple.brma <- brma.glmm(ai = tpos, bi = tneg, ci = cpos, di = cneg, data = dat.bcg, measure = "OR", seed = 1, silent = TRUE)
  fit_simple.brma <- add_marglik(fit_simple.brma)
  fit_simple.brma <- suppressWarnings(add_loo(fit_simple.brma))
  save_fit("bcg_glmm", fit_simple.brma, info = list(metafor = fit_simple.metafor))
  expect_equal(fit_simple.metafor$beta[[1]], fit_simple.brma$summary["mu", "Mean"], tolerance = 0.05)
  expect_equal(sqrt(fit_simple.metafor$tau2), fit_simple.brma$summary["tau", "Mean"], tolerance = 0.05)


  ### fit generalized meta-regression
  fit_reg.metafor <- suppressWarnings(metafor::rma.glmm(measure = "OR", ai = tpos, bi = tneg, ci = cpos, di = cneg, mods = ~ alloc, data = dat.bcg, model = "UM.FS"))

  # using RoBMA package
  fit_reg.brma <- brma.glmm(ai = tpos, bi = tneg, ci = cpos, di = cneg, mods = ~ alloc, data = dat.bcg, measure = "OR", seed = 1, silent = TRUE)
  fit_reg.brma <- add_marglik(fit_reg.brma)
  fit_reg.brma <- suppressWarnings(add_loo(fit_reg.brma))
  save_fit("bcg_glmm_reg", fit_reg.brma, info = list(mods = c("ablat"), metafor = fit_reg.metafor))
  expect_equal(fit_reg.metafor$beta[[1]], fit_reg.brma$summary["(mu) intercept", "Mean"], tolerance = 0.05)
  expect_equal(fit_reg.metafor$beta[[2]], fit_reg.brma$summary["(mu) alloc[random]", "Mean"], tolerance = 0.15)
  expect_equal(fit_reg.metafor$beta[[3]], fit_reg.brma$summary["(mu) alloc[systematic]", "Mean"], tolerance = 0.10)
  expect_equal(sqrt(fit_reg.metafor$tau2), fit_reg.brma$summary["tau", "Mean"], tolerance = 0.10)


  ### fit generalized meta-analytic model to difference in two rations
  data(dat.nielweise2008, package = "metadat")
  fit_simple.metafor <- metafor::rma.glmm(measure = "IRR", x1i = x1i, t1i = t1i, x2i = x2i, t2i = t2i, data = dat.nielweise2008, model = "UM.FS")

  fit_simple.brma <- brma.glmm(x1i = x1i, t1i = t1i, x2i = x2i, t2i = t2i, data = dat.nielweise2008, measure = "IRR", seed = 1, silent = TRUE)
  fit_simple.brma <- add_marglik(fit_simple.brma)
  fit_simple.brma <- suppressWarnings(add_loo(fit_simple.brma))
  save_fit("nielweise2008_glmm", fit_simple.brma, info = list(metafor = fit_simple.metafor))
  expect_equal(fit_simple.metafor$beta[[1]], fit_simple.brma$summary["mu", "Mean"], tolerance = 0.05)
  expect_equal(sqrt(fit_simple.metafor$tau2), fit_simple.brma$summary["tau", "Mean"], tolerance = 0.10) # the tau is very variable here
})

test_that("brma.glmm handles multilevel scale regression model", {
  # using RoBMA package
  data(dat.bcg, package = "metadat")
  fit_simple.brma <- brma.glmm(ai = tpos, bi = tneg, ci = cpos, di = cneg, scale = ~ year, cluster = alloc, data = dat.bcg, measure = "OR", seed = 1, silent = TRUE)
  fit_simple.brma <- add_marglik(fit_simple.brma)
  fit_simple.brma <- suppressWarnings(add_loo(fit_simple.brma))
  save_fit("bcg_glmm_3lvl_scale", fit_simple.brma, info = list(scale = c("year")))

})

# TODO: add single-group models
