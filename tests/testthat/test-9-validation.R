context("(9) Validation")
skip_on_cran()

test_that("Validate normal model with metafor", {

  library(metafor)
  dat <- escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg, data = dat.bcg)

  fit.metafor.fe <- rma(yi, vi, data = dat, method = "FE")
  fit.metafor.re <- rma(yi, vi, data = dat, method = "REML")

  # use RoBMA with very wide priors to get close to ML estimates
  fit.RoBMA <- NoBMA(y = dat$yi, se = sqrt(dat$vi),
                     priors_effect        = prior("normal", list(0, 5)),
                     priors_heterogeneity = prior("normal", list(0, 5), list(0, Inf)), seed = 1)


  fit.RoBMA.fe <- fit.RoBMA$models[[3]]$fit_summary
  fit.RoBMA.re <- fit.RoBMA$models[[4]]$fit_summary

  expect_equal(fit.RoBMA.fe["mu","Mean"], fit.metafor.fe$b[[1]], tolerance = 1e-3)
  expect_equal(fit.RoBMA.fe["mu","SD"],   fit.metafor.fe$se,     tolerance = 1e-3)

  expect_equal(fit.RoBMA.re["mu","Mean"], fit.metafor.re$b[[1]], tolerance = 1e-3)
  # expect_equal(fit.RoBMA.re["mu","SD"],   fit.metafor.re$se,     tolerance = 1e-2)
  # the standard error does not really match as metafor conditions on tau known (similarly for tau estimates)


  fit.metafor.fereg <- rma(yi, vi, mods = ~ factor(alloc) + scale(year), data = dat, method = "FE")
  fit.metafor.rereg <- rma(yi, vi, mods = ~ factor(alloc) + scale(year), data = dat, method = "REML")

  dat$y    <- dat$yi
  dat$se   <- sqrt(dat$vi)
  fit.RoBMA.reg <- suppressWarnings(NoBMA.reg(~ alloc + year, data = dat, autofit_control = set_autofit_control(min_ESS = 2000),
                                              prior_covariates     = prior("normal", list(0, 5)),
                                              prior_factors        = prior_factor("normal", list(0, 5), contrast = "treatment"),
                                              priors_effect        = prior("normal", list(0, 5)),
                                              priors_heterogeneity = prior("normal", list(0, 5), list(0, Inf)),
                                              test_predictors      = FALSE, priors_effect_null = NULL,
                                              prior_scale = "none", transformation = "none", standardize_predictors = TRUE, seed = 1))

  fit.RoBMA.reg.fe <- fit.RoBMA.reg$models[[1]]$fit_summary
  fit.RoBMA.reg.re <- fit.RoBMA.reg$models[[2]]$fit_summary

  expect_equal(fit.RoBMA.reg.fe["intercept","Mean"],         fit.metafor.fereg$b[[1]], tolerance = 1e-2)
  expect_equal(fit.RoBMA.reg.fe["alloc[random]","Mean"],     fit.metafor.fereg$b[[2]], tolerance = 1e-2)
  expect_equal(fit.RoBMA.reg.fe["alloc[systematic]","Mean"], fit.metafor.fereg$b[[3]], tolerance = 1e-2)
  expect_equal(fit.RoBMA.reg.fe["year","Mean"],              fit.metafor.fereg$b[[4]], tolerance = 1e-2)
  expect_equal(fit.RoBMA.reg.fe["intercept","SD"],         fit.metafor.fereg$se[[1]], tolerance = 1e-2)
  expect_equal(fit.RoBMA.reg.fe["alloc[random]","SD"],     fit.metafor.fereg$se[[2]], tolerance = 1e-2)
  expect_equal(fit.RoBMA.reg.fe["alloc[systematic]","SD"], fit.metafor.fereg$se[[3]], tolerance = 1e-2)
  expect_equal(fit.RoBMA.reg.fe["year","SD"],              fit.metafor.fereg$se[[4]], tolerance = 1e-2)

  expect_equal(fit.RoBMA.reg.re["intercept","Mean"],         fit.metafor.rereg$b[[1]], tolerance = 2e-2)
  expect_equal(fit.RoBMA.reg.re["alloc[random]","Mean"],     fit.metafor.rereg$b[[2]], tolerance = 2e-2)
  expect_equal(fit.RoBMA.reg.re["alloc[systematic]","Mean"], fit.metafor.rereg$b[[3]], tolerance = 2e-2)
  expect_equal(fit.RoBMA.reg.re["year","Mean"],              fit.metafor.rereg$b[[4]], tolerance = 2e-2)
  # the standard error does not really match as metafor conditions on tau known (similarly for tau estimates)

})

test_that("Validate binomial model with metafor", {


  library(metafor)
  fit.metafor.re <- suppressWarnings(metafor::rma.glmm(measure = "OR", ai = dat.anand1999$ai, ci = dat.anand1999$ci, n1i = dat.anand1999$n1i, n2i = dat.anand1999$n2i, model = "UM.FS"))

  # use RoBMA with very wide priors to get close to ML estimates
  fit.RoBMA <- BiBMA(x1 = dat.anand1999$ai, x2 = dat.anand1999$ci, n1 = dat.anand1999$n1i, n2 =  dat.anand1999$n2i,
                     priors_effect_null = NULL, priors_heterogeneity_null = NULL,
                     priors_effect        = prior("normal", list(0, 1)),
                     priors_heterogeneity = prior("normal", list(0, 1), list(0, Inf)), seed = 1)

  fit.RoBMA.re <- fit.RoBMA$models[[1]]$fit_summary

  expect_equal(fit.RoBMA.re["mu","Mean"], fit.metafor.re$b[[1]], tolerance = 1e-3)

})
