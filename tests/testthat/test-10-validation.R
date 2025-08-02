context("(10) Validation")
skip_on_cran()

test_that("Validate normal model with metafor", {

  library(metafor)
  dat <- escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg, data = dat.bcg)

  fit.metafor.fe <- rma(yi, vi, data = dat, method = "FE")
  fit.metafor.re <- rma(yi, vi, data = dat, method = "REML")

  # use RoBMA with very wide priors to get close to ML estimates
  fit.RoBMA <- NoBMA(y = dat$yi, se = sqrt(dat$vi),
                     priors_effect        = prior("normal", list(0, 5)),
                     priors_heterogeneity = prior("normal", list(0, 5), list(0, Inf)), seed = 1,
                     algorithm = "bridge", chains = 2, sample = 1000, burnin = 250, adapt = 250, thin = 1, parallel = FALSE)

  fit.RoBMA_ss.fe <- NoBMA(y = dat$yi, se = sqrt(dat$vi),
                           priors_effect        = prior("normal", list(0, 5)),
                           priors_heterogeneity = NULL, seed = 1,
                           priors_effect_null = NULL,
                           algorithm = "ss", chains = 2, sample = 2000, burnin = 500, adapt = 500, thin = 1, parallel = FALSE)
  fit.RoBMA_ss.re <- NoBMA(y = dat$yi, se = sqrt(dat$vi),
                     priors_effect        = prior("normal", list(0, 5)),
                     priors_heterogeneity = prior("normal", list(0, 5), list(0, Inf)), seed = 1,
                     priors_effect_null = NULL, priors_heterogeneity_null = NULL,
                     algorithm = "ss", chains = 2, sample = 2000, burnin = 500, adapt = 500, thin = 1, parallel = FALSE)


  fit.RoBMA.fe <- fit.RoBMA$models[[3]]$fit_summary
  fit.RoBMA.re <- fit.RoBMA$models[[4]]$fit_summary

  expect_equal(fit.RoBMA.fe["mu","Mean"], fit.metafor.fe$b[[1]], tolerance = 1e-3)
  expect_equal(fit.RoBMA.fe["mu","SD"],   fit.metafor.fe$se,     tolerance = 1e-3)

  expect_equal(fit.RoBMA.re["mu","Mean"], fit.metafor.re$b[[1]], tolerance = 1e-2)
  # expect_equal(fit.RoBMA.re["mu","SD"],   fit.metafor.re$se,     tolerance = 1e-2)
  # the standard error does not really match as metafor conditions on tau known (similarly for tau estimates)


  fit.metafor.fereg <- rma(yi, vi, mods = ~ factor(alloc) + scale(year), data = dat, method = "FE")
  fit.metafor.rereg <- rma(yi, vi, mods = ~ factor(alloc) + scale(year), data = dat, method = "REML")

  dat$y    <- dat$yi
  dat$se   <- sqrt(dat$vi)
  fit.RoBMA.reg.fe <- suppressWarnings(NoBMA.reg(~ alloc + year, data = dat,
                                                 prior_covariates     = prior("normal", list(0, 5)),
                                                 prior_factors        = prior_factor("normal", list(0, 5), contrast = "treatment"),
                                                 priors_effect        = prior("normal", list(0, 5)),
                                                 priors_heterogeneity = NULL,
                                                 test_predictors      = FALSE, priors_effect_null = NULL,
                                                 prior_scale = "none", transformation = "none", standardize_predictors = TRUE, seed = 1,
                                                 algorithm = "ss", chains = 2, sample = 2000, burnin = 500, adapt = 500, parallel = FALSE))
  fit.RoBMA.reg.re <- suppressWarnings(NoBMA.reg(~ alloc + year, data = dat,
                                              prior_covariates     = prior("normal", list(0, 5)),
                                              prior_factors        = prior_factor("normal", list(0, 5), contrast = "treatment"),
                                              priors_effect        = prior("normal", list(0, 5)),
                                              priors_heterogeneity = prior("normal", list(0, 5), list(0, Inf)),
                                              test_predictors      = FALSE, priors_effect_null = NULL,
                                              priors_heterogeneity_null = NULL,
                                              prior_scale = "none", transformation = "none", standardize_predictors = TRUE, seed = 1,
                                              algorithm = "ss", chains = 2, sample = 2000, burnin = 500, adapt = 500, parallel = FALSE))

  sum.RoBMA.reg.fe <- summary(fit.RoBMA.reg.fe)$estimates_predictors
  sum.RoBMA.reg.re <- summary(fit.RoBMA.reg.re)$estimates_predictors

  expect_equal(sum.RoBMA.reg.fe["intercept","Mean"],         fit.metafor.fereg$b[[1]], tolerance = 1e-2)
  expect_equal(sum.RoBMA.reg.fe["alloc[random]","Mean"],     fit.metafor.fereg$b[[2]], tolerance = 1e-2)
  expect_equal(sum.RoBMA.reg.fe["alloc[systematic]","Mean"], fit.metafor.fereg$b[[3]], tolerance = 1e-2)
  expect_equal(sum.RoBMA.reg.fe["year","Mean"],              fit.metafor.fereg$b[[4]], tolerance = 1e-2)

  expect_equal(sum.RoBMA.reg.re["intercept","Mean"],         fit.metafor.rereg$b[[1]], tolerance = 2e-2)
  expect_equal(sum.RoBMA.reg.re["alloc[random]","Mean"],     fit.metafor.rereg$b[[2]], tolerance = 2e-2)
  expect_equal(sum.RoBMA.reg.re["alloc[systematic]","Mean"], fit.metafor.rereg$b[[3]], tolerance = 2.1e-2)
  expect_equal(sum.RoBMA.reg.re["year","Mean"],              fit.metafor.rereg$b[[4]], tolerance = 2e-2)

  # Test standardized_coefficients = FALSE for raw coefficients
  fit.metafor.fereg2 <- rma(yi, vi, mods = ~ factor(alloc) + year, data = dat, method = "FE")
  fit.metafor.rereg2 <- rma(yi, vi, mods = ~ factor(alloc) + year, data = dat, method = "REML")

  sum.RoBMA.reg.fe.raw <- summary(fit.RoBMA.reg.fe, standardized_coefficients = FALSE)$estimates_predictors
  sum.RoBMA.reg.re.raw <- summary(fit.RoBMA.reg.re, standardized_coefficients = FALSE)$estimates_predictors

  # For raw coefficients, only intercept and continuous predictor (year) should be transformed
  expect_equal(sum.RoBMA.reg.fe.raw["intercept","Mean"],         fit.metafor.fereg2$b[[1]], tolerance = 1e-2)
  expect_equal(sum.RoBMA.reg.fe.raw["alloc[random]","Mean"],     fit.metafor.fereg2$b[[2]], tolerance = 1e-2)
  expect_equal(sum.RoBMA.reg.fe.raw["alloc[systematic]","Mean"], fit.metafor.fereg2$b[[3]], tolerance = 1e-2)
  expect_equal(sum.RoBMA.reg.fe.raw["year","Mean"],              fit.metafor.fereg2$b[[4]], tolerance = 1e-2)

  expect_equal(sum.RoBMA.reg.re.raw["intercept","Mean"],         fit.metafor.rereg2$b[[1]], tolerance = 2.5e-2)
  expect_equal(sum.RoBMA.reg.re.raw["alloc[random]","Mean"],     fit.metafor.rereg2$b[[2]], tolerance = 2.5e-2)
  expect_equal(sum.RoBMA.reg.re.raw["alloc[systematic]","Mean"], fit.metafor.rereg2$b[[3]], tolerance = 2.5e-2)
  expect_equal(sum.RoBMA.reg.re.raw["year","Mean"],              fit.metafor.rereg2$b[[4]], tolerance = 2.5e-2)

  ## multiple continuous predictors
  fit.metafor.fereg3std <- rma(yi, vi, mods = ~ factor(alloc) + scale(year) + scale(tpos), data = dat, method = "FE")
  fit.metafor.fereg3raw <- rma(yi, vi, mods = ~ factor(alloc) + year + tpos, data = dat, method = "FE")
  fit.RoBMA.reg3.fe     <- suppressWarnings(NoBMA.reg(~ alloc + year+ tpos, data = dat,
                                                 prior_covariates     = prior("normal", list(0, 5)),
                                                 prior_factors        = prior_factor("normal", list(0, 5), contrast = "treatment"),
                                                 priors_effect        = prior("normal", list(0, 5)),
                                                 priors_heterogeneity = NULL,
                                                 test_predictors      = FALSE, priors_effect_null = NULL,
                                                 prior_scale = "none", transformation = "none", standardize_predictors = TRUE, seed = 1,
                                                 algorithm = "ss", chains = 2, sample = 2000, burnin = 500, adapt = 500, parallel = FALSE))

  sum.RoBMA.reg3std.fe <- summary(fit.RoBMA.reg3.fe, standardized_coefficients = TRUE)$estimates_predictors
  sum.RoBMA.reg3raw.fe <- summary(fit.RoBMA.reg3.fe, standardized_coefficients = FALSE)$estimates_predictors

  expect_equal(sum.RoBMA.reg3std.fe["intercept","Mean"],         fit.metafor.fereg3std$b[[1]], tolerance = 2e-2)
  expect_equal(sum.RoBMA.reg3std.fe["alloc[random]","Mean"],     fit.metafor.fereg3std$b[[2]], tolerance = 1e-2)
  expect_equal(sum.RoBMA.reg3std.fe["alloc[systematic]","Mean"], fit.metafor.fereg3std$b[[3]], tolerance = 1e-2)
  expect_equal(sum.RoBMA.reg3std.fe["year","Mean"],              fit.metafor.fereg3std$b[[4]], tolerance = 1e-2)
  expect_equal(sum.RoBMA.reg3std.fe["tpos","Mean"],              fit.metafor.fereg3std$b[[5]], tolerance = 1e-2)

  expect_equal(sum.RoBMA.reg3raw.fe["intercept","Mean"],         fit.metafor.fereg3raw$b[[1]], tolerance = 2e-2)
  expect_equal(sum.RoBMA.reg3raw.fe["alloc[random]","Mean"],     fit.metafor.fereg3raw$b[[2]], tolerance = 1e-2)
  expect_equal(sum.RoBMA.reg3raw.fe["alloc[systematic]","Mean"], fit.metafor.fereg3raw$b[[3]], tolerance = 1e-2)
  expect_equal(sum.RoBMA.reg3raw.fe["year","Mean"],              fit.metafor.fereg3raw$b[[4]], tolerance = 1e-2)
  expect_equal(sum.RoBMA.reg3raw.fe["tpos","Mean"],              fit.metafor.fereg3raw$b[[5]], tolerance = 1e-2)

  ### validate BLUPs ----
  robma_blups   <- true_effects(fit.RoBMA_ss.fe)
  metafor_blups <- blup(fit.metafor.fe)
  expect_equal(robma_blups$estimates[,"Mean"], metafor_blups$pred, tolerance = 1e-2)

  robma_blups   <- true_effects(fit.RoBMA_ss.re)
  metafor_blups <- blup(fit.metafor.re)
  expect_equal(robma_blups$estimates[,"Mean"], metafor_blups$pred, tolerance = 4e-2)

  robma_blups   <- true_effects(fit.RoBMA.reg.re)
  metafor_blups <- blup(fit.metafor.rereg2)
  expect_equal(robma_blups$estimates[,"Mean"], metafor_blups$pred, tolerance = 2e-2)

  robma_blups   <- true_effects(fit.RoBMA.reg.fe)
  metafor_blups <- blup(fit.metafor.fereg2)
  expect_equal(robma_blups$estimates[,"Mean"], metafor_blups$pred, tolerance = 1e-2)

  ### validate predictions ----
  set.seed(1)
  robma_preds   <- predict(fit.RoBMA_ss.fe, type = "effect")
  metafor_preds <- predict(fit.metafor.fe)
  expect_equal(as.numeric(robma_preds$estimates[1,c("Mean", "0.025", "0.975")]), c(metafor_preds$pred, metafor_preds$ci.lb, metafor_preds$ci.ub), tolerance = 1e-2)

  robma_preds   <- predict(fit.RoBMA_ss.re, type = "effect")
  metafor_preds <- predict(fit.metafor.re)
  # Bayesian tau estimate is larger -> credibile/prediction intervals mismatch too
  expect_equal(as.numeric(robma_preds$estimates[1,c("Mean")]), c(metafor_preds$pred), tolerance = 1e-2)

  robma_preds   <- predict(fit.RoBMA.reg.re, type = "effect")
  metafor_preds <- predict(fit.metafor.rereg2)
  expect_equal(as.numeric(robma_preds$estimates[,c("Mean")]), c(metafor_preds$pred), tolerance = 2e-2)

  robma_preds   <- predict(fit.RoBMA.reg.fe, type = "effect")
  metafor_preds <- predict(fit.metafor.fereg2)
  expect_equal(as.numeric(robma_preds$estimates[,c("Mean")]), c(metafor_preds$pred), tolerance = 1e-2)
  expect_equal(as.numeric(robma_preds$estimates[,c("0.025")]), c(metafor_preds$ci.lb), tolerance = 1e-2)
  expect_equal(as.numeric(robma_preds$estimates[,c("0.975")]), c(metafor_preds$ci.ub), tolerance = 1e-2)

  robma_preds.terms    <- predict(fit.RoBMA.reg.re, newdata = data.frame(alloc = "random", year = 1950:2000, y = 0, se = 0.50), type = "terms")
  robma_preds.effect   <- predict(fit.RoBMA.reg.re, newdata = data.frame(alloc = "random", year = 1950:2000, y = 0, se = 0.50), type = "effect")
  robma_preds.response <- predict(fit.RoBMA.reg.re, newdata = data.frame(alloc = "random", year = 1950:2000, y = 0, se = 0.50), type = "response")

  vdiffr::expect_doppelganger("predict-distribution", function(){
    plot(NA, xlim = c(1950, 2000), ylim = c(-4, 4), xlab = "Year", ylab = "Effect", las = 1)
    polygon(c(1950:2000, rev(1950:2000)),
            c(robma_preds.terms$estimates[,"0.025"], rev(robma_preds.terms$estimates[,"0.975"])),
            col = rgb(0, 0, 1, alpha = 0.2), border = NA)
    polygon(c(1950:2000, rev(1950:2000)),
            c(robma_preds.effect$estimates[,"0.025"], rev(robma_preds.effect$estimates[,"0.975"])),
            col = rgb(0, 0, 1, alpha = 0.2), border = NA)
    polygon(c(1950:2000, rev(1950:2000)),
            c(robma_preds.response$estimates[,"0.025"], rev(robma_preds.response$estimates[,"0.975"])),
            col = rgb(0, 0, 1, alpha = 0.2), border = NA)
  })


  ### validate residuals ----
  expect_equal(unname(residuals(fit.metafor.fe)), residuals(fit.RoBMA_ss.fe)$estimates[,"Mean"], tolerance = 1e-3)
  expect_equal(unname(residuals(fit.metafor.re)), residuals(fit.RoBMA_ss.re)$estimates[,"Mean"], tolerance = 1e-3)
  expect_equal(unname(residuals(fit.metafor.fereg2)), residuals(fit.RoBMA.reg.fe)$estimates[,"Mean"], tolerance = 1e-2)
  expect_equal(unname(residuals(fit.metafor.rereg2)), residuals(fit.RoBMA.reg.re)$estimates[,"Mean"], tolerance = 2.5e-2)

})

test_that("Validate binomial model with metafor", {


  library(metafor)
  fit.metafor.re <- suppressWarnings(metafor::rma.glmm(measure = "OR", ai = dat.anand1999$ai, ci = dat.anand1999$ci, n1i = dat.anand1999$n1i, n2i = dat.anand1999$n2i, model = "UM.FS"))

  # use RoBMA with very wide priors to get close to ML estimates
  fit.RoBMA <- BiBMA(x1 = dat.anand1999$ai, x2 = dat.anand1999$ci, n1 = dat.anand1999$n1i, n2 =  dat.anand1999$n2i,
                     priors_effect_null = NULL, priors_heterogeneity_null = NULL,
                     priors_effect        = prior("normal", list(0, 1)),
                     priors_heterogeneity = prior("normal", list(0, 1), list(0, Inf)), seed = 1,
                     algorithm = "ss", chains = 2, sample = 2000, burnin = 500, adapt = 500, parallel = FALSE)

  fit.RoBMA.re <- summary(fit.RoBMA)$estimates

  expect_equal(fit.RoBMA.re["mu","Mean"], fit.metafor.re$b[[1]], tolerance = 1e-2)


  # regression
  df <- dat.anand1999
  fit.metafor.reg  <- suppressWarnings(metafor::rma.glmm(measure = "OR", ai = df$ai, ci = df$ci, n1i = df$n1i, n2i = df$n2i,
                                                       mods = ~ scale(year), data = df, model = "UM.FS"))
  fit.metafor.reg2 <- suppressWarnings(metafor::rma.glmm(measure = "OR", ai = df$ai, ci = df$ci, n1i = df$n1i, n2i = df$n2i,
                                                        mods = ~ year, data = df, model = "UM.FS"))

  dfb <- data.frame(
    x1 = df$ai,
    x2 = df$ci,
    n1 = df$n1i,
    n2 = df$n2i,
    year = df$year
  )
  fit.BiBMA.reg <- BiBMA.reg(formula = ~ year, data = dfb,
                             prior_covariates     = prior("normal", list(0, 5)),
                             priors_effect        = prior("normal", list(0, 5)),
                             priors_heterogeneity = prior("normal", list(0, 5), list(0, Inf)),
                             test_predictors      = FALSE, priors_effect_null = NULL,
                             priors_heterogeneity_null = NULL,
                             seed = 1,
                             algorithm = "ss", chains = 2, sample = 2000, burnin = 1000, adapt = 1000, parallel = FALSE)

  sum.BiBMA.reg     <- summary(fit.BiBMA.reg)$estimates_predictors
  sum.BiBMA.reg.raw <- summary(fit.BiBMA.reg, standardized_coefficients = FALSE)$estimates_predictors

  expect_equal(sum.BiBMA.reg["intercept","Mean"],  fit.metafor.reg$b[[1]], tolerance = 2.5e-2)
  expect_equal(sum.BiBMA.reg["year","Mean"],       fit.metafor.reg$b[[2]], tolerance = 2.5e-2)

  expect_equal(sum.BiBMA.reg.raw["intercept","Mean"],  fit.metafor.reg2$b[[1]], tolerance = 2e-1)
  expect_equal(sum.BiBMA.reg.raw["year","Mean"],       fit.metafor.reg2$b[[2]], tolerance = 2.5e-2)
})
