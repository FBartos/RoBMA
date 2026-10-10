test_that("as_zplot requires a positive z threshold", {

  # a replay of a plain meta-analysis: 12 estimates with fixed posterior draws
  dat <- data.frame(
    yi  = c(.35, .10, .42, .28, .05, .31, .22, .47, .16, .38, .26, .12),
    sei = rep(.15, 12L)
  )
  object <- brma(
    yi = yi, sei = sei, data = dat, measure = "GEN", prior_unit_information_sd = 1,
    prior_effect = BayesTools::prior("normal", list(mean = 0, sd = 1)),
    prior_heterogeneity = BayesTools::prior(
      "normal", list(mean = 0, sd = 1), truncation = list(lower = 0, upper = Inf)
    ),
    only_priors = TRUE, silent = TRUE
  )
  withr::local_seed(2810)
  samples <- cbind(mu = stats::rnorm(200L, .25, .04), tau = abs(stats::rnorm(200L, .1, .03)))
  fit <- review_attach_replay(review_compile_object(object), samples)

  # a zero threshold has two-sided alpha 1 and expected discovery rate 1, so
  # that the Soric false discovery rate 0 * Inf is undefined
  for (value in c(0, -1)) {
    expect_error(
      as_zplot(fit, significance_level = value),
      "'significance_level' must be higher than 0"
    )
    expect_error(
      zplot(fit, significance_level = value),
      "'significance_level' must be higher than 0"
    )
  }

  # any positive z threshold is a valid reporting rule, including a tiny one
  # whose two-sided alpha is still below 1
  for (value in c(1e-8, .05, 1, stats::qnorm(.975))) {
    zfit <- as_zplot(fit, significance_level = value, max_samples = 100)
    expect_s3_class(zfit, "zplot_brma")
    expect_identical(zfit[["zplot"]][["data"]][["significance_level"]], value)
    estimates <- summary(zfit)[["estimates"]]
    expect_false(anyNA(estimates[, "Mean"]), info = format(value))
    expect_true(all(is.finite(as.matrix(estimates[, "Mean", drop = FALSE]))), info = format(value))
  }
})
