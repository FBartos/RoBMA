context("Pooled prediction targets")


test_that("pooled effect summaries include prediction intervals", {

  object <- brma(
    yi                        = c(-0.2, 0.1, 0.4),
    sei                       = c(0.2, 0.2, 0.2),
    measure                   = "OR",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  posterior_samples <- cbind(
    mu  = c(-0.2, 0.0, 0.2, 0.4),
    tau = c(0.1, 0.2, 0.3, 0.4)
  )

  set.seed(481)
  expected_prediction <- matrix(
    stats::rnorm(
      n    = nrow(posterior_samples),
      mean = posterior_samples[, "mu"],
      sd   = posterior_samples[, "tau"]
    ),
    ncol = 1L
  )
  set.seed(481)
  pooled <- pooled_effect(
    object,
    .posterior_samples = posterior_samples
  )
  estimates <- summary(pooled)

  expect_equal(
    unname(attr(pooled, "prediction_samples")),
    unname(expected_prediction),
    tolerance = 0
  )
  expect_true(all(c("PI 0.025", "PI 0.975") %in% colnames(estimates)))
  expect_equal(
    unname(as.numeric(estimates["mu", c("PI 0.025", "PI 0.975")])),
    unname(stats::quantile(expected_prediction[, 1L], c(0.025, 0.975))),
    tolerance = 1e-12
  )

  set.seed(481)
  transformed <- pooled_effect(
    object,
    transform          = "EXP",
    .posterior_samples = posterior_samples
  )
  expect_equal(
    unname(as.matrix(transformed)),
    unname(exp(posterior_samples[, "mu", drop = FALSE])),
    tolerance = 1e-12
  )
  expect_equal(
    unname(attr(transformed, "prediction_samples")),
    unname(exp(expected_prediction)),
    tolerance = 1e-12
  )
})


test_that("pooled heterogeneity evaluates the average scale design", {

  dat <- data.frame(
    yi = c(-0.2, 0.1, 0.4),
    vi = c(0.04, 0.04, 0.04),
    x  = c(-1, 0, 2)
  )
  object <- brma(
    yi                        = yi,
    vi                        = vi,
    data                      = dat,
    scale                     = ~ x,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  posterior_samples <- cbind(
    mu                = c(0.1, 0.2),
    log_tau_intercept = c(0.2, 0.5),
    log_tau_x         = c(0.8, -0.6)
  )
  row_tau <- predict(
    object,
    type               = "terms.scale",
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )
  expected <- matrix(
    exp(rowMeans(log(as.matrix(row_tau)))),
    ncol     = 1L,
    dimnames = list(NULL, "tau")
  )
  rms <- matrix(sqrt(rowMeans(as.matrix(row_tau)^2)), ncol = 1L)

  pooled <- pooled_heterogeneity(
    object,
    .posterior_samples = posterior_samples
  )

  expect_equal(unname(as.matrix(pooled)), unname(expected), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(unname(expected), unname(rms))))
})
