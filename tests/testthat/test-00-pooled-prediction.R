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
  estimates_df <- as.data.frame(estimates)

  expect_equal(
    unname(attr(pooled, "prediction_samples")),
    unname(expected_prediction),
    tolerance = 0
  )
  expect_true(all(c("PI 0.025", "PI 0.975") %in% colnames(estimates)))
  pooled_df <- as.data.frame(pooled)
  expect_identical(
    names(pooled_df),
    c(
      "component", "parameter", "Mean", "Median",
      "CI_0.025", "CI_0.975", "PI_0.025", "PI_0.975"
    )
  )
  expect_identical(pooled_df[["component"]], "location")
  expect_identical(pooled_df[["parameter"]], "mu")
  expect_identical(data.frame(pooled), pooled_df)
  expect_equal(
    unname(as.matrix(pooled_df)),
    unname(as.matrix(estimates_df)),
    tolerance = 0
  )
  expect_identical(
    names(as.data.frame(pooled, probs = c(.05, .95))),
    c(
      "component", "parameter", "Mean", "Median",
      "CI_0.05", "CI_0.95", "PI_0.05", "PI_0.95"
    )
  )
  pooled_draws <- as.matrix(RoBMA::as_draws_matrix(pooled))
  expect_equal(
    as.vector(pooled_draws),
    as.vector(as.matrix(pooled)),
    tolerance = 0
  )

  custom_probs <- c(.6, .125, .9)
  custom_pooled <- pooled_effect(
    object,
    probs              = custom_probs,
    .posterior_samples = posterior_samples
  )
  custom_summary <- summary(custom_pooled)
  custom_df      <- as.data.frame(custom_pooled)
  expect_identical(
    names(custom_summary),
    c("Mean", "Median", "0.6", "0.125", "0.9",
      "PI 0.6", "PI 0.125", "PI 0.9")
  )
  expect_identical(
    names(custom_df),
    c(
      "component", "parameter", "Mean", "Median",
      "CI_0.6", "CI_0.125", "CI_0.9",
      "PI_0.6", "PI_0.125", "PI_0.9"
    )
  )
  expect_identical(data.frame(custom_pooled), custom_df)
  expect_equal(
    unname(as.numeric(custom_df[1L, paste0("CI_", custom_probs)])),
    unname(stats::quantile(as.matrix(custom_pooled)[, 1L], custom_probs)),
    tolerance = 0
  )
  expect_equal(
    unname(as.numeric(custom_df[1L, paste0("PI_", custom_probs)])),
    unname(stats::quantile(
      attr(custom_pooled, "prediction_samples")[, 1L],
      custom_probs
    )),
    tolerance = 0
  )
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
  pooled_df  <- as.data.frame(pooled)
  row_tau_df <- as.data.frame(row_tau)
  expect_identical(
    names(pooled_df),
    c("component", "parameter", "Mean", "Median", "CI_0.025", "CI_0.975")
  )
  expect_identical(pooled_df[["component"]], "heterogeneity")
  expect_identical(pooled_df[["parameter"]], "tau")
  expect_identical(data.frame(pooled), pooled_df)
  expect_equal(
    unname(as.matrix(pooled_df)),
    unname(as.matrix(as.data.frame(summary(pooled)))),
    tolerance = 0
  )
  expect_identical(
    names(row_tau_df),
    c("component", "parameter", "Mean", "Median", "CI_0.025", "CI_0.975")
  )
  expect_identical(row_tau_df[["component"]], rep("scale", 3L))
  expect_identical(data.frame(row_tau), row_tau_df)
  expect_equal(
    unname(as.matrix(row_tau_df)),
    unname(as.matrix(as.data.frame(summary(row_tau)))),
    tolerance = 0
  )
  expect_false(isTRUE(all.equal(unname(expected), unname(rms))))
})
