context("Summary heterogeneity helpers")

test_that("scale heterogeneity summaries aggregate variances before tau", {

  tau_within <- matrix(
    c(1, 3,
      2, 4),
    nrow = 2,
    byrow = TRUE
  )

  samples <- RoBMA:::.summary_heterogeneity_samples(
    tau_within_samples  = tau_within,
    tau_between_samples = matrix(0, nrow = 2, ncol = 2),
    v_tilde             = 1,
    is_multilevel       = FALSE
  )

  expected_tau2 <- rowMeans(tau_within^2)
  old_tau2      <- rowMeans(tau_within)^2

  expect_equal(samples[["tau2"]], expected_tau2)
  expect_equal(samples[["tau"]], sqrt(expected_tau2))
  expect_false(isTRUE(all.equal(samples[["tau2"]], old_tau2)))
  expect_equal(samples[["I2"]], rowMeans(100 * tau_within^2 / (tau_within^2 + 1)))
  expect_equal(samples[["H2"]], rowMeans(tau_within^2 + 1))
})

test_that("multilevel scale heterogeneity partitions variance and I2", {

  tau_within <- matrix(
    c(1, 2,
      3, 4),
    nrow = 2,
    byrow = TRUE
  )
  tau_between <- matrix(
    c(2, 1,
      1, 3),
    nrow = 2,
    byrow = TRUE
  )

  samples <- RoBMA:::.summary_heterogeneity_samples(
    tau_within_samples  = tau_within,
    tau_between_samples = tau_between,
    rho_samples         = c(0.8, 0.2),
    v_tilde             = 2,
    is_multilevel       = TRUE
  )

  sigma2_within  <- tau_within^2
  sigma2_between <- tau_between^2
  sigma2_total   <- sigma2_within + sigma2_between
  denominator    <- sigma2_total + 2

  expect_equal(samples[["tau2 [within]"]], rowMeans(sigma2_within))
  expect_equal(samples[["tau2 [between]"]], rowMeans(sigma2_between))
  expect_equal(samples[["tau2"]], rowMeans(sigma2_total))
  expect_equal(samples[["rho"]], c(0.8, 0.2))
  expect_equal(samples[["I2"]], rowMeans(100 * sigma2_total / denominator))
  expect_equal(samples[["I2 [within]"]], rowMeans(100 * sigma2_within / denominator))
  expect_equal(samples[["I2 [between]"]], rowMeans(100 * sigma2_between / denominator))
  expect_equal(samples[["I2"]], samples[["I2 [within]"]] + samples[["I2 [between]"]])
})
