source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-contracts.R"))

skip_if_missing_fits("bselmodel.mv_exact_random")
fit_bselmodel_mv <- load_fit("bselmodel.mv_exact_random", validate = FALSE)


test_that("bselmodel.mv presents multivariate selection summaries", {

  out <- summary(
    fit_bselmodel_mv,
    include_mcmc_diagnostics = FALSE
  )
  frame <- as.data.frame(out)

  expect_s3_class(out, "summary.brma")
  expect_match(out[["name"]], "Multivariate.*Exact Selection")
  expect_true(nrow(out[["estimates_random"]]) > 0L)
  expect_true(nrow(out[["estimates_bias"]]) > 0L)
  expect_true(all(c("random", "bias") %in% frame[["component"]]))
})


test_that("bselmodel.mv likelihood methods preserve the exact targets", {

  draws   <- .get_posterior_samples(fit_bselmodel_mv[["fit"]])
  log_lik <- log_lik(fit_bselmodel_mv, unit = "estimate")
  target  <- attr(log_lik, "RoBMA_target", exact = TRUE)

  expect_identical(dim(log_lik), c(nrow(draws), nobs(fit_bselmodel_mv)))
  expect_true(all(is.finite(log_lik)))
  expect_identical(target[["unit"]], "estimate")
  expect_identical(target[["retained_context"]], "remaining_data")
  expect_s3_class(loo(fit_bselmodel_mv, unit = "estimate"), "loo")
  expect_true(is.finite(as.numeric(logml(fit_bselmodel_mv))))
})


test_that("bselmodel.mv prediction and latent summaries remain available", {

  n     <- nobs(fit_bselmodel_mv)
  draws <- .get_posterior_samples(fit_bselmodel_mv[["fit"]])

  estimate_marginal <- predict(
    fit_bselmodel_mv,
    type               = "estimate",
    conditioning_depth = "marginal"
  )
  estimate_fitted <- predict(
    fit_bselmodel_mv,
    type               = "estimate",
    conditioning_depth = "estimate"
  )
  response <- predict(
    fit_bselmodel_mv,
    type          = "response",
    bias_adjusted = FALSE
  )
  newdata <- data.frame(
    study = factor(c("new study", "new study")),
    x     = c(0, 1)
  )
  V_new <- matrix(c(0.02, 0.006, 0.006, 0.03), nrow = 2L)
  response_new <- predict(
    fit_bselmodel_mv,
    newdata       = newdata,
    V_new         = V_new,
    type          = "response",
    bias_adjusted = FALSE
  )

  expect_brma_samples_matrix(estimate_marginal, n, "marginal estimate")
  expect_brma_samples_matrix(estimate_fitted, n, "fitted estimate")
  expect_brma_samples_matrix(response, n, "selected response")
  expect_brma_samples_matrix(response_new, nrow(newdata), "new selected response")
  expect_identical(nrow(response), nrow(draws))
  expect_brma_samples_matrix(
    ranef(fit_bselmodel_mv, component = "total", expand = TRUE),
    n,
    "fitted random effects"
  )
  expect_length(residuals(fit_bselmodel_mv), n)
  expect_equal(nrow(rstudent(fit_bselmodel_mv)), n)
})


test_that("exact random-effect posterior draws use the shared covariance plan", {

  posterior_samples <- .get_posterior_samples(fit_bselmodel_mv[["fit"]])
  posterior_samples <- posterior_samples[seq_len(20L), , drop = FALSE]
  S                 <- nrow(posterior_samples)
  K                 <- nobs(fit_bselmodel_mv)
  mu_samples        <- matrix(seq_len(S * K) / 1000, nrow = S, ncol = K)
  bias_offset       <- matrix(seq_len(S * K) / 5000, nrow = S, ncol = K)
  prior_random      <- matrix(seq_len(S * K) / 7000, nrow = S, ncol = K)
  prior_sampling    <- matrix(seq_len(S * K) / 9000, nrow = S, ncol = K)

  testthat::local_mocked_bindings(
    .predict_brma_mv_marginal_random_draws = function(...) prior_random,
    .known_v_sampling_noise = function(...) prior_sampling,
    .package = "RoBMA"
  )
  actual <- .predict_brma_mv_marginal_random_posterior_draws(
    object            = fit_bselmodel_mv,
    mu_samples        = mu_samples,
    posterior_samples = posterior_samples,
    bias_offset       = bias_offset
  )

  random_vcov <- .brma_mv_random_effects_marginal_vcov(
    object            = fit_bselmodel_mv,
    posterior_samples = posterior_samples,
    data              = fit_bselmodel_mv[["data"]]
  )[["samples"]]
  sampling_covariance <- .known_v_covariance_matrix(
    .data_known_v_data(fit_bselmodel_mv[["data"]])
  )
  residual <- matrix(
    fit_bselmodel_mv[["data"]][["outcome"]][["yi"]],
    nrow  = S,
    ncol  = K,
    byrow = TRUE
  ) - mu_samples - bias_offset - prior_random - prior_sampling
  expected <- prior_random
  for (draw in seq_len(S)) {
    Q <- matrix(random_vcov[draw, , ], nrow = K, ncol = K)
    expected[draw, ] <- prior_random[draw, ] +
      as.vector(Q %*% solve(Q + sampling_covariance, residual[draw, ]))
  }

  expect_equal(actual, expected, tolerance = 1e-10)
})


test_that("random-formula qCMDE uses the exact selection likelihood", {

  expect_s3_class(
    plot(
      fit_bselmodel_mv,
      "sd",
      component       = "random",
      density_method  = "qCMDE",
      density_control = list(
        n_points             = 20L,
        samples              = 80L,
        normalization_points = 20L
      ),
      plot_type = "ggplot"
    ),
    "ggplot"
  )
})
