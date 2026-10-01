.weighted_diagnostic_fixture <- function(multilevel = FALSE) {

  yi <- c(-0.5, 0.3, 0.2, 0.8, -0.1, 0.4)
  vi <- c(0.04, 0.09, 0.16, 0.01, 0.25, 0.06)
  weights <- c(1, 4, 0.5, 2, 1.5, 0.75)
  object <- brma.norm(yi = yi, vi = vi, weights = weights,
    cluster = if (multilevel) rep(1:3, each = 2L) else NULL,
    measure = "GEN", only_data = TRUE)
  object$priors <- list(outcome = list(tau = NULL, rho = NULL, bias = NULL), scale = NULL)
  samples <- if (multilevel) cbind(tau = 0.5, rho = 0.64) else cbind(tau = 0.3)
  object$fit <- coda::mcmc.list(coda::mcmc(samples))
  class(object) <- c("brma.norm", "brma")
  object
}

test_that("weighted residual scaling separates fit weights from outcome covariance", {

  for (multilevel in c(FALSE, TRUE)) {
    object <- .weighted_diagnostic_fixture(multilevel)
    yi <- .outcome_data_yi(object)
    vi <- .outcome_data_vi(object)
    weights <- .outcome_data_weights(object)
    between <- if (multilevel) outer(rep(1:3, each = 2L), rep(1:3, each = 2L), `==`) * 0.4^2 else matrix(0, 6, 6)
    covariance <- diag(vi + 0.3^2) + between
    fit_covariance <- diag((vi + 0.3^2) / weights) + between
    X <- matrix(1, 6L, 1L)
    precision <- solve(fit_covariance)
    H <- X %*% solve(t(X) %*% precision %*% X) %*% t(X) %*% precision
    transform <- diag(6) - H
    result <- .compute_hat_matrix_samples(object, return_se = TRUE,
      return_resid = TRUE, return_full_H = TRUE)
    expect_equal(as.vector(result$H), as.vector(H), tolerance = 1e-13)
    expect_equal(as.vector(result$resid), as.vector(transform %*% yi), tolerance = 1e-13)
    expect_equal(as.vector(result$se), sqrt(diag(transform %*% covariance %*% t(transform))), tolerance = 1e-13)
    expect_equal(as.vector(result$M_diag), diag(covariance), tolerance = 1e-13)
    expect_equal(as.vector(.pearson_residual_se_samples(object, "marginal")), sqrt(diag(covariance)), tolerance = 1e-13)
    expect_equal(as.vector(.pearson_residual_se_samples(object, "cluster")), sqrt(vi + 0.3^2), tolerance = 1e-13)
    expect_equal(as.vector(.pearson_residual_se_samples(object, "estimate")), sqrt(vi), tolerance = 1e-13)
  }
})

test_that("weighted marginal standardization agrees with metafor custom weights", {

  skip_if_not_installed("metafor")
  for (multilevel in c(FALSE, TRUE)) {
    object <- .weighted_diagnostic_fixture(multilevel)
    yi <- .outcome_data_yi(object)
    vi <- .outcome_data_vi(object)
    weights <- .outcome_data_weights(object)
    reference <- if (multilevel) {
      study <- factor(rep(1:3, each = 2L))
      between <- outer(study, study, `==`) * 0.4^2
      metafor::rma.mv(yi, V = diag(vi + 0.3^2), random = ~ 1 | study,
        sigma2 = 0.4^2, W = solve(diag((vi + 0.3^2) / weights) + between))
    } else {
      metafor::rma.uni(yi, vi, tau2 = 0.3^2, weights = weights / (vi + 0.3^2))
    }
    actual <- rstandard(object)
    expected <- rstandard(reference)
    expect_equal(actual$resid, expected$resid, tolerance = 1e-12)
    expect_equal(actual$se, expected$se, tolerance = 1e-12)
    expect_equal(actual$z, expected$z, tolerance = 1e-12)
  }
})

test_that("weighted LOO-PIT and rstudent moments use the original outcome law", {

  object <- .weighted_diagnostic_fixture()
  yi <- .outcome_data_yi(object)
  sei <- .outcome_data_sei(object)
  setup <- list(data = object$data, priors = object$priors, yi = yi, sei = sei,
    S = 1L, K = 6L, mu = matrix(0.15, 1L, 6L), tau_within = matrix(0.3, 1L, 6L),
    weights = .outcome_data_weights(object), outcome_type = "norm", is_weightfunction = FALSE)
  expected_sd <- sqrt(sei^2 + 0.3^2)
  tails <- .loo_predictive_log_tails_estimate(setup)
  expect_equal(as.vector(tails$log_lower), pnorm(yi, 0.15, expected_sd, log.p = TRUE))
  expect_equal(as.vector(tails$log_upper), pnorm(yi, 0.15, expected_sd, log.p = TRUE, lower.tail = FALSE))
  moments <- .loo_predictive_moments_estimate(object, setup, NULL, matrix(1, 1L, 6L))
  expect_equal(moments$resid, yi - 0.15)
  expect_equal(moments$se, expected_sd)
})
