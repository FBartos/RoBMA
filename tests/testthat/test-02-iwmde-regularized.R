source(testthat::test_path("common-functions.R"))

.regularized_regression_data <- function() {
  dat <- expand.grid(x = c(10, 14), z = c(80, 120), replicate = 1:4)
  dat$se <- 0.55
  dat$yi <- 0.12 * as.numeric(scale(dat$x)) + 0.68 * as.numeric(scale(dat$z)) +
    c(-.23, -.08, .08, .23)[dat$replicate]
  dat
}

.fit_regularized_regression <- function(dat, group, sample) {
  brma(yi = yi, sei = se, mods = ~ 0 + x + z, data = dat, measure = "SMD",
       prior_mods = list(shrinkage = group), prior_heterogeneity = prior("point", list(.15)),
       standardize_continuous_predictors = TRUE, chains = 2, adapt = 500,
       burnin = 500, sample = sample, seed = 154, parallel = FALSE, silent = TRUE)
}

test_that("fitted standardized ridge hypotheses and prior overlay use the shared scale", {
  skip_if_not_certification("The fitted global-scale quadrature comparison runs in iwmde-qcmde.")
  dat <- .regularized_regression_data()
  lower <- .025
  upper <- 1
  fit <- .fit_regularized_regression(dat, prior_regularized("ridge", prior("uniform", list(lower, upper))),
                                     sample = 3000)

  # For this balanced orthogonal design the two least-squares slopes are
  # independent given their true coefficients and have the same known variance.
  # Integrate both slopes out analytically to weight the global scale by the
  # likelihood of both sufficient statistics. Only a scalar integral remains.
  X <- scale(as.matrix(dat[c("x", "z")]))
  information <- sum(X[, 1L]^2)
  b <- drop(crossprod(X, dat$yi)) / information
  v <- (dat$se[[1L]]^2 + .15^2) / information
  scale_likelihood <- function(s) {
    stats::dnorm(b[[1L]], 0, sqrt(s^2 + v)) * stats::dnorm(b[[2L]], 0, sqrt(s^2 + v))
  }
  joint_at_zero <- function(s) {
    variance <- 1 / (1 / s^2 + 1 / v)
    scale_likelihood(s) * stats::dnorm(0, variance * b[[1L]] / v, sqrt(variance))
  }
  integrate_scale <- function(f) stats::integrate(f, lower, upper, rel.tol = 1e-10)$value
  posterior_at_zero <- integrate_scale(joint_at_zero) / integrate_scale(scale_likelihood)
  prior_at_zero <- log(upper / lower) / ((upper - lower) * sqrt(2 * pi))
  sd_x <- stats::sd(dat$x)
  for (method in c("qCMDE", "IWMDE")) {
    result <- hypothesis(fit, "x = 0", standardized_coefficients = TRUE,
                         density_method = method, columns = "all", seed = 23)
    row <- as.data.frame(result)
    expect_equal(row$prior, prior_at_zero, tolerance = 1e-6, info = method)
    log_error <- abs(log(row$posterior / posterior_at_zero))
    log_mcse <- sqrt(log1p((row$BF_error / 100)^2))
    expect_lte(log_error, 4 * log_mcse, label = method)
    expect_lt(log_error, .10, label = method)
    expect_equal(row$BF10, row$prior / row$posterior, info = method)
    # On the default reporting scale, beta_x = beta_standardized / sd(x).
    # The same Jacobian multiplies both ordinates, cancelling in the BF.
    original <- as.data.frame(hypothesis(fit, "x = 0", density_method = method,
                                         columns = "all", seed = 23))
    expect_equal(original$prior, sd_x * prior_at_zero, tolerance = 1e-6, info = method)
    expect_equal(original$posterior, sd_x * row$posterior, tolerance = 1e-5, info = method)
    expect_equal(original$BF10, row$BF10, tolerance = 1e-5, info = method)
  }

  expect_prior_overlay <- function(p, jacobian) {
    layers <- ggplot2::ggplot_build(p)$data
    expect_true(any(vapply(layers, function(layer) {
      if (length(layer$x) < 3L || length(layer$y) != length(layer$x)) return(FALSE)
      indices <- vapply(c(0, .3, .8) / jacobian, function(target) which.min(abs(layer$x - target)), integer(1))
      expected <- vapply(layer$x[indices], function(value) {
        jacobian * integrate_scale(function(s) stats::dnorm(jacobian * value, 0, s)) / (upper - lower)
      }, numeric(1))
      isTRUE(all.equal(layer$y[indices], expected, tolerance = 1e-5, check.attributes = FALSE))
    }, logical(1))))
  }
  expect_prior_overlay(plot(fit, parameter = "x", component = "mods", standardized_coefficients = TRUE,
                            prior = TRUE, plot_type = "ggplot", xlim = c(-1, 1)), 1)
  expect_prior_overlay(plot(fit, parameter = "x", component = "mods", prior = TRUE,
                            plot_type = "ggplot", xlim = c(-1, 1) / sd_x), sd_x)
})

test_that("fitted horseshoe point tests at zero refuse before constructing a density grid", {
  skip_if_not_certification("The fitted horseshoe point-null check runs in iwmde-qcmde.")
  fit <- .fit_regularized_regression(.regularized_regression_data(),
    prior_regularized("horseshoe", .3, slab_scale = 1), sample = 300)
  testthat::local_mocked_bindings(
    .prior_linear_combination_grid = function(...) stop("Unexpected prior grid for a singular point null."),
    .package = "BayesTools")
  for (method in c("qCMDE", "IWMDE")) {
    expect_error(hypothesis(fit, "x = 0", standardized_coefficients = TRUE, density_method = method),
                 class = "BayesTools_infinite_ordinate")
  }
})
