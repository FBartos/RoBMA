source(testthat::test_path("common-functions.R"))

test_that("mixed-contrast interaction summaries use the fitted cell design and one term gate", {
  skip_if_not_certification("The fitted interaction summary runs in iwmde-qcmde.")
  dat <- expand.grid(f = factor(c("A", "B", "C")), h = factor(c("U", "V", "W")), r = 1:2)
  dat$yi <- .12 * as.integer(dat$f) - .05 * as.integer(dat$h) + .1 * sin(seq_len(nrow(dat)))
  dat$se <- .3
  fit <- BMA(yi = yi, sei = se, data = dat, mods = ~ f * h, measure = "SMD",
             prior_mods = list(f = prior_factor("mnormal", list(0, .3), contrast = "orthonormal"),
                               h = prior_factor("normal", list(0, .3), contrast = "treatment"),
                               pool = prior_regularized("lasso", terms = "f:h", scale = .12)),
             prior_effect_null = NULL, prior_heterogeneity = prior("point", list(.2)),
             prior_heterogeneity_null = NULL, chains = 2, adapt = 300, burnin = 500,
             sample = 1000, parallel = FALSE, seed = 141, silent = TRUE)
  raw <- as.matrix(coda::as.mcmc.list(fit$fit))
  columns <- paste0("mu_f__xXx__h[", 1:4, "]")
  inclusion <- raw[, "mu_f__xXx__h_indicator"] == 2L
  expect_true(any(inclusion) && any(!inclusion))
  expect_true(all(raw[!inclusion, columns] == 0))
  design <- kronecker(stats::contr.treatment(3), BayesTools::contr.orthonormal(3))
  values <- raw[, columns] %*% t(design)
  cells <- expand.grid(f = levels(dat$f), h = levels(dat$h))
  info <- summary(fit, conditional = TRUE)
  pattern <- "^f\\[(?:dif: )?([ABC])\\]:h\\[(?:dif: )?([UVW])\\]$"
  for (conditional in c(FALSE, TRUE)) {
    table <- as.data.frame(if (conditional) info$estimates_mods_conditional else info$estimates_mods)
    rows <- grep(pattern, rownames(table), perl = TRUE)
    keys <- sub(pattern, "\\1\\2", rownames(table)[rows], perl = TRUE)
    independent <- paste0(cells$f, cells$h)
    expect_true(all(independent[cells$h != "U"] %in% keys))
    selected <- if (conditional) inclusion else rep(TRUE, nrow(raw))
    means <- colMeans(values[selected, , drop = FALSE])
    expect_equal(unname(table$Mean[rows]), unname(means[match(keys, independent)]), tolerance = 1e-11)
  }
})

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

test_that("centered shrinkage factor hypotheses use the fitted contrasts and term gate", {
  skip_if_not_certification("The fitted factor contrasts run in iwmde-qcmde.")
  dat <- expand.grid(f = factor(c("A", "B", "C")), x = c(-1, 1), replicate = 1:3)
  dat$yi <- c(-0.04, 0.08, -0.04)[dat$f] + 0.11 * dat$x + c(-.13, 0, .13)[dat$replicate]
  global <- prior("uniform", list(0.15, 0.55))
  for (coding in c("meandif", "orthonormal")) {
    C <- if (coding == "meandif") BayesTools::contr.meandif(3) else BayesTools::contr.orthonormal(3)
    for (averaged in c(FALSE, TRUE)) {
      run <- if (averaged) BMA else brma
      group <- if (averaged) prior_regularized("horseshoe", global, slab_scale = 1) else
        prior_regularized("lasso", global, contrast = coding)
      settings <- list(yi = dat$yi, sei = rep(.5, nrow(dat)), mods = ~ f + x, data = dat,
                       measure = "SMD", prior_mods = list(pool = group),
                       prior_effect = prior("normal", list(0, .4)),
                       prior_heterogeneity = prior("point", list(.2)),
                       set_contrast_factor_predictors = if (averaged) coding else "treatment",
                       chains = 2, adapt = 300, burnin = 500, sample = 1500, parallel = FALSE,
                       seed = 2249, silent = TRUE)
      if (averaged) settings[c("prior_effect_null", "prior_heterogeneity_null")] <- list(NULL, NULL)
      fit <- do.call(run, settings)
      sampled <- as.matrix(coda::as.mcmc.list(fit$fit))
      beta <- sampled[, c("mu_f[1]", "mu_f[2]")]
      included <- rep(TRUE, nrow(sampled))
      if (averaged) {
        included <- sampled[, "mu_f_indicator"] == 2L
        expect_true(any(included) && any(!included))
        expect_true(all(beta[!included, , drop = FALSE] == 0))
      }
      delta <- drop(beta[included, , drop = FALSE] %*% (C[2L, ] - C[1L, ]))
      posterior_probability <- mean(delta > 0)
      # Symmetry gives prior odds 1. The empirical posterior odds use the
      # two factor-level effects computed independently from the JAGS columns.
      res <- hypothesis(fit, "f[B] > f[A]", component = "mods", conditional = averaged,
                        columns = "all", seed = 1704)
      row <- as.data.frame(res)
      expect_equal(row$BF10, posterior_probability / (1 - posterior_probability),
                   tolerance = 1e-5)
      if (!averaged) {
        # A projection of the two-dimensional radial lasso has ordinate
        # 1/(pi*s*||C[B,]-C[A,]||) at zero given its global scale s.
        point <- hypothesis(fit, "f[B] = f[A]", component = "mods", columns = "all",
                            density_method = "qCMDE", seed = 1705)
        expected_prior <- log(.55/.15) / (.55 - .15) / (pi * sqrt(sum((C[2L, ] - C[1L, ])^2)))
        point_row <- as.data.frame(point)
        expect_equal(point_row$prior, expected_prior, tolerance = 1e-5)
        expect_true(is.finite(point_row$BF10))
      }
      expect_equal(nrow(posterior::as_draws_matrix(fit)), nrow(sampled))
      expect_false(is.null(summary(fit)))
    }
  }
})
