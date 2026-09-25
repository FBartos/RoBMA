context("Public likelihood-aware density integration")

source(testthat::test_path("common-functions.R"))


test_that("public density hypotheses and plots match conjugate Gaussian inference", {

  skip_on_cran()
  name <- "brma.mv_block_mvn_fixed_random_null"
  skip_if_missing_fits(name)
  fit <- load_fit(name)
  info <- load_info(name)
  prior <- fit[["priors"]][["outcome"]][["mu"]]
  expect_identical(prior[["distribution"]], "normal")
  expect_identical(prior[["truncation"]], list(lower = -Inf, upper = Inf))

  # The original input covariance and the declared Normal prior give the
  # posterior in closed form; no package likelihood/density helper is used.
  V <- info[["V"]]
  yi <- info[["data"]][["yi"]]
  prior_mean <- prior[["parameters"]][["mean"]]
  prior_sd <- prior[["parameters"]][["sd"]]
  posterior_variance <- 1 / (1 / prior_sd^2 + sum(solve(V, rep(1, length(yi)))))
  posterior_mean <- posterior_variance * (
    prior_mean / prior_sd^2 + sum(solve(V, yi))
  )
  posterior_sd <- sqrt(posterior_variance)
  oracle_bf <- stats::dnorm(0, prior_mean, prior_sd) /
    stats::dnorm(0, posterior_mean, posterior_sd)

  for (method in c("qCMDE", "IWMDE")) {
    control <- list(n_points = 40L, samples = 240L)
    result <- hypothesis(fit, "mu = 0", density_method = method, density_control = control)
    diagnostics <- density_diagnostics(result)
    expect_identical(diagnostics[["density_method"]], method)
    expect_true(all(diagnostics[["bf_grade_met"]]))
    expect_true(all(diagnostics[["evaluated_rows"]] > 0L))
    expect_true(all(is.finite(result[["BF_error"]])))

    # Allow four reported Monte Carlo SEs plus the deterministic normalization
    # error. Explicit log/peak-scaled errors have the same meaning in every
    # testthat edition and do not become loose absolute density tolerances.
    relative_mcse <- result[["BF_error"]][[1L]] / 100
    monte_carlo <- 4 * sqrt(log1p(relative_mcse^2))
    if (identical(method, "qCMDE")) {
      # Error model. Every row's conditional law is this Normal posterior (the
      # Normal prior's density is exact), so each row is an "exact" row. The
      # normalization range covers its central 'normalization_prob' (default
      # 0.999) interval, leaving it the reported truncation t <= 1 - 0.999.
      # A normalizer missing the fraction t overstates the density by the
      # factor 1 / (1 - t) = 1 + b with the reported bound b = t / (1 - t),
      # and the trapezoid rule adds a relative error d with |d| at most the
      # reported change c between the nested grids, which overstates the
      # selected grid's error. The posterior ordinate is f (1 + b)(1 + d), so
      # |log BF error| <= log(1 + b) - log(1 - c). Rows sharing one law share
      # one normalizer, so the whole curve carries the same relative error,
      # at most (1 + b)(1 + c) - 1 after peak scaling, and the plot's line
      # with the same control uses the same range and grids.
      bound  <- diagnostics[["truncation_ordinate_bound"]]
      change <- diagnostics[["ordinate_relative_change"]]
      expect_identical(diagnostics[["normalization_truncation_status"]], "exact")
      expect_lte(diagnostics[["normalization_truncation"]], (1 - .999) * (1 + 1e-9))
      expect_equal(bound, diagnostics[["normalization_truncation"]] /
        (1 - diagnostics[["normalization_truncation"]]))
      allowance <- log1p(bound) - log1p(-change) + monte_carlo
      curve_allowance <- (1 + bound) * (1 + change) - 1 + monte_carlo
    } else {
      # IWMDE keeps its 0.1% deterministic normalization-mass allowance.
      allowance <- curve_allowance <- 1e-3 + monte_carlo
    }
    expect_lte(abs(log(attr(result, "raw_BF")[[1L]]) - log(oracle_bf)), allowance)

    plot <- plot(fit, parameter = "mu", density_method = method,
      density_control = control, plot_type = "ggplot")
    expect_s3_class(plot, "ggplot")
    line <- which(vapply(plot[["layers"]], function(layer) {
      inherits(layer[["geom"]], "GeomLine")
    }, logical(1)))
    expect_length(line, 1L)
    curve <- plot[["layers"]][[line]][["data"]]
    oracle <- stats::dnorm(curve[["x"]], posterior_mean, posterior_sd)
    expect_true(all(is.finite(curve[["y"]])))
    expect_lte(max(abs(curve[["y"]] - oracle)) / max(oracle), curve_allowance)
  }
})


test_that("marginal means retain genuine qCMDE curves and usable point ordinates", {

  skip_on_cran()
  name <- "bcg_meta-regression2"
  skip_if_missing_fits(name)
  fit <- load_fit(name)
  # Retain the public default row budget. These calls run the genuine estimator
  # and its acceptance gates; the attachment tests elsewhere mock only routing.
  means <- marginal_means(fit, parameter = "alloc", bf = TRUE,
    n_samples = 250L, density_method = "qCMDE", density_control = list(n_points = 40L))
  cells <- means[["inference"]][["conditional"]][["mu_alloc"]]
  expect_named(cells, c("alternate", "random", "systematic"))
  for (level in names(cells)) {
    curve <- BayesTools::posterior_metadata(cells[[level]], "posterior_density")
    ordinate <- BayesTools::posterior_metadata(cells[[level]], "posterior_ordinate")
    expect_identical(curve[["density_method"]], "qCMDE", info = level)
    expect_identical(curve[["status"]], "ok", info = level)
    expect_length(curve[["x"]], 40L)
    expect_true(all(is.finite(curve[["y"]]) & curve[["y"]] >= 0), info = level)
    expect_identical(ordinate[["status"]], "ok", info = level)
    expect_true(is.finite(ordinate[["ordinate"]]) && ordinate[["ordinate"]] > 0)
    expect_true(ordinate[["diagnostics"]][["n_evaluated_rows"]] > 0L)
  }

  # Under the declared independent Normal treatment priors the random-cell
  # prior is the sum of the intercept and its one treatment coefficient.
  priors <- attr(fit[["fit"]], "prior_list", exact = TRUE)
  intercept <- priors[["mu_intercept"]]
  treatment <- priors[["mu_alloc"]]
  expect_identical(intercept[["distribution"]], "normal")
  expect_identical(treatment[["distribution"]], "normal")
  expect_true(inherits(treatment, "prior.treatment"))
  prior_mean <- intercept[["parameters"]][["mean"]] + treatment[["parameters"]][["mean"]]
  prior_sd <- sqrt(intercept[["parameters"]][["sd"]]^2 + treatment[["parameters"]][["sd"]]^2)
  ordinate <- BayesTools::posterior_metadata(cells[["random"]], "posterior_ordinate")
  expected_bf <- stats::dnorm(0, prior_mean, prior_sd) / ordinate[["ordinate"]]
  result <- hypothesis(means, "alloc[random] = 0", columns = "all")
  expect_lt(abs(log(attr(result, "raw_BF")[[1L]] / expected_bf)), 1e-3)
  expect_equal(result[["BF_error"]][[1L]], ordinate[["diagnostics"]][["BF_error_percent"]])
})
