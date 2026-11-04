context("Influence diagnostics use the fitted coordinates of factor levels")

.review_influence_mixture <- function(scale = FALSE) {

  dat <- data.frame(yi = seq(.1, .6, .1), sei = 1, x = rep(c(-1, 0, 1), 2))
  normal <- prior("normal", list(0, 1))
  point <- prior("point", list(0))
  args <- list(yi = dat$yi, sei = dat$sei, data = dat, prior_effect = normal,
    prior_effect_null = FALSE, prior_heterogeneity_null = NULL, only_priors = TRUE,
    silent = TRUE, measure = "GEN", prior_unit_information_sd = 1)
  if (scale) {
    args <- c(args, list(scale = ~ x, prior_scale = list(x = normal),
      prior_scale_null = list(x = point),
      prior_heterogeneity = prior("normal", list(0, 1), truncation = list(0, Inf))))
  } else {
    args <- c(args, list(mods = ~ x, prior_mods = list(x = normal),
      prior_mods_null = list(x = point), prior_heterogeneity = point))
  }
  setup <- review_compile_object(do.call(BMA.norm, args))
  priors <- c(do.call(c, unname(lapply(setup$compiled, `[[`, "prior_list"))),
    .create_fit_priors(setup$object$data, setup$object$priors))
  primitive <- matrix(numeric(), 24L, 0L)
  for (name in names(priors)) {
    mixture <- priors[[name]]
    stopifnot(BayesTools::is.prior.mixture(mixture))
    slope <- length(mixture) == 2L
    values <- matrix(if (slope) rep(c(1, 1, 2, 2), 6) else 1, 24L, 1L,
      dimnames = list(NULL, paste0(name, "_indicator")))
    primitive <- cbind(primitive, values)
    for (j in seq_along(mixture)) {
      value <- if (BayesTools::is.prior.point(mixture[[j]])) mixture[[j]]$parameters$location else {
        if (slope) rep(c(.25, .5, 1, 2), 6) else rep(rep(c(-1, 1), each = 4), 3) + if (scale) 2 else 0
      }
      primitive <- cbind(primitive, matrix(value, 24L, 1L,
        dimnames = list(NULL, paste0(name, "_component_", j))))
    }
  }
  seed <- review_attach_replay(setup, primitive)
  generated <- BayesTools::JAGS_evaluate_deterministic(seed$fit, draws = primitive)
  pool <- cbind(primitive, generated[, setdiff(colnames(generated), colnames(primitive)), drop = FALSE])
  monitors <- BayesTools::JAGS_to_monitor(priors)
  names <- unique(unlist(lapply(monitors, function(name) {
    if (name %in% colnames(pool)) name else colnames(pool)[startsWith(colnames(pool), paste0(name, "["))]
  })))
  review_attach_replay(setup, pool[, names, drop = FALSE])
}

test_that("mixture influence coordinates exclude declared inclusion indicators", {

  for (scale in c(FALSE, TRUE)) {
    object <- .review_influence_mixture(scale)
    component <- if (scale) "scale" else "mods"
    weights <- matrix(rep(c(3/8, 3/8, 1/8, 1/8)/6, 6), 24L, 6L)
    coefficients <- BayesTools::JAGS_estimates_table(object$fit,
      keep_formulas = if (scale) "log_tau" else "mu", remove_inclusion = TRUE,
      remove_diagnostics = TRUE, transform_factors = TRUE, transform_scaled = TRUE,
      return_samples = TRUE)
    result <- covratio(object, component = component, .weights = weights)
    expect_equal(as.numeric(result), rep(31/44, 6), tolerance = 1e-12)
    expect_identical(colnames(dfbetas(object, component = component, .weights = weights)),
      colnames(coefficients))
    expect_equal(as.numeric(covratio(object, component = component,
      .weights = matrix(1/24, 24L, 6L))), rep(1, 6), tolerance = 1e-12)
  }
})

test_that("COVRATIO of a dependent factor-level parameterization uses the fitted coordinates", {

  # The levels of a mean-difference factor are linearly dependent (they sum to
  # zero), so their covariance is singular; the fitted coordinates are not.
  n      <- 100L
  object <- review_factor_replay(c("A", "B"), "meandif", independent = TRUE, n = n)
  weights <- matrix(1 / n, nrow = n, ncol = 6L)

  levels <- .diagnostic_location_parameter_samples(object)
  expect_equal(qr(scale(levels, center = TRUE, scale = FALSE))[["rank"]], 2L)
  expect_equal(ncol(levels), 3L)

  observed <- suppressWarnings(covratio(object, .weights = weights))
  coordinates <- .diagnostic_location_parameter_samples(
    object, standardized_coefficients = TRUE, transform_factors = FALSE
  )
  expected <- .covratio_internal(coordinates, weights)

  expect_true(all(is.finite(unclass(observed))))
  expect_equal(as.numeric(observed), expected[["values"]], tolerance = 1e-12)
  # equal weights leave the posterior covariance unchanged
  expect_equal(as.numeric(observed), rep(1, 6L), tolerance = 1e-10)
  expect_null(attr(observed, "note"))
  expect_s3_class(observed, "covratio.brma")
  expect_identical(names(as.data.frame(observed)), c("component", "parameter", "value"))
  expect_identical(as.data.frame(observed), data.frame(observed))
  expect_identical(as.numeric(as.data.frame(observed)$value), as.numeric(observed))
  expect_output(print(observed), "Study 1", fixed = TRUE)
})


# A replay of a `~ 0 + g` model: the location has no intercept, so the
# structural-zero intercept is not a parameter of the posterior. The level
# draws scatter around the level means of the data.
.review_zero_intercept_factor_replay <- function(n = 100L) {

  levels <- c("A", "B", "C")
  dat    <- data.frame(
    yi  = seq(-.1, .2, length.out = 6),
    sei = .1,
    g   = factor(rep(levels, each = 2L), levels = levels)
  )
  object <- brma(
    yi = yi, sei = sei, mods = ~ 0 + g, data = dat, measure = "GEN",
    prior_unit_information_sd = 1,
    prior_mods = list(g = BayesTools::prior_factor(
      "normal", list(mean = 0, sd = 1), contrast = "independent"
    )),
    prior_heterogeneity = BayesTools::prior("point", list(location = .1)),
    only_priors = TRUE, silent = TRUE
  )
  setup <- review_compile_object(object)
  withr::local_seed(6123)
  means  <- tapply(dat[["yi"]], dat[["g"]], mean)
  slopes <- vapply(means, function(m) m + .07 * stats::rnorm(n), numeric(n))
  colnames(slopes) <- BayesTools:::.JAGS_prior_factor_names(
    review_coefficient_name(setup, "g"),
    setup[["compiled"]][["mu"]][["prior_list"]][[review_coefficient_name(setup, "g")]]
  )
  samples <- cbind(slopes, intercept = rep(0, n), tau = rep(.1, n))
  colnames(samples)[4L] <- review_coefficient_name(setup, "intercept")

  review_attach_replay(setup, samples)
}


# Normalized weights that differ between the estimates.
.review_tilted_weights <- function(object, K = 6L) {

  draws <- as.numeric(as.matrix(coda::as.mcmc.list(object[["fit"]]))[, 1L])
  out   <- vapply(seq_len(K), function(i) {
    w <- stats::dnorm(draws, mean = .1 * i - .3, sd = 1)
    w / sum(w)
  }, numeric(length(draws)))

  out
}


test_that("the influence diagnostics omit a fixed-zero intercept like the summary", {

  object <- .review_zero_intercept_factor_replay()
  levels <- c("(mu) g[A]", "(mu) g[B]", "(mu) g[C]")

  # the summary does not report the structural zero; neither do the diagnostics
  expect_true(.location_omit_fixed_zero_intercept(object, transform_scaled = TRUE))
  expect_identical(colnames(.diagnostic_location_parameter_samples(object)), levels)
  expect_identical(
    colnames(.diagnostic_location_parameter_samples(
      object, transform_factors = FALSE
    )),
    levels
  )
  expect_identical(
    colnames(.diagnostic_location_parameter_samples(
      object, standardized_coefficients = TRUE
    )),
    levels
  )

  observed <- suppressWarnings(dfbetas(object, .weights = .review_tilted_weights(object)))
  expect_identical(colnames(observed), levels)
  expect_false(anyNA(as.matrix(observed)))
  expect_null(attr(observed, "note"))

  object   <- suppressWarnings(add_loo(object, unit = "estimate"))
  observed <- suppressWarnings(influence(object))
  expect_identical(colnames(observed[["dfbs"]]), levels)
  expect_false(anyNA(as.matrix(observed[["dfbs"]])))
  expect_true(all(is.finite(observed[["inf"]][["cov.r"]])))
  expect_null(attr(observed, "note"))
})


test_that("COVRATIO of a standardized model without intercept uses the fitted coordinates", {

  # `~ 0 + x` with the standardized predictor x = 1:6: on the original scale
  # the intercept is a linear function of the slope (-mean(x) / sd(x) times
  # it), so the original-scale covariance is singular although the fitted
  # slope has a valid variance
  dat    <- data.frame(yi = seq(-.1, .2, length.out = 6), sei = .1, x = 1:6)
  object <- brma(
    yi = yi, sei = sei, mods = ~ 0 + x, data = dat, measure = "GEN",
    prior_unit_information_sd = 1,
    prior_mods = list(x = BayesTools::prior("normal", list(mean = 0, sd = 1))),
    prior_heterogeneity = BayesTools::prior("point", list(location = .1)),
    standardize_continuous_predictors = TRUE,
    only_priors = TRUE, silent = TRUE
  )
  setup <- review_compile_object(object)
  n     <- 100L
  slope <- stats::qnorm((seq_len(n) - .5) / n) + .3
  samples <- cbind(slope, 0, .1)
  colnames(samples) <- c(review_coefficient_name(setup, "x"),
                         review_coefficient_name(setup, "intercept"), "tau")
  object <- review_attach_replay(setup, samples)

  # the original-scale coordinates are linearly dependent
  original <- .diagnostic_location_parameter_samples(
    object, transform_factors = FALSE
  )
  expect_identical(colnames(original), c("(mu) x", "(mu) intercept"))
  expect_equal(qr(scale(original, center = TRUE, scale = FALSE))[["rank"]], 1L)

  # with uniform weights every deletion leaves the posterior unchanged
  weights  <- matrix(1 / n, nrow = n, ncol = 6L)
  observed <- covratio(object, .weights = weights)
  expect_equal(as.numeric(observed), rep(1, 6L), tolerance = 1e-12)
  expect_null(attr(observed, "note"))

  # one slope parameter: the ratio of the weighted (ML) variance of the slope
  # to its equally weighted variance
  weights[, 1L] <- stats::dnorm(slope, mean = .5, sd = .7)
  weights[, 1L] <- weights[, 1L] / sum(weights[, 1L])
  weighted_mean <- sum(weights[, 1L] * slope)
  reference     <- sum(weights[, 1L] * (slope - weighted_mean)^2) /
    mean((slope - mean(slope))^2)
  observed <- covratio(object, .weights = weights)

  # independent value: 0.333106352763284 (mpmath, 40 digits)
  expect_equal(reference, 0.333106352763284, tolerance = 1e-12)
  expect_equal(as.numeric(observed)[1L], reference, tolerance = 1e-12)
  expect_equal(as.numeric(observed)[-1L], rep(1, 5L), tolerance = 1e-12)
  expect_null(attr(observed, "note"))
})


test_that("VIF omits a fixed-zero intercept from the posterior correlation", {

  object <- .review_zero_intercept_factor_replay()
  levels <- c("(mu) g[A]", "(mu) g[B]", "(mu) g[C]")

  warnings <- character()
  result   <- withCallingHandlers(
    vif(object),
    warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(warnings, character())

  # the correlation of the three independent level draws
  correlation <- result[["posterior_correlation"]]
  expect_identical(rownames(correlation), levels)
  expect_identical(colnames(correlation), levels)
  expect_false(anyNA(correlation))
  draws <- as.matrix(coda::as.mcmc.list(object[["fit"]]))[, 1:3]
  expect_equal(unname(correlation), unname(stats::cor(draws)), tolerance = 1e-12)
})


test_that("scale DFBETAS follow standardized_coefficients like location DFBETAS", {

  n    <- 100L
  data <- data.frame(yi = seq(-.1, .2, length.out = 6), sei = .1, x = 1:6)
  object <- brma(
    yi = yi, sei = sei, mods = ~ x, scale = ~ x, data = data, measure = "GEN",
    prior_unit_information_sd = 1,
    prior_effect = BayesTools::prior("normal", list(mean = 0, sd = 1)),
    prior_mods = list(x = BayesTools::prior("normal", list(mean = 0, sd = 1))),
    prior_scale = list(x = BayesTools::prior("normal", list(mean = 0, sd = 1))),
    prior_heterogeneity = BayesTools::prior(
      "normal", list(mean = 0, sd = 1), truncation = list(lower = 0, upper = Inf)
    ),
    standardize_continuous_predictors = TRUE, only_priors = TRUE, silent = TRUE
  )
  setup <- review_compile_object(object)
  withr::local_seed(6123)
  samples <- cbind(
    stats::rnorm(n) + .2, stats::rnorm(n) + .3,
    exp(stats::rnorm(n, -2, .2)), stats::rnorm(n) + .4
  )
  colnames(samples) <- c(
    review_coefficient_name(setup, "intercept"), review_coefficient_name(setup, "x"),
    review_coefficient_name(setup, "intercept", "log_tau"),
    review_coefficient_name(setup, "x", "log_tau")
  )
  object  <- review_attach_replay(setup, samples)
  weights <- matrix(stats::runif(n * 6L), nrow = n, ncol = 6L)
  weights <- sweep(weights, 2L, colSums(weights), "/")

  # leave-one-out means of the scale coefficients: the fitted intercept b0 and
  # slope b1 of tau = b0 exp(b1 z), z = (x - m) / s, give tau0 = b0 exp(-b1 m / s)
  # and b1 / s on the original scale
  center  <- mean(data[["x"]])
  scale   <- stats::sd(data[["x"]])
  b0      <- samples[, 3L]
  b1      <- samples[, 4L]
  loo <- function(draws) unname(crossprod(weights, draws))

  standardized <- dfbetas(object, type = "scale", .weights = weights,
                          return_loo_estimates = TRUE,
                          standardized_coefficients = TRUE)
  original     <- dfbetas(object, type = "scale", .weights = weights,
                          return_loo_estimates = TRUE,
                          standardized_coefficients = FALSE)

  expect_equal(unname(as.matrix(standardized)), loo(cbind(b0, b1)), tolerance = 1e-10)
  expect_equal(unname(as.matrix(original)),
               loo(cbind(b0 * exp(-b1 * center / scale), b1 / scale)),
               tolerance = 1e-10)
  expect_false(isTRUE(all.equal(standardized, original, check.attributes = FALSE)))
  # the default is the original scale, as for the location coefficients
  expect_identical(
    dfbetas(object, type = "scale", .weights = weights, return_loo_estimates = TRUE),
    original
  )
})
