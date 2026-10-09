context("Influence diagnostics use the fitted coordinates of factor levels")

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
