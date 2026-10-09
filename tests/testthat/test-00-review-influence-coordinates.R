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
    object, transform_factors = FALSE
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

