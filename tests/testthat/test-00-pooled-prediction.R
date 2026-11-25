context("Pooled prediction targets")


test_that("pooled effect prediction intervals are exact predictive mixture quantiles", {

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
  # independent reference: bisection on the generalized inverse of
  # F(x) = mean(pnorm(x, mu, tau))
  reference_quantile <- function(p, location, scale) {
    lower <- -20
    upper <- 20
    for (i in 1:200) {
      middle <- (lower + upper) / 2
      if (mean(stats::pnorm(middle, location, scale)) < p) lower <- middle else upper <- middle
    }
    upper
  }
  reference <- function(probs) {
    vapply(probs, reference_quantile, numeric(1),
           location = posterior_samples[, "mu"], scale = posterior_samples[, "tau"])
  }

  set.seed(481)
  seed_before <- .Random.seed
  kind_before <- RNGkind()
  pooled <- pooled_effect(object, .posterior_samples = posterior_samples)
  expect_identical(.Random.seed, seed_before)
  expect_identical(RNGkind(), kind_before)
  mixture <- attr(pooled, "prediction_mixture", exact = TRUE)
  expect_named(mixture, c("location", "scale"))
  # The mean of equal contributions is exact in extended precision; the gap on
  # Apple silicon is at most 5.6e-17.
  .test_expect_exact_or_within(unname(mixture[["location"]][, 1L]), unname(posterior_samples[, "mu"]),
                               4 * .Machine$double.eps)
  expect_identical(unname(mixture[["scale"]][, 1L]), unname(posterior_samples[, "tau"]))
  expect_identical(colnames(mixture[["location"]]), "mu")
  expect_identical(summary(pooled_effect(object, .posterior_samples = posterior_samples)),
                   summary(pooled))

  estimates <- summary(pooled)
  expect_true(all(c("PI 0.025", "PI 0.975") %in% colnames(estimates)))
  expect_lt(max(abs(unlist(estimates["mu", c("PI 0.025", "PI 0.975")]) -
                    reference(c(0.025, 0.975)))), 1e-9)
  pooled_df <- as.data.frame(estimates)
  expect_identical(
    names(pooled_df),
    c(
      "component", "parameter", "Mean", "Median",
      "CI_0.025", "CI_0.975", "PI_0.025", "PI_0.975"
    )
  )
  expect_identical(pooled_df[["component"]], "location")
  expect_identical(pooled_df[["parameter"]], "mu")
  expect_identical(data.frame(estimates), pooled_df)
  expect_identical(
    names(as.data.frame(summary(pooled, probs = c(.05, .95)))),
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
  custom_df      <- as.data.frame(custom_summary)
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
  expect_identical(data.frame(custom_summary), custom_df)
  expect_equal(
    unname(as.numeric(custom_df[1L, paste0("CI_", custom_probs)])),
    unname(stats::quantile(as.matrix(custom_pooled)[, 1L], custom_probs)),
    tolerance = 0
  )
  expect_lt(max(abs(as.numeric(custom_df[1L, paste0("PI_", custom_probs)]) -
                    reference(custom_probs))), 1e-9)

  # strictly increasing transformations map the quantiles exactly
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
  expect_identical(attr(transformed, "prediction_mixture"), mixture)
  expect_equal(
    unlist(summary(transformed)["mu", c("PI 0.025", "PI 0.975")], use.names = FALSE),
    exp(unlist(estimates["mu", c("PI 0.025", "PI 0.975")], use.names = FALSE)),
    tolerance = 1e-14
  )

  # a mixture of identical components is that normal distribution
  identical_rows <- cbind(mu = rep(0.3, 4L), tau = rep(0.2, 4L))
  expect_lt(max(abs(
    unlist(summary(pooled_effect(object, .posterior_samples = identical_rows))[
      "mu", c("PI 0.025", "PI 0.975")]) - stats::qnorm(c(0.025, 0.975), 0.3, 0.2)
  )), 1e-9)
})


test_that("pooled prediction intervals use no random numbers", {

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }

  object <- brma(
    yi                        = c(-0.2, 0.1, 0.4),
    sei                       = c(0.2, 0.2, 0.2),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  posterior_samples <- cbind(
    mu  = c(-0.2, 0.0, 0.2, 0.4),
    tau = c(0.1, 0.2, 0.3, 0.4)
  )

  first <- summary(pooled_effect(object, .posterior_samples = posterior_samples))
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  set.seed(1)
  second <- summary(pooled_effect(object, .posterior_samples = posterior_samples))
  set.seed(2)
  third <- summary(pooled_effect(object, .posterior_samples = posterior_samples))
  expect_identical(first, second)
  expect_identical(first, third)
})


test_that("pooled prediction quantiles are generalized inverses at point masses", {

  quantile <- function(p, location, scale) {
    .prediction_mixture_quantile(p, location = location, scale = scale)
  }
  location <- c(-0.2, 0.0, 0.2, 0.4)

  # without heterogeneity the predictive distribution is the location draws
  atoms <- rep(0, 4L)
  expect_identical(quantile(0.025, location, atoms), -0.2)
  expect_identical(quantile(0.25, location, atoms), -0.2)
  expect_identical(quantile(0.5, location, atoms), 0.0)
  expect_identical(quantile(0.975, location, atoms), 0.4)
  expect_identical(quantile(0, location, atoms), -0.2)
  expect_identical(quantile(1, location, atoms), 0.4)

  # half the mass at zero: quantiles inside its jump are zero exactly, others
  # solve the continuous part
  mixed_location <- c(0, 0.5, 0, 0.8)
  mixed_scale    <- c(0, 0.2, 0, 0.3)
  below_zero     <- mean(stats::pnorm(0, mixed_location, mixed_scale)) - 0.5
  expect_gt(below_zero, 0)
  expect_lt(below_zero, 0.025)
  expect_identical(quantile(0.025, mixed_location, mixed_scale), 0)
  expect_identical(quantile(0.3, mixed_location, mixed_scale), 0)
  upper <- quantile(0.975, mixed_location, mixed_scale)
  expect_lt(abs(mean(stats::pnorm(upper, mixed_location, mixed_scale)) - 0.975), 1e-10)
  expect_identical(quantile(0, mixed_location, mixed_scale), -Inf)
  expect_identical(quantile(1, mixed_location, mixed_scale), Inf)

  expect_error(
    .new_prediction_mixture(matrix(0, 2L, 1L), matrix(-1, 2L, 1L), matrix(0, 2L, 1L)),
    "non-negative finite scale"
  )

  # many duplicated point masses (null-effect, null-heterogeneity draws)
  duplicated_location <- c(rep(0, 1e4), rep(1, 1e4))
  duplicated_scale    <- c(rep(0, 1e4), rep(0.5, 1e4))
  expect_identical(quantile(0.025, duplicated_location, duplicated_scale), 0)
  expect_identical(quantile(0.4, duplicated_location, duplicated_scale), 0)
})


test_that("pooled prediction quantiles follow the scale of the mixture", {

  quantile <- function(p, location, scale) {
    .prediction_mixture_quantile(p, location = location, scale = scale)
  }
  probs <- c(0.025, 0.5, 0.975)

  # translated and tiny-scale normals: the bracket and tolerance are relative
  translated <- vapply(probs, quantile, numeric(1), location = rep(1e7, 4L), scale = rep(0.1, 4L))
  expect_lt(max(abs(translated - (1e7 + 0.1 * stats::qnorm(probs)))), 1e-8)
  tiny <- vapply(probs, quantile, numeric(1), location = rep(0, 4L), scale = rep(1e-12, 4L))
  expect_lt(max(abs(tiny - 1e-12 * stats::qnorm(probs))), 1e-21)
  expect_true(all(diff(tiny) > 0))

  # equivariance under a + b * x for a mixture with a point mass
  location <- c(-0.2, 0.1, 0.1, 0.6)
  scale    <- c(0.3, 0, 0.1, 0.4)
  base     <- vapply(probs, quantile, numeric(1), location = location, scale = scale)
  moved    <- vapply(probs, quantile, numeric(1), location = 5 + 1e-6 * location, scale = 1e-6 * scale)
  expect_lt(max(abs((moved - 5) / 1e-6 - base)), 1e-6)

  # finite scales near the double range: representable quantiles are found,
  # unrepresentable ones fail clearly
  expect_lt(abs(quantile(0.5, rep(0, 4L), rep(1e308, 4L))) / 1e308, 1e-9)
  expect_error(quantile(0.975, rep(0, 4L), rep(1e308, 4L)), "outside the double-precision range")

  # a quantile in the low-density gap of a separated mixture, against a
  # log-scale reference that does not subtract from a rounded CDF
  separated_location <- c(-10, rep(10, 39L))
  separated_scale    <- rep(1, 40L)
  reference <- stats::uniroot(function(x) {
    stats::pnorm(x + 10, lower.tail = FALSE, log.p = TRUE) - log(39) -
      stats::pnorm(x - 10, log.p = TRUE)
  }, lower = -5, upper = 5, tol = 1e-12)[["root"]]
  expect_lt(abs(quantile(0.025, separated_location, separated_scale) - reference), 1e-6)
  symmetric_location <- rep(c(-10, 10), each = 20L)
  expect_lt(abs(quantile(0.5, symmetric_location, rep(1, 40L))), 1e-8)

  # gaps where every tail underflows: the log tails keep the sign
  expect_lt(abs(quantile(0.5, rep(c(-40, 40), each = 20L), rep(1, 40L))), 1e-8)
  expect_lt(abs(quantile(0.5, rep(c(-0.4, 0.4), each = 20L), rep(0.01, 40L))), 1e-12)
  wide_reference <- stats::uniroot(function(x) {
    stats::pnorm(x + 40, lower.tail = FALSE, log.p = TRUE) - log(39) -
      stats::pnorm(x - 40, log.p = TRUE)
  }, lower = -20, upper = 20, tol = 1e-12)[["root"]]
  expect_lt(abs(quantile(0.025, c(-40, rep(40, 39L)), rep(1, 40L)) - wide_reference), 1e-6)

  # point masses far apart stay distinct; next to a normal whose range merges
  # them in double precision, the quantile fails clearly
  expect_identical(quantile(0.5, c(1, 2, 1e20), c(0, 0, 0)), 2)
  expect_error(quantile(0.5, c(1, 2, 1e20), c(0, 0, 1)), "spacing of its point masses")

  # an affine reconstruction that would overflow on the way
  expect_lt(abs(quantile(0.975, rep(-1e308, 4L), rep(1e308, 4L)) /
                  (1e308 * (stats::qnorm(0.975) - 1)) - 1), 1e-9)
  expect_error(quantile(0.5, c(-1e308, 1e308), c(1e-308, 1e-308)), "too wide a range")

  # many duplicated point masses next to normals stay cheap
  mixed_location <- c(rep(0, 5e4), stats::qnorm(seq(0.01, 0.99, length.out = 5e4)))
  mixed_scale    <- c(rep(0, 5e4), rep(0.5, 5e4))
  expect_identical(quantile(0.3, mixed_location, mixed_scale), 0)
})


test_that("pooled prediction summaries accept empty probabilities", {

  object <- brma(
    yi                        = c(-0.2, 0.1, 0.4),
    sei                       = c(0.2, 0.2, 0.2),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  posterior_samples <- cbind(mu = c(-0.2, 0.0, 0.2, 0.4), tau = c(0.1, 0.2, 0.3, 0.4))
  for (probs in list(NULL, numeric(0))) {
    pooled <- pooled_effect(object, probs = probs, .posterior_samples = posterior_samples)
    expect_identical(colnames(summary(pooled)), c("Mean", "Median"))
  }
})


test_that("pooled heterogeneity is the mean of the observation-specific tau", {

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
  # the arithmetic mean over the fitted rows, as at 4.0.0 (not exp(mean(log(tau))))
  expected  <- matrix(
    rowMeans(as.matrix(row_tau)),
    ncol     = 1L,
    dimnames = list(NULL, "tau")
  )
  geometric <- matrix(exp(rowMeans(log(as.matrix(row_tau)))), ncol = 1L)
  rms       <- matrix(sqrt(rowMeans(as.matrix(row_tau)^2)), ncol = 1L)

  pooled <- pooled_heterogeneity(
    object,
    .posterior_samples = posterior_samples
  )

  expect_equal(unname(as.matrix(pooled)), unname(expected), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(unname(expected), unname(geometric))))
  pooled_df  <- as.data.frame(summary(pooled))
  row_tau_df <- as.data.frame(summary(row_tau))
  expect_identical(
    names(pooled_df),
    c("component", "parameter", "Mean", "Median", "CI_0.025", "CI_0.975")
  )
  expect_identical(pooled_df[["component"]], "heterogeneity")
  expect_identical(pooled_df[["parameter"]], "tau")
  expect_identical(data.frame(summary(pooled)), pooled_df)
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
  expect_identical(data.frame(summary(row_tau)), row_tau_df)
  expect_equal(
    unname(as.matrix(row_tau_df)),
    unname(as.matrix(as.data.frame(summary(row_tau)))),
    tolerance = 0
  )
  expect_false(isTRUE(all.equal(unname(expected), unname(rms))))
})


test_that("prior-only scale formulas evaluate the fitted standardization", {

  # Objects without a fit evaluate their formulas through the designs that
  # fitting builds, which standardize continuous predictors. The design-less
  # evaluation used before read the raw predictor: for these draws its tau
  # was off by factors 0.58 to 2.07.
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

  # The scale formula has a log(intercept): log(tau) = log(intercept) + b x,
  # with x standardized by its mean and SD in the data.
  x_std    <- (dat[["x"]] - mean(dat[["x"]])) / stats::sd(dat[["x"]])
  expected <- exp(
    outer(log(posterior_samples[, "log_tau_intercept"]), rep(1, 3L)) +
      outer(posterior_samples[, "log_tau_x"], x_std)
  )
  row_tau <- predict(
    object,
    type               = "terms.scale",
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )
  expect_equal(unname(as.matrix(row_tau)), expected, tolerance = 1e-12)
  expect_equal(
    expected[1L, ],
    c(0.0994863170669, 0.1679628886922, 0.4787549728007),
    tolerance = 1e-12
  )

  # The pooled tau is the arithmetic mean of the observation-specific tau, as at
  # 4.0.0. The average standardized predictor is 0, so exp(mean(log(tau))) would
  # be the intercept; the arithmetic mean is larger.
  pooled <- pooled_heterogeneity(
    object,
    .posterior_samples = posterior_samples
  )
  expect_equal(unname(as.matrix(pooled))[, 1L], rowMeans(expected), tolerance = 1e-12)
  expect_true(all(
    unname(as.matrix(pooled))[, 1L] > unname(posterior_samples[, "log_tau_intercept"])
  ))
})

test_that("pooled target metadata names its fixed reference and leaves the RNG untouched", {

  object <- brma(yi = c(-.1, 0, .1), sei = rep(.1, 3), measure = "GEN",
    prior_effect = prior("point", list(location = 0)),
    prior_heterogeneity = prior("point", list(location = .1)),
    prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE)
  posterior <- cbind(mu = rep(0, 24L), tau = rep(.1, 24L))
  withr::local_seed(481)
  before <- .Random.seed
  pooled <- pooled_effect(object, bias_adjusted = FALSE, .posterior_samples = posterior)
  target <- attr(pooled, "RoBMA_target", exact = TRUE)
  expect_identical(target$reported_target, "fixed-location average-expanded-design reference")
  expect_identical(target$prediction_target, "new unselected true effect at the pooled reference")
  expect_null(target$rng_policy)
  expect_identical(.Random.seed, before)
})
