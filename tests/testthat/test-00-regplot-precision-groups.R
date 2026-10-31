.regplot_precision_replay <- function(mods = ~ x * sei + vi, binomial = FALSE) {

  dat <- data.frame(
    yi = seq(.2, .7, by = .1), sei = rep(c(.1, .2, .3), each = 2),
    x = rep(c(0, 1), 3), ai = 4:9, ci = 3:8, n1i = 20, n2i = 20
  )
  dat$vi <- dat$sei^2
  normal <- BayesTools::prior("normal", list(mean = 0, sd = 1))
  args <- list(mods = mods, data = dat, prior_unit_information_sd = 1,
    standardize_continuous_predictors = FALSE, prior_effect = normal,
    prior_mods = normal, prior_heterogeneity = BayesTools::prior("normal",
      list(mean = 0, sd = 1), list(lower = 0, upper = Inf)),
    only_priors = TRUE, silent = TRUE)
  object <- if (binomial) {
    do.call(brma.glmm, c(args, list(ai = quote(ai), ci = quote(ci),
      n1i = quote(n1i), n2i = quote(n2i), measure = "OR",
      prior_baserate = BayesTools::prior("beta", list(alpha = 1, beta = 1)))))
  } else {
    do.call(brma, c(args, list(yi = quote(yi), sei = quote(sei), measure = "GEN")))
  }
  setup <- review_compile_object(object)
  priors <- c(do.call(c, unname(lapply(setup$compiled, `[[`, "prior_list"))),
    .create_fit_priors(object$data, object$priors))
  samples <- matrix(0, 24L, length(review_prior_columns(priors)),
    dimnames = list(NULL, review_prior_columns(priors)))
  for (name in names(priors)) {
    prior <- priors[[name]]
    stopifnot(prior$distribution %in% c("normal", "beta"))
    value <- if (prior$distribution == "beta") {
      prior$parameters$alpha / (prior$parameters$alpha + prior$parameters$beta)
    } else prior$parameters$mean
    if (name == "tau") value <- .5
    if (name == review_coefficient_name(setup, "sei")) value <- 1
    samples[, BayesTools:::.prior_linear_prior_columns(name, prior)] <- value
  }
  object <- review_attach_replay(setup, samples)
  BayesTools::JAGS_validate_fit_contract(object$fit, requires = c(
    "name_encoding", "formula_name_map", "formula_design", "parameter_map", "draw_geometry"))
  list(object = object, data = dat)
}

test_that("regplot retains declared precision grouping in normal and binomial models", {

  for (binomial in c(FALSE, TRUE)) {
    fixture <- .regplot_precision_replay(binomial = binomial)
    for (by in c("sei", "vi")) {
      values <- fixture$data[[by]]
      groups <- c(mean(values) - stats::sd(values), mean(values), mean(values) + stats::sd(values))
      expected <- rep(if (by == "sei") groups else sqrt(groups), each = 2L)
      coherent <- data.frame(x = rep(c(1, 2), 3), sei = expected, vi = expected^2)
      if (by == "vi") coherent$vi <- rep(groups, each = 2L)
      if (binomial) {
        coherent$ai <- coherent$ci <- coherent$n1i <- coherent$n2i <- 0
      } else coherent$yi <- 0
      for (sampling_bias in c(TRUE, FALSE)) {
        direct <- as.matrix(predict(fixture$object, newdata = coherent, type = "terms",
          bias_adjusted = !sampling_bias, quiet = TRUE))
        actual <- regplot(fixture$object, mod = "x", by = by, at = c(1, 2),
          sampling_bias = sampling_bias, as_data = TRUE)
        expect_identical(max(abs(sweep(direct, 2L, expected, "-"))), 0)
        expect_lte(max(abs(actual$pred$y - expected)), 8 * .Machine$double.eps)
        expect_identical(actual$sei, expected)
        expect_identical(actual$pred$x, rep(c(1, 2), 3))
        expect_identical(actual$groups, paste0(by, ": ", c("Mean - 1 SD", "Mean", "Mean + 1 SD")))
      }
    }
    explicit <- regplot(fixture$object, mod = "x", by = "sei", at = c(1, 2), as_data = TRUE)
    automatic <- regplot(fixture$object, mod = "x", at = c(1, 2), as_data = TRUE)
    expect_identical(automatic, explicit)
    message <- "'sei' cannot be supplied when the grouping moderator is 'sei' or 'vi'."
    expect_error(regplot(fixture$object, mod = "x", by = "sei", sei = .3,
      as_data = TRUE), message, fixed = TRUE)
    expect_error(regplot(fixture$object, mod = "x", by = "vi", sei = .3,
      as_data = TRUE), message, fixed = TRUE)
    expect_error(regplot(fixture$object, mod = "x", sei = .3,
      as_data = TRUE), message, fixed = TRUE)
  }
})

test_that("regplot retains ordinary precision references without precision grouping", {

  fixture <- .regplot_precision_replay(mods = ~ x + sei + vi)
  expected <- sqrt(mean(fixture$data$vi))
  for (sampling_bias in c(TRUE, FALSE)) {
    actual <- regplot(fixture$object, mod = "x", at = c(1, 2),
      sampling_bias = sampling_bias, as_data = TRUE)
    expect_null(actual$by_name)
    expect_identical(actual$sei, expected)
    expect_lte(max(abs(actual$pred$y - expected)), 8 * .Machine$double.eps)
    override <- regplot(fixture$object, mod = "x", at = c(1, 2), sei = .3,
      sampling_bias = sampling_bias, as_data = TRUE)
    expect_identical(override$sei, .3)
    expect_lte(max(abs(override$pred$y - .3)), 8 * .Machine$double.eps)
  }
})

test_that("group precision validation retains resolved fields without a second conversion", {

  sei <- c(.1, .2, .3)
  newdata <- data.frame(x = 1:3, sei = sei, vi = .075)
  actual <- .regplot_prediction_precision(newdata, "x", 3L, .25, NULL, by_name = "sei")
  expect_identical(actual$newdata$sei, sei)
  expect_identical(actual$newdata$vi, sei^2)
  fixture <- .regplot_precision_replay()
  dummy <- .regplot_add_dummy_outcome(fixture$object, actual$newdata, 3L, sei)
  expect_identical(dummy$sei, sei)
  expect_identical(dummy$vi, sei^2)
  for (invalid in c(-.1, Inf, NA_real_)) {
    newdata$sei[1L] <- invalid
    expect_error(.regplot_prediction_precision(newdata, "x", 3L, .25, NULL,
      by_name = "sei"), "Prediction standard errors must be finite and non-negative.", fixed = TRUE)
  }
  newdata <- data.frame(x = 1:3, sei = .2, vi = c(.01, .04, .09))
  for (invalid in c(-.1, Inf, NA_real_)) {
    newdata$vi[1L] <- invalid
    expect_error(.regplot_prediction_precision(newdata, "x", 3L, .25, NULL,
      by_name = "vi"), "Prediction variances must be finite and non-negative.", fixed = TRUE)
  }
})

test_that("grouped precision supplies normal sampling intervals without sampling bias", {

  fixture <- .regplot_precision_replay()
  level <- 95
  probabilities <- c((1 - level / 100) / 2, 1 - (1 - level / 100) / 2)
  for (by in c("sei", "vi")) {
    values <- fixture$data[[by]]
    groups <- c(mean(values) - stats::sd(values), mean(values), mean(values) + stats::sd(values))
    group_sei <- if (by == "sei") groups else sqrt(groups)
    actual <- regplot(fixture$object, mod = "x", by = by, at = c(1, 2),
      si = TRUE, level = level, sampling_bias = FALSE, as_data = TRUE)
    sei <- group_sei[actual$si$group_id]
    sd <- sqrt(.25 + sei^2)
    expected_lower <- stats::qnorm(probabilities[[1L]], mean = sei, sd = sd)
    expected_upper <- stats::qnorm(probabilities[[2L]], mean = sei, sd = sd)
    # The native root's local tolerance is 2*eps*abs(root)+eps*sd. This
    # allowance also covers qnorm/CDF rounding and averaging24 equal states.
    lower_scale <- pmax(1, abs(expected_lower), sd)
    upper_scale <- pmax(1, abs(expected_upper), sd)
    expect_lte(max(abs(actual$si$lower - expected_lower) / lower_scale),
      64 * .Machine$double.eps)
    expect_lte(max(abs(actual$si$upper - expected_upper) / upper_scale),
      64 * .Machine$double.eps)
    expect_identical(actual$sei, rep(group_sei, each = 2L))
  }
})

test_that("continuous regplot sorts duplicate evaluations with every interval", {

  fixture <- .regplot_precision_replay(mods = ~ sei + vi)
  x <- c(.1, .1, .2, .3)
  actual <- regplot(fixture$object, mod = "sei", at = c(.3, .1, .2, .1),
    pi = TRUE, si = TRUE, as_data = TRUE)
  expect_identical(actual$pred$x, x)
  expect_equal(actual$pred$y, x, tolerance = 0)
  expect_identical(actual$sei, x)
  for (band in c("ci", "pi", "si")) {
    expect_identical(actual[[band]]$x, c(x, rev(x)))
  }
  q <- stats::qnorm(c(.025, .975))
  expected_pi <- c(x + .5 * q[1L], rev(x + .5 * q[2L]))
  expected_si <- c(x + sqrt(.25 + x^2) * q[1L],
    rev(x + sqrt(.25 + x^2) * q[2L]))
  # The native inverse has its existing few-ulp local stopping allowance.
  expect_lte(max(abs(actual$pi$y - expected_pi)), 64 * .Machine$double.eps)
  expect_lte(max(abs(actual$si$y - expected_si)), 64 * .Machine$double.eps)
  for (mod in list(c(1, 2), NA_real_, Inf, 1.5, 0)) {
    expect_error(regplot(fixture$object, mod = mod, as_data = TRUE), "moderator index", fixed = TRUE)
  }
  for (level in c(NA_real_, Inf, -Inf)) {
    expect_error(regplot(fixture$object, level = level, as_data = TRUE), "'level'", fixed = TRUE)
  }
  for (sei in c(NA_real_, Inf, -Inf)) {
    expect_error(regplot(fixture$object, sei = sei, as_data = TRUE), "'sei'", fixed = TRUE)
  }
})

test_that("GLMM regression bands disclose their descriptive normal response law", {

  fixture <- .regplot_precision_replay(mods = ~ sei + vi, binomial = TRUE)
  expect_warning(actual <- regplot(fixture$object, mod = "sei", at = c(.1, .2),
    pi = TRUE, si = TRUE, as_data = TRUE), paste0(
      "GLMM regression prediction and sampling bands are descriptive normal ",
      "effect-size approximations; they are not coverage intervals from the ",
      "fitted discrete likelihood."), fixed = TRUE)
  expect_identical(actual$approximation, list(method = "normal_effect_size_approximation",
    descriptive = TRUE, fitted_discrete_likelihood_coverage = FALSE))
})
