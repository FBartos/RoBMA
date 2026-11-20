test_that("univariate zplot density and EDR mix eligible and exempt outcomes", {

  yi <- c(.1, .4)
  sei <- c(.2, .3)
  mu <- matrix(c(.1, -.2), 1L)
  sd <- matrix(sei, 1L)
  prior <- BayesTools::prior_weightfunction("one-sided", .05, BayesTools::wf_fixed(c(1, .25)))
  cutoff <- stats::qnorm(.95)
  z <- c(-1, 2)
  for (exempt in list(c(TRUE, FALSE), rep(TRUE, 2L))) {
    object <- bselmodel(yi = yi, sei = sei, measure = "SMD", effect_direction = "positive",
                        prior_bias = prior, bias_exempt = exempt, only_priors = TRUE)
    draws <- matrix(c(0, 0), 1L, dimnames = list(NULL, c("mu", "tau")))
    selection <- .selection_context(object, posterior_samples = draws)
    predictive <- list(mu = mu, mu_extrapolated = mu, sei = sei)
    control <- set_selection_likelihood_control()
    location <- as.numeric(mu) / sei
    low_weight <- ifelse(exempt, 1, .25)
    normalizer <- low_weight + (1 - low_weight) * stats::pnorm(cutoff - location, lower.tail = FALSE)
    expect_true(.zplot_requires_selection_marginal(object, "marginal"))
    density <- .zplot_joint_marginal(object, draws, predictive, selection, z, FALSE, control)
    expected <- vapply(z, function(value) {
      weight <- if (value >= cutoff) 1 else low_weight
      mean(stats::dnorm(value, location) * weight / normalizer)
    }, numeric(1L))
    expect_equal(as.numeric(density$fitted), expected, tolerance = 1e-12)
    expect_equal(density$weights, mean(1 / normalizer), tolerance = 1e-12)
    tail <- .zplot_joint_marginal(object, draws, predictive, selection, 2, TRUE, control)
    left <- stats::pnorm(-2 - location)
    right <- stats::pnorm(2 - location, lower.tail = FALSE)
    expect_equal(as.numeric(tail$fitted), mean((low_weight * left + right) / normalizer), tolerance = 1e-12)
    expect_equal(tail$EDR, stats::weighted.mean(left + right, 1 / normalizer), tolerance = 1e-12)
  }
})

test_that("multivariate zplots preserve correlation with an exempt outcome", {

  rho <- .6
  V <- matrix(c(1, rho, rho, 1), 2L)
  dat <- data.frame(yi = c(.1, .2), paper = c("p", "p"), registered = c(TRUE, FALSE))
  z <- c(-1, 2)
  cutoff <- stats::qnorm(.95)
  normalizer <- .05 + .95 * .25
  # Marginalizing the eligible outcome gives the selection probability
  # conditional on the exempt outcome's z-score.
  selected_density <- function(value) {
    p <- stats::pnorm(cutoff, rho * value, sqrt(1 - rho^2), lower.tail = FALSE)
    stats::dnorm(value) * (.25 + .75 * p + ifelse(value >= cutoff, 1, .25)) / (2 * normalizer)
  }
  expected_tail <- stats::integrate(selected_density, -Inf, -2)$value +
    stats::integrate(selected_density, 2, Inf)$value
  for (rule in c("product", "best")) {
    prior <- BayesTools::prior_weightfunction("one-sided", .05, BayesTools::wf_fixed(c(1, .25)),
      model = selection_model(weight_rule = rule, group = paper))
    object <- bselmodel.mv(yi = yi, V = V, data = dat, measure = "SMD", effect_direction = "positive",
                           prior_bias = prior, bias_exempt = registered, only_priors = TRUE)
    draws <- matrix(0, 1L, 1L, dimnames = list(NULL, "mu"))
    selection <- .selection_context(object, posterior_samples = draws)
    predictive <- list(mu = matrix(0, 1L, 2L), mu_extrapolated = matrix(0, 1L, 2L), sei = c(1, 1))
    control <- set_selection_likelihood_control(relative_tolerance = .001)
    density <- .zplot_joint_marginal(object, draws, predictive, selection, z, FALSE, control)
    expect_equal(as.numeric(density$fitted), selected_density(z), tolerance = .001)
    tail <- .zplot_joint_marginal(object, draws, predictive, selection, 2, TRUE, control)
    expect_equal(as.numeric(tail$fitted), expected_tail, tolerance = .001)
    expect_equal(tail$EDR, 2 * stats::pnorm(-2), tolerance = 1e-12)
  }
})
