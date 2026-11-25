.weighted_selection_fixture <- function(mode, weights = NULL, direction = "positive",
                                        weight_rule = "product",
                                        known_sampling_variance = "integrate",
                                        setup = TRUE) {

  dat <- data.frame(yi = c(.4, -.2, .3), sei = c(.8, 1.1, .6), paper = rep("a", 3L))
  object <- bselmodel(
    yi = yi, sei = sei, data = dat, weights = weights, measure = "GEN",
    effect_direction = direction, prior_unit_information_sd = 1,
    only_priors = TRUE, silent = TRUE,
    prior_effect = BayesTools::prior("normal", list(0, 1)),
    prior_heterogeneity = BayesTools::prior("normal", list(0, 1),
                                            truncation = list(lower = 0)),
    prior_bias = BayesTools::prior_weightfunction(
      "one-sided", .5, BayesTools::wf_fixed(c(1, .2)),
      model = BayesTools::selection_model(
        estimate_random_effects = mode, other_random_effects = "integrate",
        known_sampling_variance = known_sampling_variance, group = "paper",
        weight_rule = weight_rule
      )
    )
  )
  if (!setup) {
    return(list(object = object, dat = dat))
  }
  samples <- rbind(c(.2, .4, .5, -1, 1.5), c(-.1, .7, -1, .2, .6))
  colnames(samples) <- c("mu", "tau", paste0("theta[", 1:3, "]"))
  if (mode == "integrate") samples <- samples[, c("mu", "tau"), drop = FALSE]
  setup <- .log_lik_posterior_setup(
    object[["fit"]], samples, object[["data"]], object[["priors"]], "estimate", NULL
  )
  list(object = object, dat = dat, samples = samples, setup = setup)
}

# The selected-normal law of one observation: Gaussian around the candidate
# mean, one-sided selection at the origin with weights (1, .2) below it.
.weighted_selection_reference <- function(fixture, mode) {

  samples <- fixture[["samples"]]
  dat     <- fixture[["dat"]]
  mean    <- matrix(samples[, "mu"], 2L, 3L)
  sd      <- matrix(dat[["sei"]], 2L, 3L, byrow = TRUE)
  if (mode == "condition") {
    mean <- mean + unname(samples[, 3:5]) * unname(samples[, "tau"])
  } else {
    sd <- sqrt(sd^2 + unname(samples[, "tau"])^2)
  }
  yi   <- matrix(dat[["yi"]], 2L, 3L, byrow = TRUE)
  mass <- .2 + .8 * stats::pnorm(mean / sd)
  p    <- stats::pnorm((yi - mean) / sd)
  p0   <- stats::pnorm(-mean / sd)
  cdf  <- ifelse(yi < 0, .2 * p / mass, (.2 * p0 + p - p0) / mass)
  list(
    score = stats::dnorm(yi, mean, sd, log = TRUE) +
      matrix(log(ifelse(dat[["yi"]] >= 0, 1, .2)), 2L, 3L, byrow = TRUE) - log(mass),
    cdf   = cdf
  )
}

test_that("weighted joint-selection pointwise scores count the estimate effect once", {

  powers <- c(2, .5, 3)
  for (mode in c("condition", "integrate")) {
    fixture   <- .weighted_selection_fixture(mode, weights = powers)
    reference <- .weighted_selection_reference(fixture, mode)
    expected  <- sweep(reference[["score"]], 2L, powers, "*")

    # Powers multiply the selected-normal scores once; the pointwise and
    # summed routes agree.
    expect_equal(unname(.log_lik_estimate_from_setup(fixture[["setup"]])), expected,
                 tolerance = 1e-12)
    expect_equal(unname(.log_lik_estimate_sum_from_setup(fixture[["setup"]])),
                 rowSums(expected), tolerance = 1e-12)
  }
})

test_that("uniform powers of two multiply the selected-normal scores", {

  fixture   <- .weighted_selection_fixture("condition", weights = c(2, 2, 2))
  reference <- .weighted_selection_reference(fixture, "condition")
  expect_equal(unname(.log_lik_estimate_from_setup(fixture[["setup"]])),
               2 * reference[["score"]], tolerance = 1e-12)
})

test_that("weighted joint-selection LOO tails and moments follow the selected law", {

  for (mode in c("condition", "integrate")) {
    fixture   <- .weighted_selection_fixture(mode, weights = c(2, .5, 3))
    reference <- .weighted_selection_reference(fixture, mode)
    setup     <- fixture[["setup"]]

    # The law, its tails and its CDF do not depend on the powers.
    tails <- .loo_predictive_log_tails_estimate(setup)
    expect_equal(unname(exp(tails[["log_lower"]])), reference[["cdf"]], tolerance = 1e-12)
    expect_equal(unname(exp(tails[["log_upper"]])), 1 - reference[["cdf"]],
                 tolerance = 1e-12)

    # Predictive moments are those of the unweighted law.
    unweighted   <- .weighted_selection_fixture(mode)
    psis_weights <- matrix(1 / 2, 2L, 3L)
    expect_equal(
      .loo_predictive_moments_estimate(
        fixture[["object"]], setup, psis_object = NULL, psis_weights = psis_weights
      ),
      .loo_predictive_moments_estimate(
        unweighted[["object"]], unweighted[["setup"]], psis_object = NULL,
        psis_weights = psis_weights
      ),
      tolerance = 1e-12
    )
  }
})


# The PSIS-weighted LOO-PIT z, residual and SE of the selected-normal law of
# .weighted_selection_reference() in either effect direction, by integration.
.weighted_selection_loo_pit_reference <- function(fixture, mode, sign, psis_weights) {

  samples <- fixture[["samples"]]
  dat     <- fixture[["dat"]]
  mean    <- matrix(samples[, "mu"], 2L, 3L)
  sd      <- matrix(dat[["sei"]], 2L, 3L, byrow = TRUE)
  if (mode == "condition") {
    mean <- mean + unname(samples[, 3:5]) * unname(samples[, "tau"])
  } else {
    sd <- sqrt(sd^2 + unname(samples[, "tau"])^2)
  }
  cdf <- first <- variance <- matrix(NA_real_, 2L, 3L)
  for (s in 1:2) for (k in 1:3) {
    # Selection acts on the signed outcome: weight 1 above zero, .2 below.
    m       <- sign * mean[s, k]
    v       <- sd[s, k]
    density <- function(x, power = 0) x^power * stats::dnorm(x, m, v) * ifelse(x > 0, 1, .2)
    moment  <- function(power) sum(vapply(list(c(-Inf, 0), c(0, Inf)), function(range) {
      stats::integrate(density, range[1L], range[2L], power = power, rel.tol = 1e-12)$value
    }, numeric(1L)))
    mass    <- moment(0)
    signed_cdf <- function(t) {
      lower <- stats::integrate(density, -Inf, min(t, 0), rel.tol = 1e-12)$value
      upper <- if (t > 0) stats::integrate(density, 0, t, rel.tol = 1e-12)$value else 0
      (lower + upper) / mass
    }
    y <- dat[["yi"]][k]
    cdf[s, k]      <- if (sign == 1) signed_cdf(y) else 1 - signed_cdf(-y)
    first[s, k]    <- sign * moment(1) / mass
    variance[s, k] <- moment(2) / mass - (moment(1) / mass)^2
  }
  weights <- sweep(psis_weights, 2L, colSums(psis_weights), "/")
  center  <- colSums(weights * first)
  list(
    z     = stats::qnorm(colSums(weights * cdf)),
    resid = dat[["yi"]] - center,
    se    = sqrt(colSums(weights * (variance + sweep(first, 2L, center)^2)))
  )
}

test_that("the native product-weight LOO-PIT reducer follows effect direction and likelihood weights", {

  psis_weights <- matrix(c(.3, .7, .8, .2, .5, .5), 2L, 3L)
  for (direction in c("positive", "negative")) {
    for (mode in c("condition", "integrate")) {
      for (weights in list(NULL, c(2, .5, 3))) {
        fixture   <- .weighted_selection_fixture(mode, weights = weights, direction = direction)
        reference <- .weighted_selection_loo_pit_reference(
          fixture, mode, if (direction == "positive") 1 else -1, psis_weights
        )
        actual <- .loo_predictive_joint_product_summary(fixture[["setup"]], psis_weights)
        info   <- paste(direction, mode, if (is.null(weights)) "unweighted" else "weighted")
        expect_equal(actual[["z"]], reference[["z"]], tolerance = 1e-9, info = info)
        expect_equal(actual[["resid"]], reference[["resid"]], tolerance = 1e-9, info = info)
        expect_equal(actual[["se"]], reference[["se"]], tolerance = 1e-9, info = info)
      }
    }
  }
})

test_that("the native LOO-PIT reducer is taken only for product weights without retained sampling", {

  routing <- function(fixture) {
    data <- fixture[["object"]][["data"]]
    .loo_predictive_joint_product_summary(
      list(outcome_type = "norm", data = data, weights = NULL),
      psis_weights = matrix(1, 1L, 3L)
    )
  }
  product  <- .weighted_selection_fixture("integrate", setup = FALSE)
  best     <- .weighted_selection_fixture("integrate", weight_rule = "best", setup = FALSE)
  retained <- .weighted_selection_fixture(
    "integrate", known_sampling_variance = "condition", setup = FALSE
  )
  expect_true(.selection_joint_product_weights(product[["object"]][["data"]]))
  expect_false(.selection_joint_product_weights(best[["object"]][["data"]]))
  expect_true(.selection_retains_sampling(retained[["object"]][["data"]]))
  expect_null(routing(best))
  expect_null(routing(retained))
  expect_null(.loo_predictive_joint_product_summary(
    list(outcome_type = "bin", data = product[["object"]][["data"]], weights = NULL),
    psis_weights = matrix(1, 1L, 3L)
  ))
})
