.weighted_selection_fixture <- function(mode, weights = NULL, direction = "positive") {

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
        known_sampling_variance = "integrate", group = "paper", weight_rule = "product"
      )
    )
  )
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
