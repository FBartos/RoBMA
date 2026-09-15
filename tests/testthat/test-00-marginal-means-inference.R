context("Marginal-means inference")


.fixed_marginal_posterior <- function(value = 0, n = 20L) {

  posterior <- structure(
    rep(value, n),
    class = c("marginal_posterior", "numeric")
  )
  attr(posterior, "posterior_atoms") <- BayesTools::posterior_atom_attribute(
    data.frame(x = value, mass = 1)
  )
  return(posterior)
}


test_that("BF-free marginal means retain structurally fixed cells", {

  posterior <- list(zero = .fixed_marginal_posterior())
  inference <- .marginal_means_inclusion_bf(
    posterior       = posterior,
    null_hypothesis = 0,
    compute         = FALSE
  )

  expect_named(inference, "zero")
  expect_true(is.na(inference[["zero"]]))
})


test_that("structurally fixed marginal cells have unavailable BFs", {

  posterior <- list(zero = .fixed_marginal_posterior())
  testthat::local_mocked_bindings(
    Savage_Dickey_BF = function(...) {
      stop(
        paste0(
          "The posterior contains a declared point mass at the exact null ",
          "hypothesis value. The ordinary Savage-Dickey density ratio is invalid."
        ),
        call. = FALSE
      )
    },
    .package = "BayesTools"
  )
  inference <- .marginal_means_inclusion_bf(
    posterior       = posterior,
    null_hypothesis = 0,
    compute         = TRUE
  )

  expect_named(inference, "zero")
  expect_true(is.na(inference[["zero"]]))
  expect_match(
    attr(inference[["zero"]], "warnings", exact = TRUE),
    "structurally fixed"
  )
})


test_that("interaction marginals condition on every contributing coefficient", {

  terms <- c("intercept", "a", "b", "a:b")
  parameters <- c("mu_intercept", "mu_a", "mu_b", "mu_ab")
  spike_and_slab_priors <- stats::setNames(lapply(parameters, function(parameter) {

    BayesTools::prior_spike_and_slab(
      prior_parameter = BayesTools::prior("normal", list(mean = 0, sd = 1))
    )
  }), parameters)
  object <- list(fit = structure(list(), prior_list = spike_and_slab_priors))
  conditional_list <- .marginal_means_conditional_list(
    object     = object,
    terms      = terms,
    parameters = parameters
  )

  cell_weights <- rbind(
    base  = c(1, 0, 0, 0),
    a     = c(1, 1, 0, 0),
    b     = c(1, 0, 1, 0),
    `a:b` = c(1, 1, 1, 1)
  )
  colnames(cell_weights) <- parameters
  marginal <- lapply(seq_len(nrow(cell_weights)), function(i) {

    structure(
      numeric(20),
      linear_weights = cell_weights[i, ]
    )
  })
  names(marginal) <- rownames(cell_weights)
  prior_list <- stats::setNames(lapply(parameters, function(parameter) {

    BayesTools::prior("normal", parameters = list(mean = 0, sd = 1))
  }), parameters)

  effective <- BayesTools:::.marginal_inference_level_conditionals(
    marginal    = marginal,
    prior_list  = prior_list,
    conditional = conditional_list[["mu_ab"]]
  )
  expected <- lapply(seq_len(nrow(cell_weights)), function(i) {

    names(cell_weights[i, ])[cell_weights[i, ] != 0]
  })
  names(expected) <- rownames(cell_weights)
  expect_identical(effective, expected)

  # Independent finite-mixture reference for a point-null BF in the interaction
  # cell. Each active coefficient contributes an independent Normal variate.
  states <- as.matrix(expand.grid(
    mu_intercept = 0:1,
    mu_a         = 0:1,
    mu_b         = 0:1,
    mu_ab        = 0:1
  ))
  coefficient_mean <- c(
    mu_intercept = 0.4,
    mu_a         = -1.5,
    mu_b         = 1.8,
    mu_ab        = 0.8
  )
  coefficient_sd <- c(
    mu_intercept = 0.6,
    mu_a         = 0.2,
    mu_b         = 0.2,
    mu_ab        = 0.5
  )
  prior_probability <- rep(1 / nrow(states), nrow(states))
  posterior_log_weight <- as.numeric(
    states %*% c(mu_intercept = -1, mu_a = 2.5, mu_b = 2, mu_ab = -2)
  ) + 2 * states[, "mu_a"] * states[, "mu_b"]
  posterior_probability <- exp(
    posterior_log_weight - max(posterior_log_weight)
  )
  posterior_probability <- posterior_probability / sum(posterior_probability)

  conditional_density <- function(probability, conditional) {

    event <- rowSums(states[, conditional, drop = FALSE]) > 0
    state_mean <- as.numeric(states %*% coefficient_mean)
    state_sd <- sqrt(as.numeric((states^2) %*% (coefficient_sd^2)))
    event_probability <- probability[event] / sum(probability[event])
    sum(event_probability * stats::dnorm(
      0,
      mean = state_mean[event],
      sd   = state_sd[event]
    ))
  }
  mixture_bf <- function(conditional) {

    conditional_density(prior_probability, conditional) /
      conditional_density(posterior_probability, conditional)
  }

  contributing <- names(cell_weights["a:b", ])[cell_weights["a:b", ] != 0]
  actual_bf     <- mixture_bf(effective[["a:b"]])
  reference_bf  <- mixture_bf(contributing)
  incomplete_bf <- mixture_bf(c("mu_intercept", "mu_ab"))

  expect_equal(actual_bf, reference_bf, tolerance = 1e-14)
  expect_equal(actual_bf, 0.3565699, tolerance = 1e-7)
  expect_gt(abs(log(actual_bf / incomplete_bf)), log(2))
})


test_that("only coefficients with inclusion indicators are conditional candidates", {

  parameters <- c("mu_intercept", "mu_a", "mu_b")
  terms      <- c("intercept", "a", "b")
  make_object <- function(prior_list) {

    list(fit = structure(list(), prior_list = prior_list))
  }
  plain          <- BayesTools::prior("normal", list(mean = 0, sd = 1))
  spike_and_slab <- BayesTools::prior_spike_and_slab(prior_parameter = plain)
  null_mixture   <- BayesTools::prior_mixture(
    list(
      BayesTools::prior("spike",  list(location = 0)),
      BayesTools::prior("normal", list(mean = 0, sd = 1))
    ),
    components = c("null", "alternative")
  )

  # Plain priors are always included in every model, so BayesTools rejects them
  # as conditional labels; they must not be offered as candidates.
  plain_list <- .marginal_means_conditional_list(
    object     = make_object(stats::setNames(
      list(plain, plain, plain), parameters
    )),
    terms      = terms,
    parameters = parameters
  )
  expect_named(plain_list, parameters)
  expect_identical(plain_list[["mu_a"]], character(0))

  mixed_list <- .marginal_means_conditional_list(
    object     = make_object(stats::setNames(
      list(plain, spike_and_slab, null_mixture), parameters
    )),
    terms      = terms,
    parameters = parameters
  )
  expect_identical(mixed_list[["mu_intercept"]], c("mu_a", "mu_b"))
  expect_identical(mixed_list[["mu_b"]], c("mu_a", "mu_b"))

  # A parameter missing from the prior list is not a conditional candidate.
  partial_list <- .marginal_means_conditional_list(
    object     = make_object(list(mu_a = spike_and_slab)),
    terms      = terms,
    parameters = parameters
  )
  expect_identical(partial_list[["mu_a"]], "mu_a")
})
