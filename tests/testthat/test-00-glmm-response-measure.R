context("GLMM response measure prediction")


test_that("GLMM response count conversion matches continuity-corrected estimators", {

  bin_data <- data.frame(
    ai  = c(0, 0),
    ci  = c(0, 0),
    n1i = c(10, 20),
    n2i = c(12, 25)
  )
  bin_counts <- matrix(
    c(
      0, 3, 10, 25,
      4, 0, 15, 20
    ),
    nrow = 2,
    byrow = TRUE
  )

  bin_out <- .glmm_response_counts_to_measure(
    outcome_samples = bin_counts,
    outcome_data    = bin_data,
    outcome_type    = "bin"
  )

  bin_expected <- matrix(
    c(
      log(((0 + .5) * (12 - 3 + .5)) / ((10 - 0 + .5) * (3 + .5))),
      log(((10 + .5) * (25 - 25 + .5)) / ((20 - 10 + .5) * (25 + .5))),
      log(((4 + .5) * (12 - 0 + .5)) / ((10 - 4 + .5) * (0 + .5))),
      log((15 * (25 - 20)) / ((20 - 15) * 20))
    ),
    nrow = 2,
    byrow = TRUE
  )

  expect_equal(unname(bin_out), bin_expected, tolerance = 1e-12)
  expect_equal(colnames(bin_out), c("yi[1]", "yi[2]"))

  pois_data <- data.frame(
    x1i = c(0, 0),
    x2i = c(0, 0),
    t1i = c(10, 20),
    t2i = c(5, 40)
  )
  pois_counts <- matrix(
    c(
      0, 2, 4, 0,
      3, 6, 0, 8
    ),
    nrow = 2,
    byrow = TRUE
  )

  pois_out <- .glmm_response_counts_to_measure(
    outcome_samples = pois_counts,
    outcome_data    = pois_data,
    outcome_type    = "pois"
  )

  pois_expected <- matrix(
    c(
      log(((0 + .5) / 10) / ((2 + .5) / 5)),
      log(((4 + .5) / 20) / ((0 + .5) / 40)),
      log((3 / 10) / (6 / 5)),
      log(((0 + .5) / 20) / ((8 + .5) / 40))
    ),
    nrow = 2,
    byrow = TRUE
  )

  expect_equal(unname(pois_out), pois_expected, tolerance = 1e-12)
  expect_equal(colnames(pois_out), c("yi[1]", "yi[2]"))
})


test_that("GLMM response as_measure rejects undefined denominators", {

  expect_error(
    .check_glmm_response_as_measure(
      outcome_type = "bin",
      outcome_data = data.frame(n1i = c(10, 0), n2i = c(12, 12))
    ),
    "positive 'n1i' and 'n2i'"
  )

  expect_error(
    .check_glmm_response_as_measure(
      outcome_type = "pois",
      outcome_data = data.frame(t1i = c(10, 0), t2i = c(12, 12))
    ),
    "positive 't1i' and 't2i'"
  )

  expect_true(isTRUE(.check_glmm_response_as_measure(
    outcome_type = "pois",
    outcome_data = data.frame(t1i = c(10, 20), t2i = c(12, 15))
  )))
})


test_that("Poisson raw response RNG allows zero exposure as zero counts", {

  set.seed(1)

  out <- .outcome_rng.pois(
    mu_samples = matrix(c(0.2, -0.1, 0.3, 0.5), nrow = 2),
    log_phi    = matrix(0, nrow = 2, ncol = 2),
    t1i        = c(0, 3),
    t2i        = c(4, 0)
  )

  expect_equal(out[, "x1i[1]"], c(0, 0))
  expect_equal(out[, "x2i[2]"], c(0, 0))
  expect_true(all(out >= 0))
  expect_true(all(out == floor(out)))
})


test_that("raw GLMM response predictions omit effect-size metadata", {

  fitted_backend <- list()
  attr(fitted_backend, "prior_list") <- list(mu = NULL)
  outcome_data <- data.frame(
    ai  = c(0L, 0L),
    ci  = c(0L, 0L),
    n1i = c(10L, 12L),
    n2i = c(11L, 13L)
  )
  context <- list(
    object            = list(fit = fitted_backend),
    type              = "response",
    conditioning_depth = "marginal",
    as_measure        = FALSE,
    probs             = c(.025, .975),
    bias_adjusted     = TRUE,
    quiet             = TRUE,
    conditional       = FALSE,
    same_data         = FALSE,
    new_data          = list(outcome = outcome_data),
    known_V_new       = NULL,
    priors            = list(outcome = list(pi = NULL)),
    is_weightfunction = FALSE,
    is_known_v        = FALSE,
    outcome_type      = "bin",
    posterior_samples = matrix(0, nrow = 2L, ncol = 1L),
    effect_transform  = .effect_output_setup_measure(input_measure = "OR"),
    outcome_data      = outcome_data,
    K                 = 2L,
    n_chains          = 1L,
    n_iter            = 2L,
    random_mv         = FALSE
  )
  expected <- matrix(
    c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L),
    nrow = 2L,
    dimnames = list(NULL, c("ai[1]", "ci[1]", "ai[2]", "ci[2]"))
  )

  testthat::local_mocked_bindings(
    .predict_brma_estimate_draws = function(...) matrix(0, nrow = 2L, ncol = 2L),
    .evaluate.brma.baserate_newdata = function(...) matrix(0, nrow = 2L, ncol = 2L),
    .outcome_rng.binom = function(...) expected,
    .package = "RoBMA"
  )

  result <- .predict_brma_response(
    context        = context,
    location_state = list(mu = matrix(0, nrow = 2L, ncol = 2L)),
    scale_state    = list(within = matrix(0, nrow = 2L, ncol = 2L))
  )

  expect_s3_class(result, "brma_samples")
  expect_equal(unname(as.matrix(result)), unname(expected))
  expect_identical(attr(result, "title"), "Observed Counts Posterior Prediction:")
  expect_null(attr(result, "effect_transform", exact = TRUE))
  expect_null(attr(result, "footnotes", exact = TRUE))
})


test_that("GLMM response nuisance rates follow conditioning depth", {

  fitted_backend <- list()
  attr(fitted_backend, "prior_list") <- list(mu = NULL)
  outcome_data <- data.frame(
    ai  = 0L,
    ci  = 0L,
    n1i = 10L,
    n2i = 11L
  )
  context <- list(
    object             = list(
      fit  = fitted_backend,
      data = list(outcome = outcome_data)
    ),
    type               = "response",
    conditioning_depth = "marginal",
    as_measure         = FALSE,
    probs              = c(.025, .975),
    bias_adjusted      = TRUE,
    quiet              = TRUE,
    conditional        = FALSE,
    same_data          = TRUE,
    new_data           = list(outcome = outcome_data),
    known_V_new        = NULL,
    priors             = list(outcome = list(pi = NULL)),
    is_weightfunction  = FALSE,
    is_known_v         = FALSE,
    outcome_type       = "bin",
    posterior_samples  = matrix(0, nrow = 2L, ncol = 1L),
    effect_transform   = .effect_output_setup_measure(input_measure = "OR"),
    outcome_data       = outcome_data,
    K                  = 1L,
    n_chains           = 1L,
    n_iter             = 2L,
    random_mv          = FALSE
  )

  testthat::local_mocked_bindings(
    .predict_brma_estimate_draws = function(...) matrix(0, nrow = 2L, ncol = 1L),
    .evaluate.brma.baserate = function(...) matrix(11, nrow = 2L, ncol = 1L),
    .evaluate.brma.baserate_newdata = function(...) matrix(-7, nrow = 2L, ncol = 1L),
    .outcome_rng.binom = function(mu_samples, logit_baserate, ...) logit_baserate,
    .package = "RoBMA"
  )

  predict_at_depth <- function(depth) {
    context[["conditioning_depth"]] <- depth
    as.matrix(.predict_brma_response(
      context        = context,
      location_state = list(mu = matrix(0, nrow = 2L, ncol = 1L)),
      scale_state    = list(within = matrix(0, nrow = 2L, ncol = 1L))
    ))
  }

  expect_equal(unname(predict_at_depth("estimate")), matrix(11, 2L, 1L))
  expect_equal(unname(predict_at_depth("marginal")), matrix(-7, 2L, 1L))
  expect_equal(unname(predict_at_depth("cluster")), matrix(-7, 2L, 1L))
})
