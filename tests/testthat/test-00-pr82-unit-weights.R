.unit_weights_selection_object <- function(weights = NULL) {

  dat <- data.frame(yi = c(.1, .2, .3, .4), sei = rep(.2, 4L),
                    study = factor(c("a", "a", "b", "b")))
  bselmodel(
    yi = yi, sei = sei, data = dat, weights = weights, measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE,
    prior_effect = BayesTools::prior("normal", list(0, 1)),
    prior_heterogeneity = BayesTools::prior("normal", list(0, 1),
                                            truncation = list(0, Inf)),
    prior_bias = BayesTools::prior_weightfunction(
      "one-sided", .025, BayesTools::wf_fixed(c(1, .2)),
      model = BayesTools::selection_model(
        group = "study", weight_rule = "best", known_sampling_variance = "integrate"
      )
    )
  )
}

.unit_weights_sampling_fixture <- function(weights = NULL) {

  dat <- data.frame(yi = c(.4, -.2), sei = c(.8, 1.1), paper = c("a", "a"))
  object <- bselmodel(
    yi = yi, sei = sei, data = dat, weights = weights, measure = "GEN",
    effect_direction = "positive", prior_unit_information_sd = 1,
    only_priors = TRUE, silent = TRUE,
    prior_effect = BayesTools::prior("normal", list(0, 1)),
    prior_heterogeneity = BayesTools::prior("normal", list(0, 1),
                                            truncation = list(lower = 0)),
    prior_bias = BayesTools::prior_weightfunction(
      "one-sided", .5, BayesTools::wf_fixed(c(1, .2)),
      model = BayesTools::selection_model(
        estimate_random_effects = "integrate", known_sampling_variance = "condition",
        group = "paper"
      )
    ),
    selection_control = set_selection_likelihood_control(relative_tolerance = 1e-6)
  )
  samples <- rbind(c(.2, .4, .5, -1, .3, -.4), c(-.1, .7, -.2, .6, -.5, .2))
  colnames(samples) <- c("mu", "tau", paste0("theta[", 1:2, "]"),
                         paste0("sampling_z[", 1:2, "]"))
  list(object = object, samples = samples)
}

test_that("explicit unit weights construct the same data as omitted weights", {

  dat <- data.frame(yi = c(.1, .2, .3), vi = c(.04, .09, .16))
  omitted <- brma.norm(yi = yi, vi = vi, data = dat, measure = "GEN", only_data = TRUE)
  unit    <- brma.norm(yi = yi, vi = vi, data = dat, weights = rep(1, 3L),
                       measure = "GEN", only_data = TRUE)
  expect_false(.is_weights(unit))
  expect_false(attr(unit[["data"]], "weights"))
  expect_identical(unit[["data"]][["outcome"]][["weights"]],
                   omitted[["data"]][["outcome"]][["weights"]])
  expect_identical(unit[["data"]], omitted[["data"]])
  expect_true(.is_weights(brma.norm(yi = yi, vi = vi, data = dat, weights = c(1, 2, 1),
                                    measure = "GEN", only_data = TRUE)))

  glmm <- function(...) brma.glmm(
    ai = c(1, 2, 3), ci = c(2, 3, 4), n1i = rep(10, 3L), n2i = rep(12, 3L),
    measure = "OR", only_data = TRUE, ...
  )
  expect_identical(glmm(weights = rep(1, 3L))[["data"]], glmm()[["data"]])

  # The criterion applies to the final rows, after subsetting and NA removal.
  dat_final <- data.frame(yi = c(.1, .2, NA, .4), vi = c(.04, .09, .16, .25),
                          w = c(1, 2, 3, 1))
  final <- suppressWarnings(brma.norm(
    yi = yi, vi = vi, data = dat_final, weights = w, subset = c(TRUE, FALSE, TRUE, TRUE),
    measure = "GEN", only_data = TRUE
  ))
  reference <- brma.norm(
    yi = yi, vi = vi, data = dat_final[c(1L, 4L), ], measure = "GEN", only_data = TRUE
  )
  expect_false(.is_weights(final))
  expect_identical(final[["data"]][["outcome"]][["weights"]],
                   reference[["data"]][["outcome"]][["weights"]])
})

test_that("integrated best selection accepts unit weights exactly as omission", {

  omitted <- .unit_weights_selection_object()
  unit    <- .unit_weights_selection_object(weights = rep(1, 4L))
  expect_s3_class(unit, "bselmodel")
  expect_false(.is_weights(unit))
  expect_identical(unit[["data"]], omitted[["data"]])
  expect_identical(.data_selection_execution_plan(unit[["data"]]),
                   .data_selection_execution_plan(omitted[["data"]]))

  # Non-unit powers keep the integrated-best gate.
  expect_error(
    .unit_weights_selection_object(weights = c(1, 2, 1, 1)),
    "'weights' are unavailable for jointly normalized selection vectors.",
    fixed = TRUE
  )
})

test_that("sampling-conditioned unit weights leave the covariance backend as omission does", {

  omitted <- .unit_weights_sampling_fixture()
  unit    <- .unit_weights_sampling_fixture(weights = rep(1, 2L))
  expect_identical(unit[["object"]][["data"]], omitted[["object"]][["data"]])

  evaluate <- function(fixture) {
    object <- fixture[["object"]]
    setup  <- .log_lik_posterior_setup(
      object[["fit"]], fixture[["samples"]], object[["data"]], object[["priors"]],
      "estimate", NULL
    )
    list(
      backend = .estimate_normal_target_uses_covariance_backend(object[["data"]], object[["priors"]]),
      point   = .log_lik_estimate_from_setup(setup),
      sum     = .log_lik_estimate_sum_from_setup(setup),
      tails   = .loo_predictive_log_tails_estimate(setup),
      cdf     = .cdf_lik_estimate.brma(object, setup = setup)
    )
  }
  expected <- evaluate(omitted)
  actual   <- evaluate(unit)
  expect_identical(actual[["backend"]], expected[["backend"]])
  expect_identical(actual[["point"]], expected[["point"]])
  expect_identical(actual[["sum"]], expected[["sum"]])
  expect_identical(actual[["tails"]], expected[["tails"]])
  expect_identical(actual[["cdf"]], expected[["cdf"]])
})
