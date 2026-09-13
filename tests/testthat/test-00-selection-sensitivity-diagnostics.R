test_that("selection conditioning comparisons are requested explicitly", {

  calls <- list()
  warning_text <- paste0(
    "Selection conditioning sensitivity was substantial: the largest unit-level ",
    "posterior-median total variation distance was 23.5%. ",
    "Retaining and integrating the specified contexts define different reporting models. ",
    "Inspect 'selection_sensitivity_diagnostics()' for the fitted and reference specifications."
  )
  testthat::local_mocked_bindings(
    .fit = function(object, extend = FALSE) {
      list(has_posterior = TRUE, generation = if (extend) 2L else 1L)
    },
    .object_summary = function(object) data.frame(Mean = 1, row.names = "mu"),
    .autocompute_brma = function(object) object,
    .summary_estimates_pair = function(...) list(estimates = list(), conditional = list()),
    .compute_selection_sensitivity_diagnostics = function(
        object, max_posterior_samples = 256L, latent_samples = 512L, seed = 1L,
        integration_control = NULL) {

      settings <- list(max_posterior_samples = max_posterior_samples,
                       latent_samples = latent_samples, seed = seed,
                       integration_control = integration_control)
      calls[[length(calls) + 1L]] <<- list(
        settings = settings, generation = object[["fit"]][["generation"]]
      )
      result <- data.frame(block = 1L, total_variation_median = .235)
      attr(result, "settings") <- settings
      result
    },
    .package = "RoBMA"
  )
  fit_model <- function(target = "condition", only_priors = FALSE) {

    bselmodel.mv(
      yi                        = c(.1, .2, .3),
      V                         = diag(rep(.04, 3)),
      random                    = NULL,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      prior_bias = BayesTools::prior_weightfunction("one-sided", .05,
        BayesTools::wf_fixed(c(1, .5)), model = BayesTools::selection_model(
          other_random_effects = target, known_sampling_variance = target)),
      only_priors               = only_priors,
      silent                    = TRUE
    )
  }

  expect_silent(object <- fit_model())
  expect_null(object[["selection_sensitivity_diagnostics"]])
  expect_silent(without_stored <- update(object, sample_extend = 10L, recompute = "drop"))
  expect_null(without_stored[["selection_sensitivity_diagnostics"]])
  expect_length(calls, 0L)
  requested <- expect_warning(selection_sensitivity_diagnostics(object), warning_text, fixed = TRUE)
  expect_s3_class(requested, "data.frame")
  expect_null(object[["selection_sensitivity_diagnostics"]])
  expect_identical(calls[[1L]][["settings"]], list(
    max_posterior_samples = 256L, latent_samples = 512L, seed = 1L,
    integration_control = NULL
  ))
  object <- expect_warning(add_selection_sensitivity_diagnostics(object), warning_text, fixed = TRUE)
  stored <- object[["selection_sensitivity_diagnostics"]]
  expect_s3_class(stored, "data.frame")
  expect_identical(expect_warning(selection_sensitivity_diagnostics(object), warning_text, fixed = TRUE), stored)
  expect_identical(expect_warning(selection_sensitivity_diagnostics(
    object, max_posterior_samples = 256, latent_samples = 512, seed = 1
  ), warning_text, fixed = TRUE), stored)
  summarized <- summary(object)
  expect_identical(summarized[["selection_sensitivity_diagnostics"]], stored)
  expect_silent(capture.output(print(object)))
  expect_silent(capture.output(print(summarized)))
  expect_length(calls, 2L)

  custom <- expect_warning(add_selection_sensitivity_diagnostics(
    object, max_posterior_samples = 10L, latent_samples = 8L, seed = 2L
  ), warning_text, fixed = TRUE)
  expect_length(calls, 3L)
  expect_identical(expect_warning(selection_sensitivity_diagnostics(custom), warning_text, fixed = TRUE),
                   custom[["selection_sensitivity_diagnostics"]])
  expect_silent(updated <- update(custom, sample_extend = 10L, recompute = "drop"))
  expect_length(calls, 4L)
  expect_identical(calls[[4L]][["settings"]], calls[[3L]][["settings"]])
  expect_identical(calls[[4L]][["generation"]], 2L)
  expect_warning(selection_sensitivity_diagnostics(updated), warning_text, fixed = TRUE)
  relabeled <- update(updated, slab = c("A", "B", "C"))
  expect_identical(relabeled[["selection_sensitivity_diagnostics"]], updated[["selection_sensitivity_diagnostics"]])
  expect_length(calls, 4L)

  expect_silent(integrated <- fit_model("integrate"))
  expect_null(integrated[["selection_sensitivity_diagnostics"]])
  expect_silent(unfitted <- fit_model(only_priors = TRUE))
  expect_null(unfitted[["selection_sensitivity_diagnostics"]])
  expect_length(calls, 4L)
  expect_error(selection_sensitivity_diagnostics(unfitted),
               "'object' must contain a fitted posterior.", fixed = TRUE)
  expect_error(selection_sensitivity_diagnostics(custom, seed = -1),
               "'seed' must be a single non-negative integer.", fixed = TRUE)
})


test_that("explicitly stored clustered diagnostics refresh without expanding automatic scope", {

  calls <- list()
  testthat::local_mocked_bindings(
    .fit = function(object, extend = FALSE) {
      list(has_posterior = TRUE, generation = if (extend) 2L else 1L)
    },
    .object_summary = function(object) data.frame(Mean = 1, row.names = "mu"),
    .object_coefficients = function(object) c(mu = 1),
    .autocompute_brma = function(object) object,
    .compute_selection_sensitivity_diagnostics = function(
        object, max_posterior_samples = 256L, latent_samples = 512L, seed = 1L,
        integration_control = NULL) {

      settings <- list(max_posterior_samples = max_posterior_samples,
                       latent_samples = latent_samples, seed = seed,
                       integration_control = integration_control)
      generation <- object[["fit"]][["generation"]]
      calls[[length(calls) + 1L]] <<- list(settings = settings, generation = generation)
      result <- data.frame(total_variation_median = 0, generation = generation)
      attr(result, "settings") <- settings
      result
    },
    .package = "RoBMA"
  )
  object <- bselmodel(
    yi = c(.1, .2, .3), sei = rep(.2, 3), cluster = c("A", "A", "B"),
    measure = "GEN", prior_unit_information_sd = 1,
    prior_bias = BayesTools::prior_weightfunction("one-sided", .05,
      BayesTools::wf_fixed(c(1, .5)), model = BayesTools::selection_model(
        other_random_effects = "condition", known_sampling_variance = "integrate")),
    silent = TRUE
  )
  expect_false(inherits(object, "brma.mv"))
  expect_null(object[["selection_sensitivity_diagnostics"]])
  without_stored <- update(object, sample_extend = 10L, recompute = "drop")
  expect_null(without_stored[["selection_sensitivity_diagnostics"]])
  expect_length(calls, 0L)

  stored <- add_selection_sensitivity_diagnostics(
    object, max_posterior_samples = 10L, latent_samples = 8L, seed = 2L
  )
  expect_length(calls, 1L)
  updated <- update(stored, sample_extend = 10L, recompute = "drop")
  expect_length(calls, 2L)
  expect_identical(calls[[2L]][["settings"]], calls[[1L]][["settings"]])
  expect_identical(calls[[2L]][["generation"]], 2L)
  expect_identical(updated[["selection_sensitivity_diagnostics"]][["generation"]], 2L)
  expect_identical(selection_sensitivity_diagnostics(updated),
                   updated[["selection_sensitivity_diagnostics"]])
  relabeled <- update(updated, slab = c("first", "second", "third"))
  expect_identical(relabeled[["selection_sensitivity_diagnostics"]],
                   updated[["selection_sensitivity_diagnostics"]])
  expect_length(calls, 2L)
})


test_that("unavailable refreshed comparisons preserve the fitted posterior", {

  fail <- FALSE
  testthat::local_mocked_bindings(
    .fit = function(object, extend = FALSE) list(has_posterior = TRUE),
    .object_summary = function(object) data.frame(Mean = 1, row.names = "mu"),
    .autocompute_brma = function(object) object,
    .summary_estimates_pair = function(...) list(estimates = list(), conditional = list()),
    .compute_selection_sensitivity_diagnostics = function(
        object, max_posterior_samples = 256L, latent_samples = 512L, seed = 1L,
        integration_control = NULL) {

      if (fail) stop("Conditional covariance is unsupported.", call. = FALSE)
      result <- data.frame(total_variation_median = 0)
      attr(result, "settings") <- list(max_posterior_samples = max_posterior_samples,
        latent_samples = latent_samples, seed = seed, integration_control = integration_control)
      result
    },
    .package = "RoBMA"
  )
  expect_silent(object <- bselmodel.mv(
    yi = c(.1, .2, .3), V = diag(rep(.04, 3)), random = NULL,
    measure = "GEN", prior_unit_information_sd = 1, silent = TRUE
  ))
  expect_null(object[["selection_sensitivity_diagnostics"]])
  stored <- add_selection_sensitivity_diagnostics(
    object, max_posterior_samples = 10L, latent_samples = 8L, seed = 2L
  )
  settings <- attr(stored[["selection_sensitivity_diagnostics"]], "settings")
  fail <- TRUE
  expect_silent(updated <- update(stored, sample_extend = 10L, recompute = "drop"))
  expect_true(updated[["fit"]][["has_posterior"]])
  expect_s3_class(updated[["selection_sensitivity_diagnostics"]], "error")
  expect_identical(attr(updated[["selection_sensitivity_diagnostics"]], "settings"), settings)
  expect_silent(capture.output(print(updated)))
  expect_error(selection_sensitivity_diagnostics(updated),
               "Conditional covariance is unsupported.", fixed = TRUE)
  expect_error(selection_sensitivity_diagnostics(object),
               "Conditional covariance is unsupported.", fixed = TRUE)
  expect_true(object[["fit"]][["has_posterior"]])
})


test_that("diagnostics integrate fresh sampling effects and preserve the RNG", {

  samples <- cbind(mu = c(.2, .3), tau = c(0, 0))
  testthat::local_mocked_bindings(
    .known_v_diagnostic_posterior_samples = function(...) {
      list(posterior_samples = samples, n_used = 2L, n_total = 2L)
    },
    .package = "RoBMA"
  )
  bias <- BayesTools::prior_weightfunction(
    side = "one-sided", steps = .05, weights = BayesTools::wf_fixed(c(1, .5)),
    model = BayesTools::selection_model(known_sampling_variance = "condition", group = "paper")
  )
  object <- bselmodel.mv(
    yi                        = c(.1, .2, .3),
    V                         = known_v_factor(rep(.04, 3), matrix(rep(.2, 3), ncol = 1)),
    data                      = data.frame(paper = c("A", "A", "A")),
    random                    = NULL,
    prior_bias                = bias,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    effect_direction          = "positive",
    only_priors               = TRUE
  )
  object[["fit"]] <- list(has_posterior = TRUE)
  set.seed(123)
  rng_before <- .Random.seed
  baseline <- .compute_selection_sensitivity_diagnostics(object, latent_samples = 32L)
  expect_identical(.Random.seed, rng_before)
  expect_true(all(is.finite(baseline[["total_variation_median"]])))
  expect_true(all(baseline[["total_variation_median"]] > 0))
  expect_identical(attr(baseline, "target")$fitted$known_sampling_variance, "condition")
  expect_identical(attr(baseline, "target")$reference$known_sampling_variance, "integrate")
  expect_true(attr(baseline, "target")$applicable)
  expect_identical(attr(baseline, "target")$publication_groups,
    .data_selection_model(object$data)$groups)

  # Fitted sampling effects must not shift the pre-selection latent law.
  samples <- cbind(samples, "sampling_z[1]" = c(-50, 50))
  shifted <- .compute_selection_sensitivity_diagnostics(object, latent_samples = 32L)
  expect_identical(shifted, baseline)
  expect_identical(.Random.seed, rng_before)

  original_compute <- .compute_selection_sensitivity_diagnostics
  testthat::local_mocked_bindings(
    .predict_joint_selection_gaussian_parts = function(...) {
      stop("Covariance plan unavailable.", call. = FALSE)
    },
    .package = "RoBMA"
  )
  expect_error(original_compute(object), "Covariance plan unavailable.", fixed = TRUE)
  expect_identical(.Random.seed, rng_before)
})


test_that("selection reweighting metrics retain their exact definitions", {

  constant <- .selection_sensitivity_weight_metrics(rep(-1000, 32L))
  expect_equal(constant[["ess_fraction"]], 1)
  expect_equal(constant[["ess_fraction_mcse"]], 0)
  expect_equal(constant[["total_variation"]], 0)
  expect_equal(constant[["total_variation_mcse"]], 0)
  expect_equal(constant[["log_weight_iqr"]], 0)

  log_weights <- log(c(0.2, 0.5, 1, 0.7, 0.1, 0.9))
  shifted     <- .selection_sensitivity_weight_metrics(log_weights - 800)
  direct      <- exp(log_weights)
  expected_ess <- mean(direct)^2 / mean(direct^2)
  expected_tv  <- 0.5 * mean(abs(direct / mean(direct) - 1))

  expect_equal(shifted[["ess_fraction"]], expected_ess)
  expect_equal(shifted[["total_variation"]], expected_tv)
  expect_equal(
    shifted[["log_weight_iqr"]],
    unname(stats::IQR(log_weights, type = 8))
  )
  expect_true(shifted[["ess_fraction_mcse"]] > 0)
  expect_true(shifted[["total_variation_mcse"]] > 0)
})


test_that("fresh context draws retain positive-semidefinite Gaussian covariance", {

  loading <- matrix(c(1, .2, 0, .5, 1, .4, 0, .4), 4L, 2L)
  covariance <- tcrossprod(loading)
  set.seed(445)
  draws <- .selection_sensitivity_covariance_draws(covariance, 100000L)
  expect_lt(max(abs(stats::cov(draws) - covariance)), .025)
  expect_identical(.selection_sensitivity_covariance_draws(matrix(0, 2L, 2L), 3L),
    matrix(0, 3L, 2L))
})


test_that("sensitivity preserves full events and has constant no-context weights", {

  latent <- c(-1.5, -.5, 0, .5, 1.5)
  testthat::local_mocked_bindings(
    .selection_sensitivity_covariance_draws = function(covariance, latent_samples) {
      cbind(latent, latent / 2)
    }, .package = "RoBMA"
  )
  plan <- c(set_selection_likelihood_control(), list(designs = list()))
  for (rule in c("product", "best")) {
    prior <- BayesTools::prior_weightfunction("one-sided", .5, BayesTools::wf_fixed(c(1, .2)),
      model = BayesTools::selection_model(weight_rule = rule, group = "paper"))
    selection <- .selection_spec(list(outcome = list(bias = prior)), c(0, 0), c(1, 1), "positive")
    selection$omega <- matrix(c(1, .2), 1L)
    selection$kernel_mode <- SELKERNEL_STEP
    selection$vector_rule <- if (rule == "product") 0L else 1L
    selection$alpha <- 0
    selection$phack_kind <- 0L
    selection$use_normal <- FALSE
    arguments <- list(means = matrix(c(0, 0), 1L),
      covariance = array(diag(2), c(1L, 2L, 2L)),
      context_covariance = array(tcrossprod(c(1, .5)), c(1L, 2L, 2L)),
      selection_sei = c(1, 1), selection_context = selection,
      normalization_units = if (rule == "product") list(1L, 2L) else list(1:2),
      row_blocks = list(1:2), execution_plan = plan, latent_samples = 5L)
    out <- do.call(.selection_sensitivity_run, arguments)
    below <- cbind(stats::pnorm(-latent), stats::pnorm(-latent / 2))
    weights <- if (rule == "product") apply(1 - .8 * below, 1L, prod) else
      1 - .8 * apply(below, 1L, prod)
    expect_equal(out$total_variation_median, .5 * mean(abs(weights / mean(weights) - 1)), tolerance = 1e-12)
    expect_equal(out$ess_fraction_median, mean(weights)^2 / mean(weights^2), tolerance = 1e-12)
    arguments$context_covariance[] <- 0
    constant <- do.call(.selection_sensitivity_run, arguments)
    expect_identical(constant$total_variation_median, 0)
    expect_identical(constant$ess_fraction_median, 1)
    expect_identical(constant$integration_relative_mcse_max, 0)
  }
})


test_that("selection sensitivity summaries expose posterior and MC uncertainty", {

  metrics <- list(
    ess_fraction = matrix(c(0.2, 0.4, 0.6, 0.8), 2L, 2L),
    ess_fraction_mcse = matrix(c(0.01, 0.02, 0.03, 0.04), 2L, 2L),
    total_variation = matrix(c(0.1, 0.2, 0.3, 0.4), 2L, 2L),
    total_variation_mcse = matrix(c(0.02, 0.03, 0.04, 0.05), 2L, 2L),
    log_weight_iqr = matrix(c(1, 2, 3, 4), 2L, 2L)
  )
  out <- .selection_sensitivity_summarize(
    metrics,
    list(1:2, 3:5)
  )

  expect_s3_class(out, "data.frame")
  expect_identical(out[["rows"]], c("1,2", "3,4,5"))
  expect_identical(out[["n_estimates"]], c(2L, 3L))
  expect_equal(out[["ess_fraction_median"]], c(0.3, 0.7))
  expect_equal(out[["ess_fraction_mcse_max"]], c(0.02, 0.04))
  expect_equal(out[["total_variation_mcse_max"]], c(0.03, 0.05))
})


test_that("selection sensitivity notifications use posterior-median TV", {

  expect_silent(.selection_sensitivity_notify(c(0.01, 0.05)))
  expect_message(
    .selection_sensitivity_notify(c(0.01, 0.064)),
    paste0(
      "Selection conditioning sensitivity: the largest unit-level ",
      "posterior-median total variation distance was 6.4%. Inspect ",
      "'selection_sensitivity_diagnostics()' for unit-level ",
      "results and the fitted and reference specifications."
    ),
    fixed = TRUE
  )
  expect_warning(
    .selection_sensitivity_notify(c(0.01, 0.10)),
    paste0(
      "Selection conditioning sensitivity was substantial: the largest unit-level ",
      "posterior-median total variation distance was 10.0%. ",
      "Retaining and integrating the specified contexts define different reporting models. ",
      "Inspect 'selection_sensitivity_diagnostics()' for the fitted and reference specifications."
    ),
    fixed = TRUE
  )
})
