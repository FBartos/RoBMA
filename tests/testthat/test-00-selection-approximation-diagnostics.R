test_that("approximate fits cache diagnostics and print their notification", {

  calls <- list()
  warning_text <- paste0(
    "The approximate selection likelihood showed substantial ",
    "latent-distribution reweighting: the largest block-level ",
    "posterior-median total variation distance was 23.5%. The approximate ",
    "and exact likelihoods may yield different inference. Consider ",
    "refitting with 'selection_likelihood' set to 'exact'."
  )
  testthat::local_mocked_bindings(
    .fit = function(object, extend = FALSE) {
      list(has_posterior = TRUE, generation = if (extend) 2L else 1L)
    },
    .object_summary = function(object) data.frame(Mean = 1, row.names = "mu"),
    .autocompute_brma = function(object) object,
    .summary_estimates_pair = function(...) list(estimates = list(), conditional = list()),
    .compute_selection_approximation_diagnostics = function(
        object, max_posterior_samples = 256L, latent_samples = 512L, seed = 1L) {

      settings <- list(max_posterior_samples = max_posterior_samples,
                       latent_samples = latent_samples, seed = seed)
      calls[[length(calls) + 1L]] <<- list(
        settings = settings, generation = object[["fit"]][["generation"]]
      )
      result <- data.frame(block = 1L, total_variation_median = .235)
      attr(result, "settings") <- settings
      result
    },
    .package = "RoBMA"
  )
  fit_model <- function(target = "approximate", only_priors = FALSE) {

    bselmodel.mv(
      yi                        = c(.1, .2, .3),
      V                         = diag(rep(.04, 3)),
      random                    = NULL,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      selection_likelihood      = target,
      only_priors               = only_priors,
      silent                    = TRUE
    )
  }

  object <- expect_warning(fit_model(), warning_text, fixed = TRUE)
  stored <- object[["selection_approximation_diagnostics"]]
  expect_s3_class(stored, "data.frame")
  expect_length(calls, 1L)
  expect_identical(
    expect_warning(selection_approximation_diagnostics(object), warning_text, fixed = TRUE),
    stored
  )
  expect_identical(
    expect_warning(selection_approximation_diagnostics(
      object, max_posterior_samples = 256, latent_samples = 512, seed = 1
    ), warning_text, fixed = TRUE),
    stored
  )
  summarized <- summary(object)
  expect_identical(summarized[["selection_approximation_diagnostics"]], stored)
  expect_warning(capture.output(print(object)), warning_text, fixed = TRUE)
  expect_warning(capture.output(print(summarized)), warning_text, fixed = TRUE)
  expect_length(calls, 1L)

  custom <- expect_warning(add_selection_approximation_diagnostics(
    object, max_posterior_samples = 10L, latent_samples = 8L, seed = 2L
  ), warning_text, fixed = TRUE)
  expect_length(calls, 2L)
  expect_identical(
    expect_warning(selection_approximation_diagnostics(custom), warning_text, fixed = TRUE),
    custom[["selection_approximation_diagnostics"]]
  )
  expect_length(calls, 2L)
  updated <- expect_warning(
    update(custom, sample_extend = 10L, recompute = "drop"),
    warning_text, fixed = TRUE
  )
  expect_length(calls, 3L)
  expect_identical(calls[[3L]][["settings"]], calls[[2L]][["settings"]])
  expect_identical(calls[[3L]][["generation"]], 2L)
  relabeled <- update(updated, slab = c("A", "B", "C"))
  expect_identical(relabeled[["selection_approximation_diagnostics"]],
                   updated[["selection_approximation_diagnostics"]])
  expect_length(calls, 3L)

  expect_silent(exact <- fit_model("exact"))
  expect_null(exact[["selection_approximation_diagnostics"]])
  expect_silent(unfitted <- fit_model(only_priors = TRUE))
  expect_null(unfitted[["selection_approximation_diagnostics"]])
  expect_length(calls, 3L)
  expect_error(selection_approximation_diagnostics(unfitted),
               "'object' must contain a fitted posterior.", fixed = TRUE)
  expect_error(selection_approximation_diagnostics(custom, seed = -1),
               "'seed' must be a single non-negative integer.", fixed = TRUE)
})


test_that("unavailable automatic diagnostics preserve the fitted posterior", {

  testthat::local_mocked_bindings(
    .fit = function(object, extend = FALSE) list(has_posterior = TRUE),
    .object_summary = function(object) data.frame(Mean = 1, row.names = "mu"),
    .autocompute_brma = function(object) object,
    .summary_estimates_pair = function(...) list(estimates = list(), conditional = list()),
    .compute_selection_approximation_diagnostics = function(...) {
      stop("Conditional covariance is unsupported.", call. = FALSE)
    },
    .package = "RoBMA"
  )
  warning_text <- paste0(
    "Selection-approximation diagnostics are unavailable: ",
    "Conditional covariance is unsupported."
  )
  object <- expect_warning(bselmodel.mv(
    yi = c(.1, .2, .3), V = diag(rep(.04, 3)), random = NULL,
    measure = "GEN", prior_unit_information_sd = 1,
    selection_likelihood = "approximate", silent = TRUE
  ), warning_text, fixed = TRUE)
  expect_true(object[["fit"]][["has_posterior"]])
  expect_s3_class(object[["selection_approximation_diagnostics"]], "error")
  expect_warning(capture.output(print(object)), warning_text, fixed = TRUE)
  expect_error(selection_approximation_diagnostics(object),
               "Conditional covariance is unsupported.", fixed = TRUE)
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
    side = "one-sided", steps = .05, weights = BayesTools::wf_fixed(c(1, .5))
  )
  object <- bselmodel.mv(
    yi                        = c(.1, .2, .3),
    V                         = known_v_factor(rep(.04, 3), matrix(rep(.2, 3), ncol = 1)),
    random                    = NULL,
    prior_bias                = bias,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    selection_likelihood      = "approximate",
    effect_direction          = "positive",
    only_priors               = TRUE
  )
  object[["fit"]] <- list(has_posterior = TRUE)
  set.seed(123)
  rng_before <- .Random.seed
  baseline <- .compute_selection_approximation_diagnostics(object, latent_samples = 32L)
  expect_identical(.Random.seed, rng_before)
  expect_true(all(is.finite(baseline[["total_variation_median"]])))
  expect_true(all(baseline[["total_variation_median"]] > 0))

  # Fitted sampling effects must not shift the pre-selection latent law.
  samples <- cbind(samples, "sampling_z[1]" = c(-50, 50))
  shifted <- .compute_selection_approximation_diagnostics(object, latent_samples = 32L)
  expect_identical(shifted, baseline)
  expect_identical(.Random.seed, rng_before)

  original_compute <- .compute_selection_approximation_diagnostics
  testthat::local_mocked_bindings(
    .known_v_marginal_factor_plan = function(...) {
      stop("Covariance plan unavailable.", call. = FALSE)
    },
    .package = "RoBMA"
  )
  expect_error(original_compute(object), "Covariance plan unavailable.", fixed = TRUE)
  expect_identical(.Random.seed, rng_before)
})


test_that("selection reweighting metrics retain their exact definitions", {

  constant <- .selection_approximation_weight_metrics(rep(-1000, 32L))
  expect_equal(constant[["ess_fraction"]], 1)
  expect_equal(constant[["ess_fraction_mcse"]], 0)
  expect_equal(constant[["total_variation"]], 0)
  expect_equal(constant[["total_variation_mcse"]], 0)
  expect_equal(constant[["log_weight_iqr"]], 0)

  log_weights <- log(c(0.2, 0.5, 1, 0.7, 0.1, 0.9))
  shifted     <- .selection_approximation_weight_metrics(log_weights - 800)
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


test_that("selection latent draws preserve every covariance factor type", {

  latent_samples <- 100000L
  K              <- 4L
  X <- matrix(
    c(1, 0.2, 1, -0.4, 0, 0.6, 0.5, 0),
    nrow = K,
    ncol = 2L,
    byrow = TRUE
  )
  coefficient_factor <- matrix(c(1, 0.2, 0, 0.7), 2L, 2L)
  group_map           <- c(1L, 1L, 2L, 2L)
  basis               <- X %*% coefficient_factor

  plans <- list(
    group = list(
      type                  = "group",
      model_matrix          = X,
      group_map             = group_map,
      coefficient_structure = "dense"
    ),
    row_group = list(
      type                  = "row_group",
      model_matrix          = X,
      group_map             = group_map,
      coefficient_structure = "dense"
    ),
    known_group = list(
      type                  = "known_group",
      model_matrix          = X,
      group_map             = group_map,
      coefficient_structure = "dense",
      group_covariance      = matrix(c(1, 0.3, 0.3, 1.4), 2L, 2L)
    )
  )
  states <- list(
    group = list(coefficient_factor = coefficient_factor),
    row_group = list(
      coefficient_factor = coefficient_factor,
      row_scale           = c(0.5, 1, 1.5, 2)
    ),
    known_group = list(coefficient_factor = coefficient_factor)
  )
  expected <- list(
    group = tcrossprod(basis) * outer(group_map, group_map, "=="),
    row_group = {
      scaled_basis <- basis * states[["row_group"]][["row_scale"]]
      tcrossprod(scaled_basis) * outer(group_map, group_map, "==")
    },
    known_group = tcrossprod(basis) *
      plans[["known_group"]][["group_covariance"]][
        group_map,
        group_map,
        drop = FALSE
      ]
  )

  for (type in names(plans)) {
    set.seed(440 + match(type, names(plans)))
    draws <- .selection_approximation_factor_draws(
      plan           = plans[[type]],
      state          = states[[type]],
      rows           = seq_len(K),
      latent_samples = latent_samples,
      K              = K
    )
    expect_lt(max(abs(stats::cov(draws) - expected[[type]])), 0.025)
  }

  dense_covariance <- tcrossprod(matrix(
    c(1, 0.2, 0, 0.5, 1, 0.4, 0, 0.4),
    nrow = K,
    ncol = 2L
  ))
  set.seed(445)
  dense_draws <- .selection_approximation_factor_draws(
    plan           = list(type = "dense"),
    state          = list(covariance = dense_covariance),
    rows           = seq_len(K),
    latent_samples = latent_samples,
    K              = K
  )
  expect_lt(max(abs(stats::cov(dense_draws) - dense_covariance)), 0.025)
})


test_that("selection approximation summaries expose posterior and MC uncertainty", {

  metrics <- list(
    ess_fraction = matrix(c(0.2, 0.4, 0.6, 0.8), 2L, 2L),
    ess_fraction_mcse = matrix(c(0.01, 0.02, 0.03, 0.04), 2L, 2L),
    total_variation = matrix(c(0.1, 0.2, 0.3, 0.4), 2L, 2L),
    total_variation_mcse = matrix(c(0.02, 0.03, 0.04, 0.05), 2L, 2L),
    log_weight_iqr = matrix(c(1, 2, 3, 4), 2L, 2L)
  )
  out <- .selection_approximation_summarize(
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


test_that("selection approximation notifications use posterior-median TV", {

  expect_silent(.selection_approximation_notify(c(0.01, 0.05)))
  expect_message(
    .selection_approximation_notify(c(0.01, 0.064)),
    paste0(
      "Approximate selection-likelihood diagnostic: the largest block-level ",
      "posterior-median total variation distance was 6.4%. Inspect ",
      "'selection_approximation_diagnostics()' for block-level results."
    ),
    fixed = TRUE
  )
  expect_warning(
    .selection_approximation_notify(c(0.01, 0.10)),
    paste0(
      "The approximate selection likelihood showed substantial ",
      "latent-distribution reweighting: the largest block-level ",
      "posterior-median total variation distance was 10.0%. The approximate ",
      "and exact likelihoods may yield different inference. Consider ",
      "refitting with 'selection_likelihood' set to 'exact'."
    ),
    fixed = TRUE
  )
})
