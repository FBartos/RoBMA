.bridge_repetition_test_fit <- function(marglik) {

  data <- list(outcome = data.frame(yi = c(.1, -.2), sei = c(.2, .3)))
  attr(data, "outcome_type") <- "norm"
  structure(list(data = data, marglik = marglik), class = "brma")
}

.bridge_repetition_test_raw <- function(logml) {

  structure(
    list(logml = logml, niter = rep(4L, length(logml)), method = "normal",
         repetitions = length(logml)),
    class = if (length(logml) > 1L) "bridge_list" else "bridge"
  )
}

.bridge_repetition_test_stored <- function(upstream, nonfinite = "error") {

  BayesTools:::.bt_marglik_from_upstream(
    upstream, maxiter = 100, nonfinite = nonfinite,
    chain_metadata = list(count = 2L, draws = 512L, draws_per_chain = c(256L, 256L))
  )
}

test_that("add_marglik translates only thrown typed repetition failures", {

  raw <- .bridge_repetition_test_raw(c(-2, NA_real_, Inf))
  original <- tryCatch(.bridge_repetition_test_stored(raw), error = identity)
  expect_s3_class(original, "BayesTools_marglik_repetition_failure")
  fit <- .bridge_repetition_test_fit(NULL)
  signal <- function(condition) {
    testthat::local_mocked_bindings(.marglik = function(...) stop(condition),
      .package = "RoBMA")
    tryCatch(add_marglik(fit), error = identity)
  }
  actual <- signal(original)
  expect_identical(class(actual), c("RoBMA_marglik_repetition_failure",
    "RoBMA_marglik_error", "error", "condition"))
  expect_identical(conditionMessage(actual), paste0(
    "Bridge sampling was rejected by diagnostics: repetition(s) 2, 3 returned ",
    "a non-finite natural-log marginal likelihood. Increase 'maxiter' or fit ",
    "more posterior draws (for example with a larger 'sample' or update())."
  ))
  expect_null(actual$call)
  expect_identical(actual$parent, original)
  expect_identical(actual$parent$repetitions$repetition, c(2L, 3L))
  expect_identical(actual$parent$repetitions, original$repetitions)
  expect_identical(conditionMessage(actual$parent), conditionMessage(original))
  expect_identical(class(actual$parent), class(original))

  unrelated <- list(
    errorCondition("Unrelated typed input failure.",
      class = "BayesTools_marglik_input", diagnostic = list(nodes = "mu")),
    errorCondition("Refit required.", class = "BayesTools_marglik_refit_required"),
    errorCondition("Required monitor is missing.", class = "BayesTools_marglik_missing_monitor"),
    simpleError("Unrelated callback or programming failure."),
    simpleError(conditionMessage(original))
  )
  for (condition in unrelated) {
    expect_identical(signal(condition), condition)
  }
})

test_that("add_marglik warns about unused dots after existing argument validation", {

  stored <- .bridge_repetition_test_stored(.bridge_repetition_test_raw(-2))
  fit <- .bridge_repetition_test_fit(NULL)
  testthat::local_mocked_bindings(.marglik = function(...) stored,
    .package = "RoBMA")
  expect_warning(result <- add_marglik(fit, unknown = TRUE),
    "Unused argument in add_marglik(): 'unknown'", fixed = TRUE)
  expect_identical(result$marglik, stored)
  expect_warning(add_marglik(fit, nonfinite = "drop"),
    "Unused argument in add_marglik(): 'nonfinite'", fixed = TRUE)
  # Supply an unnamed entry to dots after every ordinary formal argument.
  expect_warning(do.call(add_marglik, c(list(object = fit, parallel = FALSE,
    cores = 1L, repetitions = 1L, method = "normal", maxiter = 100L,
    silent = TRUE), list(TRUE))),
    "Unused argument in add_marglik(): '<unnamed>'", fixed = TRUE)
  expect_warning(add_marglik(fit, parallel = FALSE, cores = 1L,
    repetitions = 3L, method = "warp3", maxiter = 100L, silent = FALSE), NA)
  expect_warning(expect_error(add_marglik(fit, method = "bad", unknown = TRUE),
    "'arg' should be one of", fixed = TRUE), NA)
  expect_warning(expect_error(add_marglik(fit, cores = 0L, unknown = TRUE),
    "cores", fixed = TRUE), NA)
  product <- fit
  class(product) <- c("RoBMA", "brma")
  expect_warning(expect_error(add_marglik(product, unknown = TRUE),
    "product-space", fixed = TRUE), NA)
})

test_that("public posterior model probabilities normalize finite prior weights safely", {

  controls <- list(c(1e308, 1e308), c(1e308, 5e307), c(1e-300, 2e-300),
    c(0, 1), c(1, 3), c(1e100, 3e100))
  probabilities <- list(c(.5, .5), c(2 / 3, 1 / 3), c(1 / 3, 2 / 3),
    c(0, 1), c(.25, .75), c(.25, .75))
  for (repetitions in c(1L, 3L)) {
    fit <- .bridge_repetition_test_fit(.bridge_repetition_test_stored(
      .bridge_repetition_test_raw(rep(0, repetitions))))
    for (i in seq_along(controls)) {
      actual <- post_prob(fit, fit, prior_prob = controls[[i]],
        model_names = c("A", "B"))
      expect_equal(as.numeric(actual), rep(probabilities[[i]], each = repetitions),
        tolerance = 1e-14)
      expect_true(all(is.finite(actual)))
      if (repetitions == 1L) {
        expect_identical(names(actual), c("A", "B"))
      } else {
        expect_identical(colnames(actual), c("A", "B"))
        expect_equal(rowSums(actual), rep(1, repetitions), tolerance = 1e-14)
      }
    }
    for (invalid in list(c(0, 0), c(1, Inf), 1)) {
      expected <- if (length(invalid) == 1L) {
        "'prior_prob' must have the same length as the number of models."
      } else if (any(!is.finite(invalid))) {
        "'prior_prob' must contain only finite values."
      } else "'prior_prob' must contain at least one positive value."
      error <- tryCatch(post_prob(fit, fit, prior_prob = invalid), error = identity)
      expect_s3_class(error, "error")
      expect_identical(conditionMessage(error), expected)
    }
  }
  expect_identical(names(.post_prob_normalize_prior(c(A = 1e308, B = 1e308), 2L)),
    c("A", "B"))
})

test_that("real bridge repetitions survive extraction and model comparisons", {

  withr::local_seed(360)
  draws <- matrix(stats::qnorm(stats::ppoints(512)), ncol = 1L,
                  dimnames = list(NULL, "theta"))
  raw <- bridgesampling::bridge_sampler(
    draws, log_posterior = function(samples, data) {
      stats::dnorm(samples[["theta"]], log = TRUE) + data
    },
    data = -3, lb = c(theta = -Inf), ub = c(theta = Inf),
    repetitions = 3, cores = 1, use_neff = FALSE, silent = TRUE
  )
  expect_s3_class(raw, "bridge_list")
  stored <- .bridge_repetition_test_stored(raw)
  attr(stored, "RoBMA_target") <- list(reported_target = "full joint fitted likelihood")
  fit <- .bridge_repetition_test_fit(stored)
  extracted <- bridge_sampler(fit)
  expect_s3_class(extracted, "bridge_list")
  expect_equal(extracted$logml, raw$logml)
  expect_identical(attr(extracted, "RoBMA_target"), attr(stored, "RoBMA_target"))
  expect_identical(logml(fit), stored$logml)
  expect_length(logml(fit), 1L)

  other <- raw
  other$logml <- raw$logml + c(.1, -.2, .4)
  other_fit <- .bridge_repetition_test_fit(.bridge_repetition_test_stored(other))
  expected <- bridgesampling::bf(raw, other)
  actual <- bf(fit, other_fit)
  expect_equal(actual$bf, expected$bf)
  expect_equal(actual$bf_median_based, expected$bf_median_based)
  expect_equal(post_prob(fit, other_fit, model_names = c("one", "two")),
               bridgesampling::post_prob(raw, other, model_names = c("one", "two")))
  expect_identical(fit$marglik, stored)
})

test_that("repeated Bayes factors use separate medians before upstream recycling", {

  first_raw <- .bridge_repetition_test_raw(c(-3, 0, 2))
  second_raw <- .bridge_repetition_test_raw(c(-2, 1))
  first <- .bridge_repetition_test_fit(.bridge_repetition_test_stored(first_raw))
  second <- .bridge_repetition_test_fit(.bridge_repetition_test_stored(second_raw))

  for (log in c(FALSE, TRUE)) {
    expect_warning(expected <- bridgesampling::bf(first_raw, second_raw, log = log),
                   "Not all objects provide 3 logmls. Some values are recycled.", fixed = TRUE)
    expect_warning(actual <- bf(first, second, log = log),
                   "Not all objects provide 3 logmls. Some values are recycled.", fixed = TRUE)
    expect_s3_class(actual, "bf_bridge_list")
    expect_equal(actual$bf, expected$bf)
    expect_equal(actual$bf_median_based, expected$bf_median_based)
    expect_output(print(actual), "based on medians")

    tabular <- as.data.frame(actual)
    expect_identical(data.frame(actual), tabular)
    expect_named(tabular, c("component", "parameter", "value"))
    expect_equal(tabular$value[tabular$component == "repetitions"], actual$bf)
    expect_equal(tabular$value[tabular$component == "summary"],
                 c(actual$bf_median_based, range(actual$bf), stats::IQR(actual$bf)))
  }
})

test_that("a scalar first model still retains all repeated posterior comparisons", {

  repeated_raw <- .bridge_repetition_test_raw(c(-2, -1, 0))
  scalar_raw <- .bridge_repetition_test_raw(-1.5)
  repeated <- .bridge_repetition_test_fit(.bridge_repetition_test_stored(repeated_raw))
  scalar <- .bridge_repetition_test_fit(.bridge_repetition_test_stored(scalar_raw))
  expect_warning(reference <- bridgesampling::post_prob(
    repeated_raw, scalar_raw, prior_prob = c(.25, .75), model_names = c("repeat", "scalar")
  ), "Some values are recycled.", fixed = TRUE)
  expect_warning(actual <- post_prob(
    scalar, repeated, prior_prob = c(3, 1), model_names = c("scalar", "repeat")
  ), "Some values are recycled.", fixed = TRUE)
  expect_equal(actual, reference[, c("scalar", "repeat")])
  expect_equal(rowSums(actual), rep(1, 3))

  expect_warning(actual_bf <- bf(scalar, repeated, log = TRUE),
                 "Some values are recycled.", fixed = TRUE)
  expect_equal(actual_bf$bf, scalar_raw$logml - repeated_raw$logml)
  expect_equal(actual_bf$bf_median_based, scalar_raw$logml - stats::median(repeated_raw$logml))
  expect_warning(aliased <- bayes_factor(scalar, repeated, log = TRUE),
                 "Some values are recycled.", fixed = TRUE)
  expect_equal(aliased$bf, actual_bf$bf)
  expect_identical(attr(aliased, "model_names"), c("scalar", "repeated"))
})

test_that("standardized comparisons refuse nonfinite repeats and legacy propagates NA", {

  failed_raw <- .bridge_repetition_test_raw(c(-2, NA_real_, -1))
  # A finite estimate that exhausted maxiter remains included under "drop".
  failed_raw$niter[1L] <- 101L
  finite_raw <- .bridge_repetition_test_raw(c(0, 1, 2))
  legacy_failed <- .bridge_repetition_test_fit(failed_raw)
  legacy_finite <- .bridge_repetition_test_fit(finite_raw)
  expect_identical(bridge_sampler(legacy_finite), finite_raw)
  expect_identical(logml(legacy_finite), 1)
  expected <- bridgesampling::bf(failed_raw, finite_raw)
  actual <- bf(legacy_failed, legacy_finite)
  expect_equal(actual$bf, expected$bf)
  expect_equal(actual$bf_median_based, expected$bf_median_based)
  expect_warning(print(actual), "estimate(s) are NAs.", fixed = TRUE)
  expect_warning(reference <- bridgesampling::post_prob(failed_raw, finite_raw,
                   model_names = c("failed", "finite")),
                 "NAs in logml values.", fixed = TRUE)
  expect_warning(probabilities <- post_prob(legacy_failed, legacy_finite,
                   model_names = c("failed", "finite")),
                 "NAs in logml values.", fixed = TRUE)
  expect_equal(probabilities, reference)
  expect_true(all(is.na(probabilities[2, ])))

  expect_warning(dropped <- .bridge_repetition_test_stored(failed_raw, "drop"),
                 "Dropped non-finite", fixed = TRUE)
  dropped_fit <- .bridge_repetition_test_fit(dropped)
  comparison_fit <- .bridge_repetition_test_fit(.bridge_repetition_test_raw(c(0, 1)))
  expect_error(bf(dropped_fit, comparison_fit, log = TRUE),
    "Stored log marginal likelihoods must contain only finite values.", fixed = TRUE)
  expect_error(post_prob(dropped_fit, comparison_fit),
    "Stored log marginal likelihoods must contain only finite values.", fixed = TRUE)
  expect_identical(logml(dropped_fit), dropped$logml)
  expect_equal(dropped_fit$marglik$repetitions$logml, c(-2, NA, -1))
  dropped_fit$marglik$aggregation$nonfinite_policy <- "error"
  expect_error(bf(dropped_fit, comparison_fit), "must contain only finite values", fixed = TRUE)

  expect_warning(one_left <- .bridge_repetition_test_stored(
    .bridge_repetition_test_raw(c(NA_real_, -1)), "drop"
  ), "Dropped non-finite", fixed = TRUE)
  one_left <- .bridge_repetition_test_fit(one_left)
  scalar <- .bridge_repetition_test_fit(.bridge_repetition_test_raw(0))
  expect_error(bf(one_left, scalar), "must contain only finite values", fixed = TRUE)
  expect_error(post_prob(one_left, scalar), "must contain only finite values", fixed = TRUE)

  finite_raw$niter[1L] <- 101L
  over_limit <- .bridge_repetition_test_fit(.bridge_repetition_test_stored(finite_raw))
  expect_false(over_limit$marglik$repetitions$success[1L])
  expect_identical(logml(over_limit), 1)
  expect_equal(bf(over_limit, legacy_finite, log = TRUE)$bf, c(0, 0, 0))
  expect_equal(post_prob(over_limit, legacy_finite), matrix(.5, 3L, 2L,
    dimnames = list(NULL, c("over_limit", "legacy_finite"))))
})

test_that("single and exact marginal likelihood comparisons retain their contracts", {

  first <- .bridge_repetition_test_fit(.bridge_repetition_test_stored(
    .bridge_repetition_test_raw(-2)
  ))
  exact <- BayesTools:::.bt_marglik_exact_result(
    -1, list(count = 1L, draws = 1L, draws_per_chain = 1L)
  )
  second <- .bridge_repetition_test_fit(exact)
  expect_identical(class(bf(first, second)), "bf_default")
  expect_equal(bf(first, second)$bf, exp(-1))
  expect_length(post_prob(first, second), 2L)
  expect_null(dim(post_prob(first, second)))
  expect_error(bridge_sampler(second), class = "RoBMA_exact_marglik_no_bridge")
  expect_identical(logml(second), -1)

  repeated <- .bridge_repetition_test_fit(.bridge_repetition_test_raw(c(-2, -1)))
  expect_error(post_prob(first, repeated, prior_prob = c(1, -1)), "prior_prob", fixed = TRUE)
  expect_error(post_prob(first, repeated, prior_prob = c(0, 0)), "positive value", fixed = TRUE)
  expect_error(post_prob(first, repeated, model_names = "one"), "same length", fixed = TRUE)
  different <- repeated
  different$data$outcome$yi[1L] <- 9
  expect_error(bf(first, different), "same outcome data", fixed = TRUE)
  expect_error(post_prob(first, different), "same outcome data", fixed = TRUE)
})

test_that("finite model weights retain likelihood-compensated small probabilities", {

  weights <- c(1e-300, 1e300)
  # Actual binary input references evaluated at 100-digit precision.
  reference <- c(9.7213220602525587e-9, 1.3156378855174478e-9,
    9.7213220602525587e-9)
  for (repetitions in c(1L, 3L)) {
    first <- .bridge_repetition_test_fit(.bridge_repetition_test_stored(
      .bridge_repetition_test_raw(c(1400, 1402, 1400)[seq_len(repetitions)])))
    second <- .bridge_repetition_test_fit(.bridge_repetition_test_stored(
      .bridge_repetition_test_raw(rep(0, repetitions))))
    actual <- post_prob(first, second, prior_prob = weights,
      model_names = c("same", "same"))
    small <- if (repetitions == 1L) actual[2L] else actual[, 2L]
    expect_true(all(abs(small / reference[seq_len(repetitions)] - 1) <= 1e-10))
    expect_equal(as.numeric(actual), c(1 - reference[seq_len(repetitions)],
      reference[seq_len(repetitions)]), tolerance = 1e-14)
    if (repetitions == 1L) {
      expect_identical(names(actual), c("same", "same"))
      expect_equal(sum(actual), 1, tolerance = 1e-14)
    } else {
      expect_identical(dim(actual), c(3L, 2L))
      expect_identical(colnames(actual), c("same", "same.1"))
      expect_equal(rowSums(actual), rep(1, repetitions), tolerance = 1e-14)
    }
  }
})

test_that("common extreme log likelihoods preserve normalized model probabilities", {

  for (repetitions in c(1L, 3L)) {
    fit <- .bridge_repetition_test_fit(.bridge_repetition_test_stored(
      .bridge_repetition_test_raw(rep(1e308, repetitions))))
    expect_warning(actual <- post_prob(fit, fit, prior_prob = c(1, 3)), NA)
    expect_equal(as.numeric(actual), rep(c(.25, .75), each = repetitions),
      tolerance = 1e-14)
    if (repetitions == 1L) {
      expect_identical(names(actual), c("fit", "fit"))
      expect_equal(sum(actual), 1, tolerance = 1e-14)
    } else {
      expect_identical(colnames(actual), c("fit", "fit.1"))
      expect_equal(rowSums(actual), rep(1, repetitions), tolerance = 1e-14)
    }
  }
})

test_that("finite log model probabilities retain zero and extreme input contracts", {

  expect_warning(integer <- .post_prob_from_logml(
    c(-2147483647L, 2147483647L)), NA)
  expect_equal(integer, c(0, 1), tolerance = 1e-14)
  expect_equal(.post_prob_from_logml(c(1e308, 0), c(0, 1)), c(0, 1))
  expect_equal(.post_prob_from_logml(c(-1.7e308, 1.7e308)), c(0, 1))
  expect_equal(.post_prob_from_logml(c(0, 0), c(1e-300, 1e300)), c(0, 1))
  expect_equal(.post_prob_from_logml(c(0, 0)), c(.5, .5), tolerance = 1e-14)
  smallest <- .Machine$double.xmin * .Machine$double.eps
  compensated <- .post_prob_from_logml(c(740, 0), c(smallest, 1))
  expect_lte(abs(compensated[1L] / .011657587869975697 - 1), 1e-10)
  expect_equal(sum(compensated), 1, tolerance = 1e-14)
  expect_identical(names(.post_prob_from_logml(c(0, 0), c(1, 3),
    model_names = c("A", "A"))), c("A", "A"))
})

test_that("legacy nonfinite repeated comparisons delegate upstream contracts", {

  capture <- function(expression) {
    warnings <- character()
    value <- tryCatch(withCallingHandlers(expression, warning = function(warning) {
      warnings <<- c(warnings, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }), error = identity)
    list(value = value, warnings = warnings)
  }
  weights <- c(1, 3)
  # Independent legacy normalization, deliberately not the repaired helper.
  normalized <- weights / max(weights)
  normalized <- normalized / sum(normalized)
  for (bad in c(NA_real_, NaN, Inf, -Inf)) {
    for (position in c(1L, 2L)) {
      for (model in c(1L, 2L)) {
        values <- list(c(-2, -1, 0), c(0, 1, 2))
        values[[model]][position] <- bad
        fits <- lapply(values, function(value) {
          .bridge_repetition_test_fit(.bridge_repetition_test_raw(value))
        })
        actual <- capture(post_prob(fits[[1L]], fits[[2L]], prior_prob = weights,
          model_names = c("A", "A")))
        expected <- capture(t(vapply(seq_len(3L), function(row) {
          bridgesampling::post_prob(vapply(values, "[[", numeric(1), row),
            prior_prob = normalized, model_names = c("A", "A"))
        }, numeric(2L))))
        expect_identical(actual$warnings, expected$warnings)
        expect_identical(dim(actual$value), dim(expected$value))
        expect_identical(colnames(actual$value), c("A", "A.1"))
        expect_identical(is.na(actual$value), is.na(expected$value))
        expect_identical(is.nan(actual$value), is.nan(expected$value))
        expect_identical(actual$value[position, ], expected$value[position, ])
        expect_equal(actual$value[-position, , drop = FALSE],
          expected$value[-position, , drop = FALSE], tolerance = 1e-14)
      }
    }
  }
  for (position in c(1L, 2L)) {
    values <- c(-1, 0, 1)
    values[position] <- Inf
    fit <- .bridge_repetition_test_fit(.bridge_repetition_test_raw(values))
    actual <- capture(post_prob(fit, fit, prior_prob = weights))
    expected <- capture(bridgesampling::post_prob(c(Inf, Inf),
      prior_prob = normalized, model_names = c("fit", "fit.1")))
    expect_s3_class(actual$value, "error")
    expect_identical(class(actual$value), class(expected$value))
    expect_identical(conditionMessage(actual$value), conditionMessage(expected$value))
    expect_identical(actual$warnings, expected$warnings)
  }
})
