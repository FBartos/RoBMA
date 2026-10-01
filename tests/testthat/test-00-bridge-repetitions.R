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

test_that("recorded nonfinite drop policy differs from legacy NA propagation", {

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
  expect_equal(bf(dropped_fit, comparison_fit, log = TRUE)$bf, c(-2, -2))
  expect_identical(logml(dropped_fit), dropped$logml)
  expect_equal(dropped_fit$marglik$repetitions$logml, c(-2, NA, -1))
  dropped_fit$marglik$aggregation$nonfinite_policy <- "error"
  expect_error(bf(dropped_fit, comparison_fit), "must contain only finite values", fixed = TRUE)

  expect_warning(one_left <- .bridge_repetition_test_stored(
    .bridge_repetition_test_raw(c(NA_real_, -1)), "drop"
  ), "Dropped non-finite", fixed = TRUE)
  one_left <- .bridge_repetition_test_fit(one_left)
  scalar <- .bridge_repetition_test_fit(.bridge_repetition_test_raw(0))
  expect_s3_class(bf(one_left, scalar), "bf_bridge_list")
  expect_equal(dim(post_prob(one_left, scalar)), c(1L, 2L))
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
