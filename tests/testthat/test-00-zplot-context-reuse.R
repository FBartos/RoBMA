.zplot_context_reuse_selection <- function(sei, weights = c(1, .2), cutoff = 0) {

  prior <- BayesTools::prior_weightfunction("one-sided", stats::pnorm(cutoff, lower.tail = FALSE),
    BayesTools::wf_fixed(weights),
    model = BayesTools::selection_model(weight_rule = "product"))
  selection <- .selection_spec(list(outcome = list(bias = prior)), rep(0, length(sei)), sei, "positive")
  selection$omega <- matrix(selection$fixed_omega[1L, ], 1L)
  selection$alpha <- 0
  selection$phack_kind <- 0L
  selection$kernel_mode <- SELKERNEL_STEP
  selection$vector_rule <- 0L
  selection$use_normal <- FALSE
  BayesTools::selection_context_validate(selection, n_samples = 1L)
}


.zplot_context_reuse_args <- function() {

  sei <- c(.7, 1.3)
  means <- rbind(c(-.4, .1), c(0, 0), c(.3, .5), c(.8, -.2),
    c(-.1, .9), c(.2, -.3), c(-.2, -.1), c(.6, .3))
  list(z = c(-.7, .3, 1.6, 2.4), mean = colMeans(means),
    covariance = matrix(c(.61, .7, .7, 1.81), 2L),
    context_factor = matrix(0, 1L, 2L), sei = sei,
    selection = .zplot_context_reuse_selection(sei),
    control = set_selection_likelihood_control(relative_tolerance = .001),
    context_means = means)
}


# Disable only the new fusion assignment for the original two-pass controller.
# Its denominator and all error-ledger inputs then come from the mass-only
# pass, including when collected and mass-only log arithmetic round differently.
.zplot_context_reuse_run <- function(args, separate = FALSE, alter = NULL, capture_error = FALSE) {

  calls <- list()
  projection <- .zplot_context_projection
  if (separate) {
    replaced <- 0L
    disable_reuse <- function(expression) {
      if (!is.call(expression)) return(expression)
      if (identical(expression[[1L]], as.name("<-")) &&
          identical(expression[[2L]], as.name("own_normalizer"))) {
        expression[[3L]] <- FALSE
        replaced <<- replaced + 1L
        return(expression)
      }
      as.call(lapply(as.list(expression), disable_reuse))
    }
    body(projection) <- disable_reuse(body(projection))
    stopifnot(replaced == 1L)
  }
  env <- new.env(parent = environment(projection))
  env$.Call <- function(.NAME, ..., PACKAGE) {

    native_args <- list(...)
    if (.NAME != "RoBMA_selnorm_context_star_zplot") {
      return(base::.Call(.NAME, ..., PACKAGE = PACKAGE))
    }
    collecting <- any(is.finite(native_args[[10L]]))
    own <- collecting && is.null(native_args[[9L]])
    value <- do.call(base::.Call, c(list(.NAME), native_args, list(PACKAGE = PACKAGE)))
    if (!is.null(alter)) value <- alter(value, native_args)
    calls[[length(calls) + 1L]] <<- list(order = length(native_args[[11L]]),
      own = own, collecting = collecting, denominator = native_args[[9L]],
      weights = native_args[[10L]], raw = value$compact_log_normalizers, means = native_args[[1L]],
      bounded_cdf = isTRUE(attr(native_args[[11L]], "bounded_cdf")))
    value
  }
  environment(projection) <- env
  value <- if (capture_error) tryCatch(do.call(projection, args), error = identity) else
    do.call(projection, args)
  list(value = value, calls = calls)
}


.zplot_context_reuse_expect_equal <- function(actual, expected) {

  expect_identical(is.null(actual$value), is.null(expected$value))
  if (is.null(actual$value)) return(invisible(NULL))
  expect_identical(actual$value$quadrature, expected$value$quadrature)
  expect_identical(actual$value$bounded_cdf, expected$value$bounded_cdf)
  expect_identical(names(actual$value$integration_error), names(expected$value$integration_error))
  # Only coefficient-subtraction roundoff differs between the two routes.
  expect_equal(actual$value, expected$value, tolerance = 1e-10)
}


test_that("same-rule context reuse preserves independent checks and the full error ledger", {

  args <- .zplot_context_reuse_args()
  actual <- .zplot_context_reuse_run(args)
  expected <- .zplot_context_reuse_run(args, separate = TRUE)
  .zplot_context_reuse_expect_equal(actual, expected)
  expect_type(actual$value, "list")
  own <- which(vapply(actual$calls, function(call) call$own, logical(1L)))
  expect_gt(length(own), 0L)
  for (index in own) {
    reference <- Filter(function(call) !call$collecting && call$order == actual$calls[[index]]$order &&
      identical(call$means, actual$calls[[index]]$means), expected$calls)
    expect_length(reference, 1L)
    expect_identical(actual$calls[[index]]$raw, reference[[1L]]$raw)
    anchor <- as.double(actual$calls[[index]]$raw)
    # A distinct numerator rule uses the exact same denominator. Its mass
    # and curve comparisons remain independent of the fine own-A mass of one.
    expect_true(any(vapply(actual$calls, function(call) {
      call$collecting && call$order < actual$calls[[index]]$order &&
        identical(call$denominator, anchor)
    }, logical(1L))))
  }
  expect_gt(actual$value$integration_error[["denominator"]], 0)
  expect_gt(actual$value$mass_error, 0)

  # Full-covariance bivariate orthants, independent of the native factor rule.
  covariance <- args$covariance
  sd <- sqrt(diag(covariance))
  schur <- diag(covariance) - covariance[1L, 2L]^2 / rev(diag(covariance))
  weight <- .2
  reference <- vapply(args$z, function(z) mean(apply(args$context_means, 1L, function(mean) {
    joint <- as.numeric(mvtnorm::pmvnorm(c(0, 0), c(Inf, Inf), mean, sigma = covariance,
      algorithm = mvtnorm::TVPACK()))
    denominator <- weight^2 + weight * (1 - weight) * sum(stats::pnorm(mean / sd)) +
      (1 - weight)^2 * joint
    x <- z * args$sei
    companion <- rev(mean) + covariance[1L, 2L] / diag(covariance) * (x - mean)
    selected <- weight + (1 - weight) * stats::pnorm(companion / sqrt(rev(schur)))
    observed <- if (z >= 0) 1 else weight
    mean(args$sei * stats::dnorm(x, mean, sd) * observed * selected) / denominator
  })), numeric(1L))
  expect_lte(max(abs(as.numeric(actual$value$density) - reference)),
    actual$value$integration_error[["absolute"]] + 1e-8)
})


test_that("unequal context rules keep the supplied denominator route", {

  args <- .zplot_context_reuse_args()
  args$covariance <- diag(args$sei^2) + .01 * tcrossprod(args$sei)
  actual <- .zplot_context_reuse_run(args)
  expected <- .zplot_context_reuse_run(args, separate = TRUE)
  .zplot_context_reuse_expect_equal(actual, expected)
  expect_type(actual$value, "list")
  expect_false(actual$value$quadrature[["numerator"]] == actual$value$quadrature[["denominator"]])
  expect_false(any(vapply(actual$calls, function(call) call$own, logical(1L))))
  expect_true(any(vapply(actual$calls, function(call) call$collecting &&
    !is.null(call$denominator), logical(1L))))
})


test_that("context reuse retains CDF fallback and diagnostic refusal", {

  args <- .zplot_context_reuse_args()
  uncertain_cdf <- function(value, native_args) {
    if (isTRUE(attr(native_args[[11L]], "bounded_cdf"))) {
      raw <- value$compact_log_normalizers
      attr(raw, "cdf_factor_log_error") <- rep(.005, length(raw))
      attr(raw, "cdf_log_error") <- rep(.01, length(raw))
      attr(raw, "cdf_numerator_log_error") <- rep(.005, length(raw))
      value$compact_log_normalizers <- raw
    }
    value
  }
  actual <- .zplot_context_reuse_run(args, alter = uncertain_cdf)
  expected <- .zplot_context_reuse_run(args, separate = TRUE, alter = uncertain_cdf)
  .zplot_context_reuse_expect_equal(actual, expected)
  expect_type(actual$value, "list")
  expect_false(actual$value$bounded_cdf)
  expect_identical(actual$value$integration_error[["cdf_curve"]], 0)
  for (separate in c(FALSE, TRUE)) {
    expect_error(.zplot_context_reuse_run(args, separate, alter = function(value, native_args) {
      attr(value$compact_log_normalizers, "cdf_log_error") <- NULL
      value
    }), "CDF integration diagnostics are unavailable", fixed = TRUE)
    unavailable <- .zplot_context_reuse_run(args, separate, alter = function(value, native_args) {
      value$available <- FALSE
      value
    })
    expect_null(unavailable$value)
  }
})


test_that("context reuse preserves unsupported and zero-weight refusals", {

  args <- .zplot_context_reuse_args()
  for (case in c("zero", "covariance", "retained")) {
    current <- args
    if (case == "zero") current$selection <- .zplot_context_reuse_selection(args$sei, c(1, 0))
    if (case == "covariance") current$covariance[1L, 2L] <- current$covariance[2L, 1L] <- -.364
    if (case == "retained") {
      current$context_means <- NULL
      current$context_factor <- matrix(.2, 2L, 2L)
    }
    actual <- .zplot_context_reuse_run(current)
    expected <- .zplot_context_reuse_run(current, separate = TRUE)
    .zplot_context_reuse_expect_equal(actual, expected)
    expect_null(actual$value)
  }
})


test_that("continuous context pruning keeps the same checked result", {

  args <- .zplot_context_reuse_args()
  args$context_means <- NULL
  args$context_factor <- matrix(.8, 1L, 2L)
  actual <- .zplot_context_reuse_run(args)
  expected <- .zplot_context_reuse_run(args, separate = TRUE)
  .zplot_context_reuse_expect_equal(actual, expected)
  expect_type(actual$value, "list")
  expect_gt(actual$value$integration_error[["pruning"]], 0)
  expect_true(any(vapply(actual$calls, function(call) call$collecting &&
    any(call$weights == -Inf), logical(1L))))
})


test_that("log-product inputs retain the existing two-pass refusal", {

  args <- .zplot_context_reuse_args()
  # Disabling the direct telescoping product exercises the separate log-scale
  # mass-only and collected integrals, whose summation orders differ.
  args$selection$telescope_probabilities <- FALSE
  # The old collected A can replace a mass-only cache entry by a few ulps.
  # On affected platforms this fails the existing shared-D identity check.
  # Fusion must not bypass that refusal or alter a platform's accepted result.
  actual <- .zplot_context_reuse_run(args, capture_error = TRUE)
  expected <- .zplot_context_reuse_run(args, separate = TRUE, capture_error = TRUE)
  if (inherits(expected$value, "error")) {
    expect_s3_class(actual$value, "error")
    expect_identical(conditionMessage(actual$value), conditionMessage(expected$value))
  } else {
    .zplot_context_reuse_expect_equal(actual, expected)
  }
  expect_false(any(vapply(actual$calls, function(call) call$own, logical(1L))))
})


test_that("extreme absolute weight scales stay on the original route", {

  args <- .zplot_context_reuse_args()
  for (scale in c(1e-200, 1e200)) {
    # The native context seam accepts unnormalized weights as well.
    args$selection$omega <- matrix(c(1, .2) * scale, 1L)
    actual <- .zplot_context_reuse_run(args)
    expected <- .zplot_context_reuse_run(args, separate = TRUE)
    .zplot_context_reuse_expect_equal(actual, expected)
    expect_type(actual$value, "list")
    expect_false(any(vapply(actual$calls, function(call) call$own, logical(1L))))
  }
})


test_that("increasing and nonmonotone weights retain the original controller", {

  args <- .zplot_context_reuse_args()
  for (weights in list(c(.2, 1), c(1, .2, .7))) {
    if (length(weights) == 3L) {
      args$selection <- .zplot_context_reuse_selection(args$sei, weights, cutoff = c(1, -1))
    } else {
      args$selection$omega <- matrix(weights, 1L)
    }
    actual <- .zplot_context_reuse_run(args)
    expected <- .zplot_context_reuse_run(args, separate = TRUE)
    .zplot_context_reuse_expect_equal(actual, expected)
    expect_false(any(vapply(actual$calls, function(call) call$own, logical(1L))))
  }
})


test_that("own-A context collection is invariant to thread count", {

  args <- .zplot_context_reuse_args()
  previous <- RoBMA.get_option("native_threads")
  on.exit(RoBMA.options(native_threads = previous), add = TRUE)
  RoBMA.options(native_threads = 1L)
  reference <- .zplot_context_reuse_run(args)
  for (threads in c(2L, 3L)) {
    RoBMA.options(native_threads = threads)
    actual <- .zplot_context_reuse_run(args)
    expect_identical(actual, reference)
  }
})


test_that("extreme tail contexts keep relative refusal and checked absolute acceptance", {

  mean <- c(2.35938851263664, 2.0511642238707, 1.62233919836236)
  sei <- c(.107238052947636, .0871779788708135, .0806225774829855)
  covariance <- matrix(c(.0115, .00654415769981133, .0060520657630267,
    .00654415769981133, .0076, .00491995934942556,
    .0060520657630267, .00491995934942556, .0065), 3L)
  args <- list(z = seq(-6, 10, length.out = 81L), mean = mean, covariance = covariance,
    context_factor = matrix(0, 1L, 3L), sei = sei,
    selection = .zplot_context_reuse_selection(sei, c(1, .914973906385372), cutoff = 1.96),
    control = set_selection_likelihood_control(), context_means = matrix(rep(mean, each = 8L), 8L))
  relative <- .zplot_context_reuse_run(args)
  .zplot_context_reuse_expect_equal(relative, .zplot_context_reuse_run(args, separate = TRUE))
  expect_null(relative$value)
  args$absolute_tolerance <- 1.99110200352278e-5
  actual <- .zplot_context_reuse_run(args)
  expected <- .zplot_context_reuse_run(args, separate = TRUE)
  .zplot_context_reuse_expect_equal(actual, expected)
  expect_type(actual$value, "list")
  expect_lt(max(actual$value$density), 1e-20)
  expect_lte(actual$value$integration_error[["absolute"]], args$absolute_tolerance)
  expect_lte(actual$value$mass_error, args$control$relative_tolerance)
})
