test_that("repeated context-rule probabilities retain their finite-rule density", {

  # Repeated nodes deliberately make every row probability identical within
  # its rule. At node zero the common factor is fixed, leaving independent
  # unit-variance rows whose selected marginal densities are analytic.
  means <- rbind(c(-.4, .7), c(.3, -.2), c(-3, 2), c(4, -4),
    c(-.1, .8), c(.2, -.5), c(-40, 40), c(40, -40))
  contexts <- nrow(means)
  z <- seq(-4, 4, length.out = 17L)
  order <- 31L
  for (sign in c(-1L, 1L)) {
    for (telescope in c(FALSE, TRUE)) {
      omega <- c(1, .2)
      nodes <- rep(0, order)
      attr(nodes, "bounded_cdf") <- FALSE
      actual <- .Call("RoBMA_selnorm_context_star_zplot", means,
        c(1.25, .25, 1.25), c(1, 1), omega, c(0, -Inf), c(Inf, 0),
        sign, telescope, NULL, rep(-log(contexts), contexts), nodes,
        rep(-log(order), order), z, 0, PACKAGE = "RoBMA")
      normalizer <- omega[[2L]] + (omega[[1L]] - omega[[2L]]) *
        stats::pnorm(sign * means)
      expected <- vapply(z, function(point) {
        weight <- if (sign * point >= 0) omega[[1L]] else omega[[2L]]
        mean(exp(stats::dnorm(point, means, log = TRUE) + log(weight) - log(normalizer)))
      }, numeric(1L))
      expect_true(actual$available)
      expect_equal(actual$density, expected, tolerance = 1e-12)
      expect_equal(actual$mass, 1, tolerance = 1e-14)
      expect_equal(as.double(actual$compact_log_normalizers), rowSums(log(normalizer)),
        tolerance = 1e-12)
      expect_identical(actual$mixture_components, as.double(contexts * ncol(means) * order))
      expect_identical(attr(actual$compression_error, "omitted_mass_error"), 0)
    }
  }
})


.zplot_context_workspace_inputs <- function() {

  groups <- c(1L, 1L, 2L, 3L)
  covariance <- diag(1, 4L) + .25 + .5 * outer(groups, groups, "==")
  rule <- .gauss_hermite_nodes(15L)
  list(means = matrix(seq(-.8, .9, length.out = 32L), 8L),
    covariance = as.double(covariance[lower.tri(covariance, diag = TRUE)]),
    sei = rep(1, 4L), omega = c(1, .2), lower = c(0, -Inf), upper = c(Inf, 0),
    sign = 1L, telescope = TRUE, denominator = NULL, weights = rep(-log(8), 8L),
    nodes = as.double(rule$nodes), log_weights = as.double(rule$log_weights),
    z = seq(-3, 3, length.out = 17L), allowance = 1e-6)
}


.zplot_context_workspace_call <- function(args, cache = NULL, budget = 0) {

  attr(args$nodes, "context_rule_cache") <- cache
  attr(args$nodes, "context_rule_cache_budget") <- as.double(budget)
  do.call(.Call, c(list("RoBMA_selnorm_context_star_zplot"), unname(args), list(PACKAGE = "RoBMA")))
}


.zplot_context_workspace_value <- function(value) {

  attr(value, "context_rule_cache") <- NULL
  attr(value, "context_rule_cache_hit") <- NULL
  value
}


test_that("raw context rules replay new denominators and retain the complete native ledger", {

  args <- .zplot_context_workspace_inputs()
  for (weights in list(c(1, .2), c(1, .2, .6))) {
    args$omega <- weights
    if (length(weights) == 3L) {
      args$lower <- c(1, -1, -Inf)
      args$upper <- c(Inf, 1, -1)
      args$telescope <- FALSE
    }
    for (bounded in c(FALSE, TRUE)) {
      attr(args$nodes, "bounded_cdf") <- bounded
      initial <- .zplot_context_workspace_call(args, budget = 1024^2)
      cache <- attr(initial, "context_rule_cache")
      expect_true(initial$available)
      expect_type(cache, "list")
      expect_null(attributes(cache[[1L]]))
      expect_null(attributes(cache[[2L]]))
      expect_lte(as.double(utils::object.size(cache)), 1024^2)
      expect_false(attr(initial, "context_rule_cache_hit"))
      changed <- args
      changed$denominator <- as.double(initial$compact_log_normalizers) + seq(0, .01, length.out = 8L)
      changed$z <- seq(-3.5, 3.5, length.out = 25L)
      changed$allowance <- 0
      cached <- .zplot_context_workspace_call(changed, cache)
      uncached <- .zplot_context_workspace_call(changed)
      expect_true(attr(cached, "context_rule_cache_hit"))
      expect_identical(.zplot_context_workspace_value(cached), uncached)
    }
  }
})


test_that("raw context keys reject stale inputs and preserve uncached refusals", {

  args <- .zplot_context_workspace_inputs()
  initial <- .zplot_context_workspace_call(args, budget = 1024^2)
  cache <- attr(initial, "context_rule_cache")
  variants <- list()
  changed <- args; changed$means[1L, 1L] <- changed$means[1L, 1L] + .1
  variants$means <- changed
  changed <- args; changed$covariance[[1L]] <- changed$covariance[[1L]] + .05
  variants$covariance <- changed
  changed <- args; changed$sei[[1L]] <- 1.01
  variants$scales <- changed
  changed <- args; changed$omega[[2L]] <- .3
  variants$weights <- changed
  changed <- args; changed$lower[[1L]] <- changed$upper[[2L]] <- .1
  variants$cutoffs <- changed
  changed <- args; changed$sign <- -1L
  variants$direction <- changed
  changed <- args; changed$telescope <- FALSE
  variants$telescope <- changed
  changed <- args; changed$weights[[1L]] <- -Inf
  variants$mask <- changed
  changed <- args; changed$weights <- rep(-Inf, 8L)
  variants$mass_only <- changed
  changed <- args; changed$nodes <- rev(changed$nodes)
  variants$nodes <- changed
  changed <- args; changed$log_weights <- changed$log_weights + log(.9)
  variants$rule_mass <- changed
  changed <- args; attr(changed$nodes, "bounded_cdf") <- TRUE
  variants$cdf <- changed
  changed <- args; changed$means <- changed$means[1:4, , drop = FALSE]; changed$weights <- changed$weights[1:4]
  variants$contexts <- changed
  changed <- args; rule <- .gauss_hermite_nodes(7L)
  changed$nodes <- as.double(rule$nodes); changed$log_weights <- as.double(rule$log_weights)
  variants$order <- changed
  changed <- args; changed$means <- changed$means[, 1:3, drop = FALSE]; changed$sei <- rep(1, 3L)
  covariance <- matrix(c(1.75, .75, .25, .75, 1.75, .25, .25, .25, 1.75), 3L)
  changed$covariance <- as.double(covariance[lower.tri(covariance, diag = TRUE)])
  variants$dimension <- changed
  for (name in names(variants)) {
    actual <- .zplot_context_workspace_call(variants[[name]], cache)
    expected <- .zplot_context_workspace_call(variants[[name]])
    expect_null(attr(actual, "context_rule_cache"), info = name)
    expect_identical(actual, expected, info = name)
  }
  # The saved key owns numeric values, not the input objects or their attributes.
  key <- cache[[1L]]
  args$means[1L, 1L] <- 10
  attr(args$nodes, "bounded_cdf") <- TRUE
  expect_identical(cache[[1L]], key)
})


test_that("raw context cache bounds and malformed payloads are explicit", {

  args <- .zplot_context_workspace_inputs()
  uncached <- .zplot_context_workspace_call(args)
  expect_identical(.zplot_context_workspace_call(args, budget = 1), uncached)
  initial <- .zplot_context_workspace_call(args, budget = 1024^2)
  cache <- attr(initial, "context_rule_cache")
  expect_error(.zplot_context_workspace_call(args, list(1L, 0)), "numeric key and payload", fixed = TRUE)
  bad <- cache; bad[[2L]] <- head(bad[[2L]], -1L)
  expect_error(.zplot_context_workspace_call(args, bad), "payload is invalid", fixed = TRUE)
  bad <- cache; bad[[2L]][[2L]] <- NaN
  expect_error(.zplot_context_workspace_call(args, bad), "payload is invalid", fixed = TRUE)
  bad <- cache; bad[[2L]][[1L]] <- -Inf
  expect_error(.zplot_context_workspace_call(args, bad), "raw normalizers are invalid", fixed = TRUE)
  expect_error(.zplot_context_workspace_call(args, budget = -1), "finite and nonnegative", fixed = TRUE)
  args$weights[c(1L, 4L)] <- -Inf
  initial <- .zplot_context_workspace_call(args, budget = 1024^2)
  actual <- .zplot_context_workspace_call(args, attr(initial, "context_rule_cache"))
  expect_true(attr(actual, "context_rule_cache_hit"))
  expect_identical(.zplot_context_workspace_value(actual), .zplot_context_workspace_call(args))
})


test_that("context fragments preserve child and absorbed-singleton component order", {

  previous_threads <- RoBMA.get_option("native_threads")
  withr::defer(RoBMA.options(native_threads = previous_threads))
  groups <- c(1L, 1L, 2L, 3L)
  covariance <- diag(1, 4L) + .25 + .5 * outer(groups, groups, "==")
  means <- matrix(seq(-.8, .9, length.out = 32L), 8L)
  # Child axes of the first two rows are fixed at zero by this finite rule;
  # the two singleton child variances are analytically absorbed by the kernel.
  sd <- c(1, 1, sqrt(1.5), sqrt(1.5))
  z <- seq(-3, 3, length.out = 17L)
  order <- 15L
  normalizer <- .2 + .8 * stats::pnorm(sweep(means, 2L, sd, "/"))
  expected <- vapply(z, function(point) {
    weight <- if (point >= 0) 1 else .2
    mean(exp(stats::dnorm(point, means, matrix(sd, 8L, 4L, byrow = TRUE), log = TRUE) +
      log(weight) - log(normalizer)))
  }, numeric(1L))
  baseline <- NULL
  for (threads in c(1L, 2L)) {
    RoBMA.options(native_threads = threads)
    actual <- .Call("RoBMA_selnorm_context_star_zplot", means,
      as.double(covariance[lower.tri(covariance, diag = TRUE)]), rep(1, 4L),
      c(1, .2), c(0, -Inf), c(Inf, 0), 1L, TRUE, NULL, rep(-log(8), 8L),
      rep(0, order), rep(-log(order), order), z, 0, PACKAGE = "RoBMA")
    expect_true(actual$available)
    expect_equal(actual$density, expected, tolerance = 1e-12)
    expect_equal(actual$mass, 1, tolerance = 1e-14)
    expect_equal(as.double(actual$compact_log_normalizers), rowSums(log(normalizer)),
      tolerance = 1e-12)
    expect_identical(actual$mixture_components, as.double(8L * (2L * order^2 + 2L * order)))
    if (is.null(baseline)) baseline <- actual else expect_identical(actual, baseline)
  }
})
