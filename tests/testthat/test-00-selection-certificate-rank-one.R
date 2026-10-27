context("Rank-one product interval refinement")
skip_on_cran()
source(testthat::test_path("helper-selection-kernel.R"))

.test_rank_one_fallback <- function(block, ..., sign = 1L) {

  .test_cert_native_cluster(block, ..., sign = sign,
    rules = .selection_joint_cluster_quadrature_rules(c(1L, 3L)))
}

.test_rank_one_refusal <- function(result) {

  expect_true(is.nan(result$log_normalizer))
  expect_identical(result$log_density, -Inf)
  expect_identical(result$relative_mcse, Inf)
  expect_identical(result$relative_change, Inf)
  expect_false(.test_cert_deterministic(result))
}

test_that("all fourteen independently checked QMC outliers satisfy their requested error", {

  fixture <- readRDS(test_path("fixtures", "selection-rank-one-outliers.rds"))
  for (name in names(fixture$blocks)) {
    block <- fixture$blocks[[name]]
    at <- match(name, paste(fixture$outliers$set, fixture$outliers$block, sep = "_"))
    tolerance <- fixture$outliers$tolerance[[at]]
    reference <- fixture$references[[name]]$reference
    result <- .test_cert_native_cluster(block, tolerance = tolerance)
    expect_true(.test_cert_deterministic(result), info = name)
    expect_identical(result$relative_mcse, 0, info = name)
    expect_lte(result$relative_change, tolerance, label = name)
    expect_lte(abs(expm1(result$log_normalizer - reference)),
      result$relative_change + 1e-10, label = name)
    expect_identical(.test_cert_native_cluster(block, tolerance = tolerance, seed = 91L), result)
  }
})

test_that("signed bins, mixed scales and loading signs preserve independent integrals", {

  group <- function(mean, sd, sei, loading, count = 1L) {
    .test_cert_parallel_group(count, loading, sd, mean, sei)
  }
  # One Gaussian convolution, including a narrow finite band.
  for (loading in c(-1, 1)) {
    block <- .test_cert_block(list(group(.7, .3, 1.4, loading)), c(0, 1, 0), c(1.1, -.2, -Inf))
    marginal_sd <- sqrt(.3^2 + loading^2)
    reference <- log(pnorm((1.1 * 1.4 - .7) / marginal_sd) -
      pnorm((-.2 * 1.4 - .7) / marginal_sd))
    actual <- .test_rank_one_fallback(block)
    expect_true(.test_cert_deterministic(actual))
    expect_lte(abs(expm1(actual$log_normalizer - reference)), actual$relative_change + 1e-11)
    mirrored <- block
    mirrored$mean <- -block$mean
    mirrored$loading <- -block$loading
    expect_identical(.test_rank_one_fallback(mirrored, sign = -1L), actual)
  }
  block <- .test_cert_block(list(group(-.2, .6, .8, .7, 3L),
    group(.5, .4, 1.3, -.4, 2L), group(.1, .9, .7, 0)), c(.2, 1, .4), c(1, -.5, -Inf))
  reference <- .test_cert_reference_parallel(block) + .test_cert_log_probability(.1, .9, .7, block$omega, block$z_lower)
  result <- .test_rank_one_fallback(block)
  expect_lte(abs(expm1(result$log_normalizer - reference)), result$relative_change + 1e-10)
  # Zero columns compact to the same rank-one target.
  compacted <- block
  compacted$loading <- cbind(0, block$loading, 0)
  compacted$rank <- 3L
  factor <- .test_cert_native_factor(compacted,
    rules = .selection_joint_cluster_quadrature_rules(c(1L, 3L, 5L)))
  expect_lte(abs(expm1(factor$log_normalizer - reference)), .005)
})

test_that("arithmetic ambiguity and inadequate work refuse rather than accept", {

  block <- .test_cert_block(list(.test_cert_parallel_group(1L, 1, .01, 0)),
    c(0, 1, 0), c(2, 0, -Inf))
  centered <- .test_rank_one_fallback(block)
  expect_true(.test_cert_deterministic(centered))
  overflow <- block
  overflow$sei <- 4
  overflow$z_lower <- c(1e308, 0, -Inf)
  overflow$z_upper <- c(Inf, 1e308, 0)
  .test_rank_one_refusal(.test_rank_one_fallback(overflow))
  .test_rank_one_refusal(.test_rank_one_fallback(block, initial = 4L, points = 4L, scrambles = 2L))
  initial <- .test_cert_block(lapply(seq_len(8L), function(j) {
    .test_cert_parallel_group(1L, 1, .2, .1 * j)
  }), c(1, .2))
  .test_rank_one_refusal(.test_rank_one_fallback(initial, initial = 4L, points = 4L, scrambles = 2L))
  zero <- block; zero$omega[] <- 0
  .test_rank_one_refusal(.test_rank_one_fallback(zero))
})

test_that("constant and factor-free exact returns retain their identities", {

  block <- .test_cert_block(list(.test_cert_parallel_group(3L, 1, .4, .2)), c(.7, .7))
  actual <- .test_rank_one_fallback(block)
  expect_identical(actual$log_normalizer, 3 * log(.7))
  expect_identical(actual$relative_change, 0)
  block$loading[] <- 0
  block$omega <- c(1, .2)
  actual <- .test_rank_one_fallback(block)
  expect_lte(abs(actual$log_normalizer - 3 * log(.2 + .8 * pnorm(.2 / .4))), 1e-13)
  expect_identical(actual$relative_change, 0)
})

test_that("unavailable diagnostics are rejected by the shared factor and conditioned judges", {

  unavailable <- list(log_density = -Inf, log_normalizer = NaN,
    relative_mcse = Inf, relative_change = Inf)
  plan <- list(relative_tolerance = .005)
  expect_error(.selection_joint_factor_check_diagnostics(unavailable, plan, "cluster"),
    "normalizer was rejected by diagnostics", fixed = TRUE)
  expect_error(.selection_conditioned_sampling_diagnostics(unavailable, plan),
    "normalizer was rejected by diagnostics", fixed = TRUE)
  expect_error(.selection_joint_factor_check_diagnostics(unavailable, plan, "cluster"),
    "deterministic work budget", fixed = TRUE)
})


test_that("subnormal positive excess weights preserve the scale-equivalent target", {

  block <- .test_cert_block(list(.test_cert_parallel_group(2L, 1, .3, -.7)), c(1, .5))
  reference <- .test_cert_reference_parallel(block)
  tiny <- .Machine$double.xmin * .Machine$double.eps
  scaled <- block
  scaled$omega <- c(2 * tiny, tiny)
  ordinary <- .test_rank_one_fallback(block)
  result <- .test_rank_one_fallback(scaled)
  shift <- block$rows * log(2 * tiny)
  expect_true(.test_cert_deterministic(result))
  expect_identical(result$relative_mcse, 0)
  expect_lte(abs(expm1(result$log_normalizer - shift - reference)), result$relative_change + 1e-10)
  expect_lt(abs(result$log_normalizer - ordinary$log_normalizer - shift), 1e-8)
  expect_lt(result$log_normalizer, -1400)
})

test_that("full translated sign and loading grids refuse unresolved mean arithmetic", {

  reference <- log(pnorm(2 / sqrt(1 + .01^2)) - .5)
  for (origin in c(0, 1e16, 2^53 - 2, 2^53, -2^53)) for (loading in c(-1, 1)) {
    for (sign in c(-1L, 1L)) {
      block <- .test_cert_block(list(.test_cert_parallel_group(1L, sign * loading, .01, sign * origin)),
        c(0, 1, 0), c(origin + 2, origin, -Inf))
      if (origin == 0) {
        result <- .test_rank_one_fallback(block, sign = sign)
        expect_true(.test_cert_deterministic(result))
        expect_lte(abs(expm1(result$log_normalizer - reference)), result$relative_change + 1e-10)
      } else for (points in c(4L, 8192L)) {
        .test_rank_one_refusal(.test_rank_one_fallback(block, initial = 4L, points = points, sign = sign))
      }
    }
  }
})

test_that("rank-one fallback batches retain exact serial and parallel identity", {

  .with_kernel_thread_worker({
    S <- 32L
    K <- 2L
    means <- matrix(rep(seq(-.6, -.5, length.out = S), K), S)
    rules <- .selection_joint_cluster_quadrature_rules(c(1L, 3L, 5L))
    evaluate <- function(threads, factor) {
      RoBMA::RoBMA.options(native_threads = threads)
      rank <- if (factor) 3L else 1L
      design <- BayesTools::selection_qmc_design(2L * rank, 2048L, 8L, seed = 1L)
      loading <- if (factor) matrix(rep(c(rep(0, K), rep(1, K), rep(0, K)), each = S), S) else matrix(1, S, K)
      common <- list(rep(0, K), means, matrix(.02, S, K), loading)
      args <- c(common, list(rep(1, K), matrix(rep(c(1, 0), each = S), S), c(0, -Inf), c(Inf, 0),
        rep(1L, K), 1L, TRUE, 1L, rules$nodes, rules$log_weights, as.double(rules$orders),
        as.double(design), 4L, 2048L, 8L, .005, TRUE, 0L))
      do.call('.Call', c(list(if (factor) 'RoBMA_selnorm_factor_step_loglik_batch' else 'RoBMA_selnorm_cluster_step_loglik_batch'), args, list(PACKAGE = 'RoBMA')))
    }
    for (factor in c(FALSE, TRUE)) {
      serial <- evaluate(1L, factor)
      threaded <- evaluate(2L, factor)
      expect_identical(threaded, serial)
      expect_true(all(attr(serial, 'deterministic')))
      expect_true(all(serial$relative_mcse == 0))
    }
  })
})

test_that("finite-bin peaks next to the zero knot respect sign and selection SE", {

  for (se in c(.7, 1.3)) for (delta in c(-64, -1, 0, 1, 64) * .Machine$double.eps) {
    mean <- se + delta
    marginal_sd <- sqrt(.4^2 + .8^2)
    reference <- log(pnorm((2 * se - mean) / marginal_sd) - pnorm(-mean / marginal_sd))
    for (loading in c(-.8, .8)) for (sign in c(-1L, 1L)) {
      # The finite-bin peak moves just outside/inside the initial zero knot.
      block <- .test_cert_block(list(.test_cert_parallel_group(1L,
        sign * loading, .4, sign * mean, se)), c(0, 1, 0), c(2, 0, -Inf))
      actual <- .test_rank_one_fallback(block, sign = sign)
      expect_true(.test_cert_deterministic(actual))
      expect_lte(abs(expm1(actual$log_normalizer - reference)), actual$relative_change + 1e-11)
    }
  }
})

test_that("large and small scales retain the declared convolution and log normalizer", {

  reference <- log(pnorm(1.5 / sqrt(1 + .3^2)) - pnorm(-.5 / sqrt(1 + .3^2)))
  for (scale in c(1e-100, 1e100)) {
    block <- .test_cert_block(list(.test_cert_parallel_group(1L, scale, .3 * scale,
      .5 * scale, scale)), c(0, 1, 0), c(2, 0, -Inf))
    result <- .test_rank_one_fallback(block)
    expect_true(.test_cert_deterministic(result))
    expect_lte(abs(expm1(result$log_normalizer - reference)), result$relative_change + 1e-10)
  }
  block <- .test_cert_block(list(.test_cert_parallel_group(2L, 1, .3, -.7)), c(1, .5))
  reference <- .test_cert_reference_parallel(block)
  for (scale in c(1e-250, 1e250)) {
    scaled <- block
    scaled$omega <- scale * block$omega
    result <- .test_rank_one_fallback(scaled)
    expect_true(.test_cert_deterministic(result))
    expect_lte(abs(expm1(result$log_normalizer - block$rows * log(scale) - reference)),
      result$relative_change + 1e-10)
    expect_gt(abs(result$log_normalizer), 1100)
  }
})
