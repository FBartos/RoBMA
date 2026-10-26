context("Rank-one product interval refinement")
skip_on_cran()

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
  for (origin in c(1e16, 2^53)) {
    translated <- block
    translated$mean <- origin
    translated$z_lower <- c(origin + 2, origin, -Inf)
    translated$z_upper <- c(Inf, origin + 2, origin)
    for (points in c(4L, 8192L)) {
      .test_rank_one_refusal(.test_rank_one_fallback(translated, initial = 4L, points = points))
    }
  }
  .test_rank_one_refusal(.test_rank_one_fallback(block, initial = 4L, points = 4L, scrambles = 2L))
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


