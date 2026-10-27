context("Selection certificate QMC proposals")
skip_on_cran()

# The Gaussian part of the product rule's randomized fallback proposes from the
# prior, a Gaussian with the tilted moments of the finest deterministic rung (mean and covariance
# of the latent coordinates under the selected integrand) and a copy with 1.5
# times their standard deviations. The moments are used for tensors, forests of
# depth at most one when at least three rungs were finite
# and the covariance is positive definite; every other block keeps the
# optimized-mode proposals. A finest rung that misses the exterior mass gives no
# moments either (test-00-selection-factor-tail-coverage.R keeps that case). The
# fallback also covers selection transitions; test-00-selection-qmc-tail-coverage.R
# checks the distant mass that the Gaussian proposals can all miss. The stopping
# rule and its diagnostics are unchanged.

# The native kernel of a deferred-geometry case: factor kernel, or the cluster
# kernel with its two-rung ladder.
.test_cert_proposals_run <- function(case, seed, initial, points) {

  if (case[["kind"]] == "factor") {
    .test_cert_native_factor(case[["block"]], seed = seed, initial = initial, points = points)
  } else {
    rules <- .selection_joint_cluster_quadrature_rules(c(15L, 31L))
    .test_cert_native_cluster(case[["block"]], seed = seed, initial = initial,
                              points = points, rules = rules)
  }
}

test_that("the tilted moments of a rung equal the plain tensor moments in loading-column order", {

  blocks <- .test_cert_small_blocks()
  for (name in names(blocks)) for (telescope in c(TRUE, FALSE)) {
    block  <- blocks[[name]]
    info   <- paste(name, "telescope", telescope)
    probe  <- .test_cert_probe(block, 9L, telescope = telescope, refine = FALSE)
    oracle <- .test_cert_tensor_moments(block, 9L)
    expect_true(probe[["has_moments"]], info = info)
    expect_lt(max(abs(probe[["moment_mean"]] - oracle[["mean"]])), 1e-12, label = info)
    expect_lt(max(abs(probe[["moment_covariance"]] - as.vector(oracle[["covariance"]]))), 1e-12,
              label = info)
    # not vacuous: the tilted density is not the prior
    expect_gt(max(abs(oracle[["covariance"]] - diag(block[["rank"]]))), .05, label = info)
  }
  # the structures reach the moments through the nested rule, the direct
  # tensor sum, or (without the telescoping products) the log-scale tensor sum
  expect_identical(.test_cert_probe(blocks[["two_leaves"]], 9L, refine = FALSE)[["path"]], .test_cert_paths[["nested"]])
  expect_identical(.test_cert_probe(blocks[["crossing"]], 9L, refine = FALSE)[["path"]], .test_cert_paths[["tensor_direct"]])
  expect_identical(.test_cert_probe(blocks[["two_leaves"]], 9L, telescope = FALSE, refine = FALSE)[["path"]],
                   .test_cert_paths[["tensor_log"]])
})

test_that("the rank-one cluster rule reports the mean and variance of its latent coordinate", {

  block <- .test_cert_small_blocks()[["rank_one"]]
  for (order in c(15L, 63L)) {
    probe  <- .test_cert_probe(block, order, cluster = TRUE, refine = FALSE)
    oracle <- .test_cert_tensor_moments(block, order)
    expect_true(probe[["has_moments"]], info = as.character(order))
    expect_lt(abs(probe[["moment_mean"]] - oracle[["mean"]]), 1e-12, label = order)
    expect_lt(abs(probe[["moment_covariance"]] - oracle[["covariance"]][1L, 1L]), 1e-12, label = order)
  }
})

test_that("permuting the loading columns of a forest permutes its moments", {

  block <- .test_cert_small_blocks()[["two_leaves"]]
  base  <- .test_cert_probe(block, 11L, refine = FALSE)
  for (permutation in list(c(3L, 1L, 2L), c(2L, 3L, 1L), c(3L, 2L, 1L))) {
    permuted <- block
    permuted[["loading"]] <- block[["loading"]][, permutation]
    moved  <- .test_cert_probe(permuted, 11L, refine = FALSE)
    oracle <- .test_cert_tensor_moments(permuted, 11L)
    info   <- paste(permutation, collapse = "")
    # loading column j of the permuted block is column permutation[j] of the original
    expect_lt(max(abs(moved[["moment_mean"]] - base[["moment_mean"]][permutation])), 1e-12, label = info)
    expect_lt(max(abs(matrix(moved[["moment_covariance"]], 3L) -
      matrix(base[["moment_covariance"]], 3L)[permutation, permutation])), 1e-12, label = info)
    expect_lt(max(abs(moved[["moment_mean"]] - oracle[["mean"]])), 1e-12, label = info)
    expect_lt(max(abs(moved[["moment_covariance"]] - as.vector(oracle[["covariance"]]))), 1e-12,
              label = info)
  }
  # the leaves are conditionally independent given the root, the root couples them
  covariance <- matrix(base[["moment_covariance"]], 3L)
  expect_gt(abs(covariance[2L, 3L]), 1e-4)
})

test_that("forests deeper than one level keep the optimized-mode proposals", {

  chain <- .test_cert_block(list(
    .test_cert_small_group(2L, c(.6, .4, .3), .33, .2),
    .test_cert_small_group(2L, c(.5, .4, 0), .36, -.1),
    .test_cert_small_group(2L, c(.7, 0, 0), .3, .1)), c(1, .2))
  probe <- .test_cert_probe(chain, 9L, refine = FALSE)
  expect_identical(probe[["path"]], .test_cert_paths[["nested"]])
  expect_false(probe[["has_moments"]])
  # a nested rule that left the direct range is a tensor rule, which has them
  expect_true(.test_cert_probe(chain, 9L, telescope = FALSE, refine = FALSE)[["has_moments"]])
})

# The deferred geometries (.test_cert_proposal_cases) live in helper-selection-certificate.R,
# shared with test-00-selection-qmc-stopping.R.

test_that("the rank-two reference agrees with the two-group reference", {

  block <- .test_cert_block(list(
    list(count = 40L, loading = c(-0.91463423771209651, 0.78986739681623297), sd = 0.28436590998117955,
         mean = -2.3718989917066349, sei = 1),
    list(count = 9L, loading = c(0.015058325926115159, -0.30020325671357212), sd = 0.019704307952072554,
         mean = 1.8604332933661574, sei = 1)),
    c(1, 0.05), c(1.96, -Inf))
  expect_lt(abs(.test_cert_reference_rank_two(block) - .test_cert_reference_two_groups(block)), 1e-5)
})

test_that("the fallbacks return finite estimates on deferred geometries", {

  cases <- .test_cert_proposal_cases()
  for (name in names(cases)) {
    case      <- cases[[name]]
    reference <- case[["reference"]](case[["block"]])
    if (!is.null(case[["pinned"]])) expect_lt(abs(reference - case[["pinned"]]), 1e-8, label = name)
    for (seed in 1:2) {
      actual <- if (case[["kind"]] == "factor") {
        .test_cert_native_factor(case[["block"]], seed = seed, initial = 256L, points = 4096L)
      } else {
        .test_cert_native_cluster(case[["block"]], seed = seed, initial = 256L, points = 4096L)
      }
      info <- paste(name, "seed", seed)
      if (case$block$rank == 1L) {
        expect_true(.test_cert_deterministic(actual), info = info)
        expect_identical(actual$relative_mcse, 0, info = info)
        expect_lte(abs(expm1(actual$log_normalizer - reference)), actual$relative_change + 1e-10,
          label = info)
      } else {
        expect_false(.test_cert_deterministic(actual), info = info)
        expect_gt(actual[["relative_mcse"]], 0, label = info)
      }
      error <- abs(expm1(actual[["log_normalizer"]] - reference))
      expect_true(all(is.finite(c(error, actual[["relative_mcse"]], actual[["relative_change"]]))),
                  info = info)
    }
  }
})

test_that("randomized uncertainty and mean accuracy calibrate across independent seeds", {

  # Empirical MCSE is not a bound on each design. These seeds and criteria were
  # fixed before their results were inspected; deeper multi-case verification
  # uses separate seeds. The former isolated forest_84 outlier motivates this
  # representative ordinary-test panel without selecting individual outcomes.
  case <- .test_cert_proposal_cases()[["forest_84"]]
  reference <- case[["reference"]](case[["block"]])
  seeds <- 101:140
  actual <- lapply(seeds, function(seed) {
    .test_cert_native_factor(case[["block"]], seed = seed, initial = 256L, points = 4096L)
  })
  log_ratio <- vapply(actual, function(result) result[["log_normalizer"]] - reference, numeric(1L))
  mcse <- vapply(actual, function(result) result[["relative_mcse"]], numeric(1L))
  error <- expm1(log_ratio)
  # Native MCSE is SE(Z_hat) / Z_hat, so this is the studentized discrepancy
  # (Z_hat - Z_reference) / SE(Z_hat), with the correct denominator.
  studentized <- -expm1(-log_ratio) / mcse
  expect_true(all(is.finite(c(error, mcse, studentized))))
  expect_true(all(mcse > 0))
  expect_false(any(vapply(actual, .test_cert_deterministic, logical(1L))))

  # Eight scrambles give the nominal t(7) benchmark. A one-sided exact 99%
  # binomial interval screens ensemble undercoverage, rather than requiring
  # every run to fall inside three estimated standard errors.
  exceedances <- sum(abs(studentized) > 3)
  nominal <- 2 * stats::pt(-3, df = 7)
  coverage <- stats::binom.test(exceedances, length(seeds), p = nominal,
    alternative = "greater", conf.level = .99)
  expect_lte(unname(coverage[["conf.int"]][1L]), nominal,
    label = paste("99% lower exceedance limit from", exceedances, "of", length(seeds), "designs"))

  # Independently of the reported MCSEs, a 99% across-seed t interval for mean
  # true-relative error must fit inside the existing .005 accuracy tolerance.
  # Both Student references are empirical approximations, not QMC guarantees.
  mean_error <- mean(error)
  half_width <- stats::qt(.995, df = length(seeds) - 1L) * stats::sd(error) / sqrt(length(seeds))
  mean_interval <- mean_error + c(-1, 1) * half_width
  expect_lt(max(abs(mean_interval)), .005,
    label = "largest absolute endpoint of the 99% mean-error interval")
  message(sprintf(paste0(
    "QMC calibration: %d/%d beyond 3 SE; 99%% exceedance lower limit %.6g; ",
    "mean relative error %.6g, 99%% CI [%.6g, %.6g]; RMSE %.6g; median MCSE %.6g."),
    exceedances, length(seeds), coverage[["conf.int"]][1L], mean_error,
    mean_interval[1L], mean_interval[2L], sqrt(mean(error^2)), stats::median(mcse)))
})

test_that("the fallbacks retain their uncertainty and accuracy on deferred geometries", {

  # Same policy as the forest_84 panel above. The reported MCSE is an empirical
  # standard error, not a bound, so accuracy is judged across independent seeds
  # fixed in advance. Every reference is independent of the production
  # proposals: the helpers' composite Gauss-Legendre integrals, and for the
  # stored geometries a transfer-matrix chain (chain3, chain4), a cyclic
  # transfer matrix (cycle7, cycle8), the change of variables to the two
  # row locations (band_a, band_b) and the one-dimensional reduction (cluster2).
  # The 99% exceedance screen is that of the forest_84 panel; the budgets are
  # those of the tests above.
  tolerance <- .005
  seeds <- 401:412
  nominal <- 2 * stats::pt(-3, df = 7)
  cases <- .test_cert_proposal_cases()
  cases[["forest_84"]] <- NULL
  stored <- readRDS(testthat::test_path("fixtures", "selection-certificate-laplace-reference.rds"))[["cases"]]
  references <- list(
    cycle7 = .test_cert_reference_cycle, cycle8 = .test_cert_reference_cycle,
    chain3 = .test_cert_reference_chain, chain4 = .test_cert_reference_chain,
    band_a = .test_cert_reference_band_pair, band_b = .test_cert_reference_band_pair,
    cluster2 = .test_cert_reference_parallel)
  stored <- Map(function(case, name) c(case, list(reference = references[[name]])),
                stored, names(stored))
  budget <- c(rep(list(4096L), length(cases)), rep(list(2048L), length(stored)))
  cases <- c(cases, stored)
  # At 2048 points per scramble the rank-seven and -eight cycles, the depth-four
  # chain and the narrow band cannot reach the tolerance (their MCSE is larger
  # than 0.4%); only their calibration is screened.
  imprecise <- c("cycle7", "cycle8", "chain4", "band_b")

  for (index in seq_along(cases)) {
    name <- names(cases)[index]
    case <- cases[[index]]
    reference <- case[["reference"]](case[["block"]])
    actual <- lapply(seeds, function(seed) {
      .test_cert_proposals_run(case, seed, 256L, budget[[index]])
    })
    if (case$block$rank == 1L) {
      expect_true(all(vapply(actual, .test_cert_deterministic, logical(1L))), info = name)
      expect_true(all(vapply(actual, identical, logical(1L), actual[[1L]])), info = name)
      expect_identical(actual[[1L]]$relative_mcse, 0, info = name)
      expect_lte(abs(expm1(actual[[1L]]$log_normalizer - reference)),
        actual[[1L]]$relative_change + 1e-10, label = name)
      expect_lte(actual[[1L]]$relative_change, tolerance, label = name)
      next
    }
    expect_false(any(vapply(actual, .test_cert_deterministic, logical(1L))), info = name)
    log_ratio <- vapply(actual, function(result) result[["log_normalizer"]] - reference, numeric(1L))
    mcse <- vapply(actual, function(result) result[["relative_mcse"]], numeric(1L))
    change <- vapply(actual, function(result) result[["relative_change"]], numeric(1L))
    error <- expm1(log_ratio)
    studentized <- -expm1(-log_ratio) / mcse
    expect_true(all(is.finite(c(error, mcse, change, studentized))), info = name)
    expect_true(all(mcse > 0), info = name)

    # False confidence is a run far from the reference under a small MCSE.
    exceedances <- sum(abs(studentized) > 3)
    coverage <- stats::binom.test(exceedances, length(seeds), p = nominal,
      alternative = "greater", conf.level = .99)
    expect_lte(unname(coverage[["conf.int"]][1L]), nominal,
      label = paste(name, ": 99% lower exceedance limit from", exceedances, "of",
        length(seeds), "designs"))
    if (name %in% imprecise) next

    # Runs that the scaled stopping rule accepts are within three tolerances, and
    # the mean error over the seeds lies inside the tolerance with 99% confidence.
    accepted <- pmax(.selection_qmc_z() * mcse, change) <= tolerance
    expect_lt(max(0, abs(error[accepted])), 3 * tolerance,
      label = paste(name, ": largest accepted error"))
    half_width <- stats::qt(.995, df = length(seeds) - 1L) * stats::sd(error) / sqrt(length(seeds))
    expect_lt(abs(mean(error)) + half_width, tolerance,
      label = paste(name, ": mean error and its 99% half-width"))
  }
})

test_that("blocks without usable moments retain finite diagnostics and nested stopping", {

  # Geometries from robma-claude-final7 (the build before the certificate) that
  # reach the fallback there and, here, with fewer than three finite rungs
  # (rank seven and eight tensors: two and one affordable rungs; a rank-one
  # cluster with a two-rung ladder), with a forest of depth two or three, or
  # with a moment covariance that is not positive definite. These geometries
  # still use the optimized-mode Gaussian component. Transition coverage changes
  # the randomized values, so check finite diagnostics and the actual nested
  # stopping rule instead of preserving the old proposal's numerical outputs.
  fixture <- readRDS(testthat::test_path("fixtures", "selection-certificate-laplace-reference.rds"))
  cases <- fixture[["cases"]]
  budget <- SELNORM_FACTOR_NODE_BUDGET
  for (name in c("cycle7", "cycle8")) {
    affordable <- sum(SELNORM_FACTOR_QUADRATURE_ORDERS^cases[[name]][["block"]][["rank"]] <= budget)
    expect_lt(affordable, 3L, label = name)
  }
  expect_false(.test_cert_probe(cases[["chain4"]][["block"]], 9L, refine = FALSE)[["has_moments"]])
  expect_false(.test_cert_probe(cases[["chain3"]][["block"]], 9L, refine = FALSE)[["has_moments"]])
  for (name in c("band_a", "band_b")) {
    probe <- .test_cert_probe(cases[[name]][["block"]], 95L, refine = FALSE)
    expect_true(probe[["has_moments"]], info = name)
    covariance <- matrix(probe[["moment_covariance"]], cases[[name]][["block"]][["rank"]])
    expect_error(chol(covariance), info = name)
  }
  run <- .test_cert_proposals_run
  for (name in names(cases)) for (seed in if (name %in% c("band_a", "band_b")) 1L else 1:2) {
    case <- cases[[name]]
    actual <- run(case, seed, 256L, 2048L)
    info <- paste(name, "seed", seed)
    if (case$block$rank == 1L) {
      expect_true(.test_cert_deterministic(actual), info = info)
      expect_identical(actual$relative_mcse, 0, info = info)
      expect_lte(abs(expm1(actual$log_normalizer - .test_cert_reference_parallel(case$block))),
        actual$relative_change + 1e-10, label = info)
      next
    }
    expect_false(.test_cert_deterministic(actual), info = info)
    expect_gt(actual[["relative_mcse"]], 0, label = info)
    expect_true(all(is.finite(unlist(actual[c("log_normalizer", "relative_mcse", "relative_change")]))),
                info = info)
    for (points in c(256L, 512L, 1024L, 2048L)) {
      stage <- run(case, seed, points, points)
      criterion <- max(.selection_qmc_z() * stage[["relative_mcse"]], stage[["relative_change"]])
      if (criterion <= .005) break
    }
    for (field in c("log_normalizer", "relative_mcse", "relative_change")) {
      expect_equal(actual[[field]], stage[[field]], tolerance = 1e-12, info = paste(info, field))
    }
  }
})
