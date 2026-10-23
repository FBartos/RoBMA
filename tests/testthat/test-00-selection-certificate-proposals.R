context("Selection certificate QMC proposals")
skip_on_cran()

# The randomized fallback of the product rule proposes from the prior, a Gaussian
# with the tilted moments of the finest deterministic rung (mean and covariance
# of the latent coordinates under the selected integrand) and a copy with 1.5
# times their standard deviations. The moments are used for tensors, forests of
# depth at most one and rank-one clusters when at least three rungs were finite
# and the covariance is positive definite; every other block keeps the
# optimized-mode proposals. A finest rung that misses the exterior mass gives no
# moments either (test-00-selection-factor-tail-coverage.R keeps that case). The
# stopping rule and its diagnostics are unchanged.

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

test_that("the moment proposals calibrate the randomized fallback on deferred geometries", {

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
      # the randomized fallback ran, with a positive spread across the scrambles
      expect_gt(actual[["relative_mcse"]], 0, label = info)
      error <- abs(expm1(actual[["log_normalizer"]] - reference))
      expect_lte(error, 3 * actual[["relative_mcse"]], label = info)
    }
  }
})

test_that("blocks without usable moments keep the optimized-mode proposals and their outputs", {

  # Outputs of robma-claude-final7 (the build before the certificate) for blocks
  # that reach the fallback there and, here, with fewer than three finite rungs
  # (rank seven and eight tensors: two and one affordable rungs; a rank-one
  # cluster with a two-rung ladder), with a forest of depth two or three, or
  # with a moment covariance that is not positive definite. The proposals are
  # those of that build, so every output agrees, except that a run whose stored
  # stage fails the scaled stopping rule (2.36 times the standard error) goes on to a
  # later stage: the stored output is then the stage where the old rule stopped
  # (.test_cert_expect_stored_fallback).
  reference <- readRDS(testthat::test_path("fixtures", "selection-certificate-laplace-reference.rds"))
  cases <- reference[["cases"]]
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
  run <- function(case, seed, initial, points) {
    if (case[["kind"]] == "factor") {
      .test_cert_native_factor(case[["block"]], seed = seed, initial = initial, points = points)
    } else {
      rules <- .selection_joint_cluster_quadrature_rules(c(15L, 31L))
      block <- case[["block"]]
      K <- block[["rows"]]
      qmc <- BayesTools::selection_qmc_design(dimensions = 2L, points = points, scrambles = 8L, seed = seed)
      .Call("RoBMA_selnorm_cluster_step_loglik_batch", rep(0, K), matrix(block[["mean"]], 1L),
        matrix(block[["sd"]], 1L), matrix(as.double(block[["loading"]]), 1L), block[["sei"]],
        matrix(block[["omega"]], 1L), block[["z_lower"]], block[["z_upper"]], rep(1L, K), 1L, TRUE, 1L,
        rules[["nodes"]], rules[["log_weights"]], as.double(rules[["orders"]]), as.double(qmc),
        as.integer(initial), as.integer(points), 8L, .005, TRUE, 0L, PACKAGE = "RoBMA")
    }
  }
  for (name in names(cases)) for (seed in if (name %in% c("band_a", "band_b")) 1L else 1:2) {
    case <- cases[[name]]
    actual <- run(case, seed, 256L, 2048L)
    expected <- reference[["results"]][[paste(name, seed)]]
    info <- paste(name, "seed", seed)
    expect_gt(actual[["relative_mcse"]], 0, label = info)
    .test_cert_expect_stored_fallback(actual, expected,
      function(points) run(case, seed, points, points),
      stages = c(256L, 512L, 1024L, 2048L), info = info)
  }
})
