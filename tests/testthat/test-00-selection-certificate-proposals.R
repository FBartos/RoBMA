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
    expect_true(probe[["has_moments"]], info = order)
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

# Native fallback runs. Every block here is certified by neither the ladder nor the
# refinement tier (the clusters have no tier: by no rung of their ladder), so it
# reaches the randomized fallback; the fixed design seeds make the runs repeatable.
# The parameters come from the research battery (seeds 41 and 21) and random searches
# over forests (a shared axis and two leaves) and rank-two partial overlaps (tensor
# rule). The optimized-mode proposals of the previous build accept a minority of
# them and err by 10-1000 MCSE on several. Blocks whose groups all load on both axes
# (parallel_80, band_125, crossing_20) are forests of depth one: a root that owns no
# rows and a leaf.
.test_cert_proposal_cases <- function() {

  g <- function(count, loading, sd, mean, sei = 1) list(count = count, loading = loading, sd = sd, mean = mean, sei = sei)
  list(
    parallel_80 = list(kind = "factor", reference = .test_cert_reference_parallel, block = .test_cert_block(list(
      g(12L, c(0.014692097199723289, -0.8616553723389665), 0.13653164601012455, -5.6000270018551461, 0.34873594478420489)),
      c(1, 0.050000000000000003), c(0, -Inf))),
    band_125 = list(kind = "factor", reference = .test_cert_reference_parallel, block = .test_cert_block(list(
      g(44L, c(1.4323285021004943, 0.083058174726865039), 0.060454506795780072, 1.8231652047112927, 0.55342308823237685)),
      c(1, 1.1010918456055847, 1), c(1, -1, -Inf))),
    crossing_20 = list(kind = "factor", reference = .test_cert_reference_two_groups, block = .test_cert_block(list(
      g(77L, c(-0.91463423771209651, 0.78986739681623297), 0.28436590998117955, -2.3718989917066349),
      g(18L, c(0.015058325926115159, -0.30020325671357212), 0.019704307952072554, 1.8604332933661574)),
      c(1, 0.050000000000000003), c(1.96, -Inf))),
    # partial overlaps (a group on both axes, one group on each axis alone): the tensor rule
    tensor_45 = list(kind = "factor", reference = .test_cert_reference_rank_two, pinned = -13.835385492311222,
      block = .test_cert_block(list(
        g(27L, c(0.82706647334157291, 0.65077649884120636), 0.054795022412045077, -0.043517232008220219),
        g(11L, c(-1.3674417541470836, 0), 0.064991112511571233, 1.1165970716872602),
        g(22L, c(0, 1.3796914536846054), 0.057839652420356456, 0.27615990856989825)),
        c(1, 0.20000000000000001), c(1.96, -Inf))),
    tensor_77 = list(kind = "factor", reference = .test_cert_reference_rank_two, pinned = -2.5592954958880094,
      block = .test_cert_block(list(
        g(12L, c(-0.34458291062427782, -0.83820528250345905), 0.41598330888330703, 1.3861580619958969),
        g(25L, c(-0.69611133285596238, 0), 0.081507582329821651, 2.515650777258446),
        g(17L, c(0, -1.8009762478430924), 0.13199939157250279, 0.82995553693087865)),
        c(1, 0.01), c(1.96, -Inf))),
    forest_74 = list(kind = "factor", reference = .test_cert_reference_forest, block = .test_cert_block(list(
      g(65L, c(0.060605698492417107, 0.9655385554788849, 0), 0.57052084244334667, -2.3719209529207341),
      g(56L, c(0.060605698492417107, 0, 1.6252115421914184), 0.040203888654504902, -4.6092623941950031)),
      c(1, 0.20000000000000001), c(0, -Inf))),
    forest_84 = list(kind = "factor", reference = .test_cert_reference_forest, block = .test_cert_block(list(
      g(127L, c(0.12399440036971428, 1.9003074081263616, 0), 0.0989997504638504, -2.5188992921523741),
      g(30L, c(0.12399440036971428, 0, 1.1540407331499025), 0.030234154411339583, 1.6193386859839336)),
      c(1, 0.050000000000000003), c(1.645, -Inf))),
    # a cluster whose Laplace proposal sees only the bulk of the prior: all eight scrambles
    # then agree, the MCSE is zero and the estimate sits 7% low
    cluster_733 = list(kind = "cluster", reference = .test_cert_reference_parallel, block = .test_cert_block(list(
      g(4L, 1.8707955473336113, 0.042233265807645956, -6.8303662251056725, 0.22565989842880826)),
      c(1, 0.20000000000000001), c(0, -Inf))),
    cluster_189 = list(kind = "cluster", reference = .test_cert_reference_parallel, block = .test_cert_block(list(
      g(100L, 0.65477525608481191, 0.026923363058843731, -2.5196777242881367, 0.20064200965816639)),
      c(1, 0.050000000000000003), c(1.645, -Inf)))
  )
}

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
  # those of that build, so every output agrees.
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
  for (name in names(cases)) for (seed in if (name %in% c("band_a", "band_b")) 1L else 1:2) {
    case <- cases[[name]]
    actual <- if (case[["kind"]] == "factor") {
      .test_cert_native_factor(case[["block"]], seed = seed, initial = 256L, points = 2048L)
    } else {
      rules <- .selection_joint_cluster_quadrature_rules(c(15L, 31L))
      block <- case[["block"]]
      K <- block[["rows"]]
      qmc <- BayesTools::selection_qmc_design(dimensions = 2L, points = 2048L, scrambles = 8L, seed = seed)
      .Call("RoBMA_selnorm_cluster_step_loglik_batch", rep(0, K), matrix(block[["mean"]], 1L),
        matrix(block[["sd"]], 1L), matrix(as.double(block[["loading"]]), 1L), block[["sei"]],
        matrix(block[["omega"]], 1L), block[["z_lower"]], block[["z_upper"]], rep(1L, K), 1L, TRUE, 1L,
        rules[["nodes"]], rules[["log_weights"]], as.double(rules[["orders"]]), as.double(qmc), 256L,
        2048L, 8L, .005, TRUE, 0L, PACKAGE = "RoBMA")
    }
    expected <- reference[["results"]][[paste(name, seed)]]
    info <- paste(name, "seed", seed)
    expect_gt(actual[["relative_mcse"]], 0, label = info)
    for (field in c("log_normalizer", "relative_mcse", "relative_change")) {
      expect_equal(as.numeric(actual[[field]]), as.numeric(expected[[field]]), tolerance = 1e-12,
                   info = paste(info, field))
    }
  }
})
