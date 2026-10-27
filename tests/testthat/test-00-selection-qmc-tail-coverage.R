context("Selection QMC tail coverage")
skip_on_cran()

.test_qmc_tail_block <- function(loading, omega = c(1, .5), cuts = c(0, -Inf)) {

  scale <- sqrt(sum(loading^2))
  .test_cert_block(list(.test_cert_parallel_group(
    10L, loading, .04 * scale, -4 * scale)), omega, cuts)
}

.test_qmc_tail_native <- function(block, rules, kind = "factor", seed = 1L, sign = 1L) {

  rows <- block[["rows"]]
  design <- BayesTools::selection_qmc_design(
    dimensions = 2L * block[["rank"]], points = 8192L, scrambles = 8L, seed = seed)
  symbol <- paste0("RoBMA_selnorm_", kind, "_step_loglik_batch")
  .Call(symbol, rep(0, rows), matrix(sign * block[["mean"]], 1L),
    matrix(block[["sd"]], 1L), matrix(as.double(sign * block[["loading"]]), 1L),
    block[["sei"]], matrix(block[["omega"]], 1L), block[["z_lower"]], block[["z_upper"]],
    rep(1L, rows), as.integer(sign), TRUE, 1L,
    rules[["nodes"]], rules[["log_weights"]], as.double(rules[["orders"]]),
    as.double(design), 512L, 8192L, 8L, .005, TRUE, 0L, PACKAGE = "RoBMA")
}

test_that("the randomized fallback resolves a tail missed by every Gaussian proposal", {

  # Coincident rows reduce exactly to one standard normal coordinate, whatever
  # the number or orientation of the loading columns. Splitting around its
  # transition makes QUADPACK independent of the production proposal and of the
  # composite Gauss-Legendre reference in helper-selection-certificate.R.
  breaks <- c(-Inf, 0, 3.6, 4, 4.4, Inf)
  pieces <- lapply(seq_len(length(breaks) - 1L), function(index) {
    stats::integrate(function(u) stats::dnorm(u) *
      (.5 + .5 * stats::pnorm((u - 4) / .04))^10,
      breaks[[index]], breaks[[index + 1L]], rel.tol = 1e-12,
      abs.tol = 1e-15, subdivisions = 1000L)
  })
  mass <- sum(vapply(pieces, `[[`, numeric(1L), "value"))
  reference <- log(mass)
  expect_lt(sum(vapply(pieces, `[[`, numeric(1L), "abs.error")) / mass, 1e-10)
  expect_lt(abs(reference - -6.9048465796091056), 1e-12)
  expect_gt(1 - .5^10 / mass, .02)
  expect_lt(abs(.test_cert_reference_parallel(.test_qmc_tail_block(c(.6, .8))) - reference), 1e-11)

  factor_rules <- .selection_joint_factor_quadrature_rules()
  cases <- list(
    rank_two = list(loading = c(.6, .8), kind = "factor",
      rules = .selection_joint_cluster_quadrature_rules(c(3L, 5L, 7L))),
    cluster = list(loading = 1, kind = "cluster",
      rules = .selection_joint_cluster_quadrature_rules(c(15L, 31L))),
    rank_seven = list(loading = rep(1 / sqrt(7), 7L), kind = "factor", rules = factor_rules),
    rank_eight = list(loading = rep(1 / sqrt(8), 8L), kind = "factor", rules = factor_rules)
  )
  for (name in names(cases)) {
    case <- cases[[name]]
    for (seed in 1:2) {
      actual <- .test_qmc_tail_native(.test_qmc_tail_block(case[["loading"]]),
        case[["rules"]], case[["kind"]], seed = seed)
      info <- paste(name, "seed", seed)
      expect_identical(.test_cert_deterministic(actual), name == "cluster", info = info)
      expect_lt(abs(expm1(actual[["log_normalizer"]] - reference)), .005, label = info)
      # The old fallback reported about 3e-16 despite missing 2.6% of the mass.
      if (name == "cluster") {
        expect_identical(actual$relative_mcse, 0)
        expect_lte(abs(expm1(actual$log_normalizer - reference)), actual$relative_change + 1e-10,
          label = info)
      } else expect_gt(actual[["relative_mcse"]], 1e-10, label = info)
      expect_true(is.finite(actual[["relative_mcse"]]), info = info)
    }
  }
})

test_that("a rare intersection of distinct transition directions is accurate or rejected", {

  # Orthogonal loading rows give independent N(-4, 1 + .04^2) outcomes. Their
  # weighted probabilities therefore multiply exactly, although no loading row
  # is aligned with an input coordinate and the high-weight intersection is rare.
  loading <- qr.Q(qr(matrix(c(1, 2, 3, 4, 1, 2, 3, 2, 1), 3L)))
  block <- .test_cert_block(lapply(seq_len(3L), function(row) {
    .test_cert_parallel_group(1L, loading[row, ], .04, -4)
  }), c(1, 1e-4))
  reference <- 3 * log(1e-4 + (1 - 1e-4) * stats::pnorm(-4 / sqrt(1 + .04^2)))
  actual <- .test_qmc_tail_native(block, .selection_joint_factor_quadrature_rules(), seed = 2L)
  error <- abs(expm1(actual[["log_normalizer"]] - reference))
  if (error > .005) {
    expect_error(.selection_joint_factor_check_diagnostics(actual,
      list(relative_tolerance = .005), "factor"), "rejected by diagnostics")
  } else {
    expect_lte(error, .005)
  }
})

test_that("an almost unreachable low-weight bin cannot hide the ordinary selection tail", {

  ordinary <- .test_qmc_tail_block(c(.6, .8))
  block <- .test_qmc_tail_block(c(.6, .8), omega = c(1, .5, 1e-6),
    cuts = c(0, -100, -Inf))
  reference <- .test_cert_reference_parallel(ordinary)
  # Changing a row's weight below -100 can alter the product by at most the
  # probability that any row enters that bin. The union bound is negligible
  # relative to the ordinary normalizer, so its independent reference applies.
  log_omitted_bound <- log(block[["rows"]]) +
    stats::pnorm((-100 + 4) / sqrt(1 + .04^2), log.p = TRUE)
  expect_lt(log_omitted_bound - reference, log(1e-12))
  rules <- .selection_joint_cluster_quadrature_rules(c(3L, 5L, 7L))
  actual <- .test_qmc_tail_native(block, rules)
  error <- abs(expm1(actual[["log_normalizer"]] - reference))
  if (error > .005) {
    expect_error(.selection_joint_factor_check_diagnostics(actual,
      list(relative_tolerance = .005), "factor"), "rejected by diagnostics")
  } else {
    expect_lte(error, .005)
  }
})

test_that("deeper simultaneous tails cannot pass with an underestimated integration error", {

  # A dense Householder reflection preserves four independent standard normals.
  # At a five-SD transition, matching the low weight to the marginal tail keeps
  # simultaneous high-weight intersections important to the product normalizer.
  # This regression was identified during independent stress checks.
  direction <- seq_len(4L)
  loading <- diag(4L) - 2 * tcrossprod(direction) / sum(direction^2)
  probability <- stats::pnorm(-5 / sqrt(diag(tcrossprod(loading)) + .04^2))
  low_weight <- mean(probability) / .3
  block <- .test_cert_block(lapply(seq_len(4L), function(row) {
    .test_cert_parallel_group(1L, loading[row, ], .04, -5)
  }), c(1, low_weight))
  reference <- sum(log(low_weight + (1 - low_weight) * probability))
  actual <- .test_qmc_tail_native(block, .selection_joint_factor_quadrature_rules(), seed = 13L)
  error <- abs(expm1(actual[["log_normalizer"]] - reference))
  if (error > .005) {
    expect_error(.selection_joint_factor_check_diagnostics(actual,
      list(relative_tolerance = .005), "factor"), "rejected by diagnostics")
  } else {
    expect_lte(error, .005)
  }
})

test_that("tail proposals preserve nonmonotone weights, their scale and effect direction", {

  # A distant narrow elevated band has two transitions with opposite jumps.
  block <- .test_qmc_tail_block(c(.6, .8), omega = c(.5, 1, .5), cuts = c(.2, 0, -Inf))
  rules <- .selection_joint_cluster_quadrature_rules(c(3L, 5L, 7L))
  reference <- .test_cert_reference_parallel(block)
  actual <- .test_qmc_tail_native(block, rules)
  expect_false(.test_cert_deterministic(actual))
  expect_lt(abs(expm1(actual[["log_normalizer"]] - reference)), .005)
  reflected <- .test_qmc_tail_native(block, rules, sign = -1L)
  # Reflection can change the fixed QMC coordinates. Both estimates target the
  # same independent integral, and the Gaussian numerator is invariant exactly.
  expect_lt(abs(expm1(reflected[["log_normalizer"]] - reference)), .005)
  expect_lt(abs(sum(unlist(reflected[c("log_density", "log_normalizer")])) -
    sum(unlist(actual[c("log_density", "log_normalizer")]))), 1e-10)
  for (scale in c(3, 1e-20, 1e20)) {
    scaled <- block
    scaled[["omega"]] <- scale * block[["omega"]]
    rescaled <- .test_qmc_tail_native(scaled, rules)
    # Extreme scaling shifts log targets by hundreds; numerical optimization
    # of the Gaussian proposals then incurs small finite-difference roundoff.
    tolerance <- if (scale == 3) 1e-10 else 1e-8
    expect_lt(abs(rescaled[["log_normalizer"]] - actual[["log_normalizer"]] -
      block[["rows"]] * log(scale)), tolerance, label = scale)
    expect_lt(abs(rescaled[["log_density"]] - actual[["log_density"]]), tolerance, label = scale)
    expect_lt(abs(rescaled[["relative_mcse"]] - actual[["relative_mcse"]]), tolerance, label = scale)
    expect_lt(abs(rescaled[["relative_change"]] - actual[["relative_change"]]), tolerance, label = scale)
  }
})

test_that("zero weights and rows without latent loadings retain their probability factors", {

  block <- .test_qmc_tail_block(c(.6, .8), omega = c(1, 0))
  rules <- .selection_joint_cluster_quadrature_rules(c(3L, 5L, 7L))
  reference <- .test_cert_reference_parallel(block)
  actual <- .test_qmc_tail_native(block, rules)
  expect_lt(abs(expm1(actual[["log_normalizer"]] - reference)), .005)

  # An independent row with no latent loading contributes its exact selection
  # probability. It introduces no transition direction for the proposal.
  extended <- .test_cert_block(c(block[["groups"]], list(
    .test_cert_parallel_group(1L, c(0, 0), .7, -.2))), block[["omega"]])
  extra <- .test_qmc_tail_native(extended, rules)
  expected <- reference + stats::pnorm(-.2 / .7, log.p = TRUE)
  expect_lt(abs(expm1(extra[["log_normalizer"]] - expected)), .005)

  inactive <- .test_qmc_tail_block(c(.6, .8, 0))
  inactive_result <- .test_qmc_tail_native(inactive, rules)
  expect_lt(abs(expm1(inactive_result[["log_normalizer"]] -
    .test_cert_reference_parallel(inactive))), .005)
})

test_that("an unrepresentable event probability cannot silently become zero likelihood", {

  ordinary <- .test_cert_block(list(
    .test_cert_parallel_group(10L, c(1, 0), .04, -4)), c(1, .5))
  block <- .test_cert_block(c(ordinary[["groups"]], list(
    .test_cert_parallel_group(1L, c(0, 1), .04, -1000))), c(1, .5))
  # The two latent coordinates are independent. The distant row's positive
  # minimum weight keeps its factor near .5 even though its excess-bin mass
  # cannot be represented relative to the ordinary transition's event mass.
  reference <- .test_cert_reference_parallel(ordinary) +
    log(.5 + .5 * stats::pnorm(-1000 / sqrt(1 + .04^2)))
  expect_gt(reference, block[["rows"]] * log(.5))
  rules <- .selection_joint_cluster_quadrature_rules(c(3L, 5L, 7L))
  actual <- .test_qmc_tail_native(block, rules)
  # The distant event has no sampling probability, but its row stays in the
  # target. The fallback must neither reject the block nor lose the factor.
  expect_true(all(is.finite(unlist(actual[c("log_density", "log_normalizer",
    "relative_mcse", "relative_change")]))))
  expect_lte(abs(expm1(actual[["log_normalizer"]] - reference)), .005)
  expect_error(.selection_joint_factor_check_diagnostics(actual,
    list(relative_tolerance = .005), "factor"), NA)
})

test_that("a negligible event component is dropped from the proposal without rejecting the block", {

  # Ten coincident rows put a transition at one or two SD on the first
  # coordinate. One row on the independent second coordinate sits so far below
  # its transition that the relative mass of its event underflows even long
  # double (about exp(-11400) at -160, exactly zero probability at -1e200), so
  # it is never sampled. The target factor of that row is still
  # .5 + .5 pnorm(far / sqrt(1 + .04^2)), the product of independent
  # coordinates gives the reference, and every run is accepted.
  rules <- .selection_joint_cluster_quadrature_rules(c(3L, 5L, 7L))
  for (transition in c(1, 2)) {
    ordinary <- .test_cert_block(list(
      .test_cert_parallel_group(10L, c(1, 0), .04, -transition)), c(1, .5))
    for (far in c(-150, -160, -1000, -1e200)) {
      block <- .test_cert_block(c(ordinary[["groups"]], list(
        .test_cert_parallel_group(1L, c(0, 1), .04, far))), c(1, .5))
      reference <- .test_cert_reference_parallel(ordinary) +
        log(.5 + .5 * stats::pnorm(far / sqrt(1 + .04^2)))
      for (seed in 1:2) {
        info <- paste("transition", transition, "far", far, "seed", seed)
        actual <- .test_qmc_tail_native(block, rules, seed = seed)
        expect_false(.test_cert_deterministic(actual), info = info)
        # The data density of a row at -1e200 is itself exactly zero; only the
        # normalizer and its diagnostics must be finite.
        expect_true(all(is.finite(unlist(actual[c("log_normalizer",
          "relative_mcse", "relative_change")]))), info = info)
        expect_lte(abs(expm1(actual[["log_normalizer"]] - reference)), .005, label = info)
        expect_error(.selection_joint_factor_check_diagnostics(actual,
          list(relative_tolerance = .005), "factor"), NA, info = info)
      }
    }
  }
})

test_that("an internal proposal boundary in a supplied QMC design has finite weights", {

  block <- .test_cert_block(list(
    .test_cert_parallel_group(1L, c(1, 0), .04, -4),
    .test_cert_parallel_group(1L, c(0, 1), .04, -4)), c(1, .5))
  # Two equally weighted transition directions put the categorical boundary at
  # .5. The origin and their two structural cutoffs also give three locations
  # per coordinate in the boundary proposal. Deliberately identical designs
  # test valid-coordinate boundaries, not Monte Carlo calibration; the short
  # ladder forces the fallback.
  rules <- .selection_joint_cluster_quadrature_rules(c(1L, 2L, 3L))
  for (uniform in c(.5, 1 / 3)) {
    design <- array(uniform, dim = c(8L, 16L, 4L))
    actual <- .test_cert_native_factor(block, initial = 16L, points = 16L,
      rules = rules, design = design)
    expect_false(.test_cert_deterministic(actual), info = as.character(uniform))
    expect_true(all(is.finite(unlist(actual[c("log_density", "log_normalizer",
      "relative_mcse", "relative_change")]))), info = as.character(uniform))
  }
})
