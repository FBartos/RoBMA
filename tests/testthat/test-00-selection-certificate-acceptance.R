context("Selection certificate acceptance")
skip_on_cran()

# The product deterministic ladder of the factor kernel accepts the first rung
# that is finite, has both successive changes within the tolerance, passes the
# tail gate, and whose shift-amplitude certificate and combined allowance fit
# the tolerance. The tests compare the native kernel with a plain-R emulation of
# that ladder (tensor Gauss-Hermite rungs, independent of the native code) and,
# for the accepted rung, with an independent integral.

test_that("random two-group blocks follow the R certificate ladder at default and relaxed tolerances", {

  set.seed(4711)
  routes <- character()
  for (i in seq_len(8L)) {
    group <- function() list(count = sample(1:8, 1L), loading = stats::rnorm(2L, 0, .7),
      sd = stats::runif(1L, .5, 1), mean = stats::rnorm(1L, 0, 1.5), sei = 1)
    block <- .test_cert_block(list(group(), group()), c(1, stats::runif(1L, 0, .8)))
    reference <- if (i <= 3L) .test_cert_reference_two_groups(block)
    for (tolerance in c(.005, .02, .001)) {
      actual   <- .test_cert_native_factor(block, tolerance = tolerance)
      expected <- .test_cert_emulate_ladder(block, tolerance)
      info     <- paste(i, tolerance)
      if (is.null(expected)) {
        routes <- c(routes, "fallback")
        expect_gt(actual[["relative_mcse"]], 0, label = info)
      } else {
        routes <- c(routes, "deterministic")
        expect_identical(actual[["relative_mcse"]], 0, label = info)
        expect_equal(actual[["log_normalizer"]], expected[["value"]], tolerance = 1e-10, info = info)
        # the accepted rung is within the (non-default) tolerance of an independent integral
        if (!is.null(reference)) {
          expect_lt(.test_cert_relative_error(actual[["log_normalizer"]], reference), tolerance, label = info)
        }
      }
    }
  }
  expect_true("deterministic" %in% routes)
})

test_that("the combined exterior and amplitude allowance decides at its exact boundary", {

  # 20 rows selected deep in a tail: at order 9 the exterior bound and the
  # amplitude are comparable, so each separate gate accepts the rung over a
  # range of tolerances that the combined allowance refuses. Neither value
  # depends on the tolerance.
  block <- .test_cert_block(list(
    list(count = 20L, loading = c(.4, .8), sd = .3, mean = -3, sei = 1)), c(1, 0))
  probe <- .test_cert_probe(block, 9L, tolerance = 1)
  tail <- probe[["tail_bound"]]
  amplitude <- probe[["amplitude_sum"]]
  allowance <- tail + .test_cert_allowance_factor * amplitude
  threshold <- max(2 * amplitude, allowance / (1 - allowance))
  expect_lt(allowance, 1)
  # the combined allowance is the binding condition, not the half-tolerance bound
  expect_gt(allowance / (1 - allowance), 2 * amplitude)
  expect_gt(tail, 0.7 * amplitude)
  in_between <- mean(c(max(2 * amplitude, tail), threshold))
  expect_true(2 * amplitude <= in_between && tail <= in_between && in_between < threshold)
  expect_false(.test_cert_probe(block, 9L, tolerance = in_between)[["accepted"]])
  for (side in c(-1, 1)) {
    tolerance <- threshold * (1 + side * 1e-6)
    expect_identical(.test_cert_probe(block, 9L, tolerance = tolerance)[["accepted"]], side > 0,
                     info = side)
    # the kernel moves up the ladder exactly when the rung is refused
    orders <- SELNORM_FACTOR_QUADRATURE_ORDERS
    values <- vapply(orders, function(order) .test_cert_probe(block, order)[["log_normalizer"]], numeric(1L))
    actual <- .test_cert_native_factor(block, tolerance = tolerance)
    accepted <- orders[which.min(abs(values - actual[["log_normalizer"]]))]
    expect_identical(actual[["relative_mcse"]], 0, info = side)
    expect_identical(accepted, if (side > 0) 9L else 11L, info = side)
  }
})

test_that("a block with a constant selection weight is accepted without a certificate", {

  # No transition: no aliasing to certify. Equal weights leave only the exterior
  # mass, held to the combined allowance.
  block  <- .test_cert_block(list(list(count = 4L, loading = c(.6, .5), sd = 1.1, mean = .2, sei = 1)),
                             c(1, 1))
  probe  <- .test_cert_probe(block, 11L)
  actual <- .test_cert_native_factor(block)
  expect_false(probe[["transitions"]])
  expect_identical(actual[["relative_mcse"]], 0)
  # four rows of constant weight one: the normalizer is exactly one
  expect_equal(actual[["log_normalizer"]], 0, tolerance = 1e-12)
  # the first rung whose exterior mass fits the allowance (order 7)
  expect_equal(.test_cert_emulate_ladder(block)[["order"]], 7L)
  expect_equal(actual[["relative_change"]], .test_cert_emulate_ladder(block)[["change"]], tolerance = 1e-8)
})

test_that("the cluster kernel follows the one-change R certificate ladder", {

  set.seed(99)
  routes <- character()
  rules <- SELNORM_CLUSTER_QUADRATURE_ORDERS
  for (i in seq_len(6L)) {
    count <- sample(c(3L, 12L, 40L), 1L)
    block <- .test_cert_block(list(list(count = count, loading = stats::runif(1L, .4, 1.6),
      sd = stats::runif(1L, .3, 1), mean = stats::rnorm(1L, 0, 2), sei = 1)),
      c(1, stats::runif(1L, 0, .6)))
    for (tolerance in c(.005, .02)) {
      actual   <- .test_cert_native_cluster(block, tolerance = tolerance)
      expected <- .test_cert_emulate_ladder(block, tolerance, orders = rules, changes = 1L)
      info     <- paste(i, tolerance)
      if (is.null(expected)) {
        routes <- c(routes, "fallback")
        expect_gt(actual[["relative_mcse"]], 0, label = info)
      } else {
        routes <- c(routes, "deterministic")
        expect_identical(actual[["relative_mcse"]], 0, label = info)
        expect_equal(actual[["log_normalizer"]], expected[["value"]], tolerance = 1e-10, info = info)
        expect_equal(actual[["relative_change"]], expected[["change"]], tolerance = 1e-8, info = info)
        expect_lt(.test_cert_relative_error(actual[["log_normalizer"]],
          .test_cert_reference_parallel(block)), tolerance, label = info)
      }
    }
  }
  expect_true("deterministic" %in% routes)
})

test_that("the best weight rule keeps its geometric acceptance, decisions and values", {

  # Seeded random best-rule blocks (one- and two-sided; rank-one clusters and
  # rank-two and rank-three factors) with the outputs of the build before the
  # certificate: the best rule keeps the precheck, the exact ladder start and
  # the geometric acceptance, so decisions and values are unchanged.
  reference <- readRDS(testthat::test_path("fixtures", "selection-certificate-best-reference.rds"))
  factor_rules  <- .selection_joint_factor_quadrature_rules()
  cluster_rules <- .selection_joint_cluster_quadrature_rules(SELNORM_CLUSTER_QUADRATURE_ORDERS)
  routes <- character()
  for (i in seq_along(reference[["blocks"]])) {
    b <- reference[["blocks"]][[i]]
    actual <- if (b[["rank"]] == 1L) {
      qmc <- BayesTools::selection_qmc_design(2L, 512L, 8L, seed = 173L)
      .Call("RoBMA_selnorm_cluster_step_loglik_batch", b[["yi"]], matrix(b[["mean"]], 1L),
        matrix(b[["sd"]], 1L), matrix(b[["loading"]][, 1L], 1L), b[["sei"]], matrix(b[["omega"]], 1L),
        b[["z_lower"]], b[["z_upper"]], b[["obs_bin"]], 1L, TRUE, 1L, cluster_rules[["nodes"]],
        cluster_rules[["log_weights"]], as.double(cluster_rules[["orders"]]), as.double(qmc),
        64L, 512L, 8L, .005, TRUE, as.integer(b[["rule"]]), PACKAGE = "RoBMA")
    } else {
      qmc <- BayesTools::selection_qmc_design(2L * b[["rank"]], 512L, 8L, seed = 173L)
      .Call("RoBMA_selnorm_factor_step_loglik_batch", b[["yi"]], matrix(b[["mean"]], 1L),
        matrix(b[["sd"]], 1L), matrix(as.double(b[["loading"]]), 1L), b[["sei"]], matrix(b[["omega"]], 1L),
        b[["z_lower"]], b[["z_upper"]], b[["obs_bin"]], 1L, TRUE, 1L, factor_rules[["nodes"]],
        factor_rules[["log_weights"]], as.double(factor_rules[["orders"]]), as.double(qmc),
        64L, 512L, 8L, .005, TRUE, as.integer(b[["rule"]]), PACKAGE = "RoBMA")
    }
    expected <- reference[["results"]][[i]]
    info <- paste("block", i)
    expect_identical(actual[["relative_mcse"]] == 0, expected[["relative_mcse"]] == 0, info = info)
    routes <- c(routes, if (actual[["relative_mcse"]] == 0) "deterministic" else "fallback")
    for (field in c("log_density", "relative_mcse", "relative_change", "log_normalizer")) {
      expect_equal(as.numeric(actual[[field]]), expected[[field]], tolerance = 1e-12, info = paste(info, field))
    }
  }
  expect_true(all(c("deterministic", "fallback") %in% routes))
})
