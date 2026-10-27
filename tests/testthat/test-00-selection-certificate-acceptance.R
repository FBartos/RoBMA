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
                     info = as.character(side))
    # The probe exposes the raw rung criterion; the full ladder also caps a
    # loose request at the tolerance of its calibration.
    orders <- SELNORM_FACTOR_QUADRATURE_ORDERS
    values <- vapply(orders, function(order) .test_cert_probe(block, order)[["log_normalizer"]], numeric(1L))
    actual <- .test_cert_native_factor(block, tolerance = tolerance)
    accepted <- orders[which.min(abs(values - actual[["log_normalizer"]]))]
    expect_identical(actual[["relative_mcse"]], 0, info = as.character(side))
    expected <- .test_cert_emulate_ladder(block, tolerance)
    expect_identical(accepted, expected[["order"]], info = as.character(side))
  }
})

test_that("the allowance boundary decides at the calibrated cap of loose requests", {

  # A request above .005 is held to .005 when the weight has transitions, so
  # the kernel's allowance rule has its boundary at a threshold of .005 for
  # every such request; the test above sits far from it. Twenty rows with a
  # latent mean of -1 and the residual SD as the free parameter: the combined
  # allowance of order 31 falls through .005 between SDs of .275 and .3, and
  # the other gates of that rung hold on both sides. The threshold belongs to
  # the block, not to the tolerance.
  order <- 31L
  block_at <- function(sd) .test_cert_block(list(
    list(count = 20L, loading = c(.4, .8), sd = sd, mean = -1, sei = 1)), c(1, 0))
  threshold <- function(sd) {
    probe <- .test_cert_probe(block_at(sd), order, tolerance = 1)
    allowance <- probe[["tail_bound"]] + .test_cert_allowance_factor * probe[["amplitude_sum"]]
    max(2 * probe[["amplitude_sum"]], allowance / (1 - allowance))
  }
  expect_gt(threshold(.275), .005)
  expect_lt(threshold(.3), .005)
  orders <- SELNORM_FACTOR_QUADRATURE_ORDERS
  for (side in c(1, -1)) {
    # side 1: the threshold just inside the cap; -1: just outside
    sd <- stats::uniroot(function(x) log(threshold(x)) - log(.005 * (1 - side * 1e-6)),
                         c(.275, .3), tol = 1e-13)$root
    block <- block_at(sd)
    expect_lt(abs(threshold(sd) / .005 - (1 - side * 1e-6)), 1e-9)
    values <- vapply(orders, function(o) .test_cert_probe(block, o)[["log_normalizer"]], numeric(1L))
    strict <- .test_cert_native_factor(block, tolerance = .005)
    for (tolerance in c(.005, .02, .05)) {
      info   <- paste(side, tolerance)
      actual <- .test_cert_native_factor(block, tolerance = tolerance)
      accepted <- orders[which.min(abs(values - actual[["log_normalizer"]]))]
      expect_identical(actual[["relative_mcse"]], 0, info = info)
      # the rung is kept inside the cap; outside it the kernel moves up
      expect_identical(accepted, if (side > 0) order else 43L, info = info)
      expect_identical(accepted, .test_cert_emulate_ladder(block, tolerance)[["order"]], info = info)
      expect_identical(actual[["log_normalizer"]], strict[["log_normalizer"]], info = info)
      # the raw criterion of a loose request accepts the rung on both sides:
      # outside, the cap is what refuses it
      expect_identical(.test_cert_probe(block, order, tolerance = tolerance)[["accepted"]],
                       side > 0 || tolerance > .005, info = info)
    }
  }
})

test_that("a block with a constant selection weight is accepted without a certificate", {

  # No transition: no aliasing to certify. Equal weights leave only the exterior
  # mass, held to the combined allowance tail <= tolerance / (1 + tolerance).
  block  <- .test_cert_block(list(list(count = 4L, loading = c(.6, .5), sd = 1.1, mean = .2, sei = 1)),
                             c(1, 1))
  orders <- SELNORM_FACTOR_QUADRATURE_ORDERS
  probes <- lapply(orders, function(order) .test_cert_probe(block, order, refine = FALSE))
  values <- vapply(probes, function(probe) probe[["log_normalizer"]], numeric(1L))
  tails  <- vapply(probes, function(probe) probe[["tail_bound"]], numeric(1L))
  expect_false(probes[[1L]][["transitions"]])
  # The rung the kernel returned is the one whose value it carries: the values of
  # constant-weight rungs differ by rounding only, and each is distinct.
  native_rung <- function(actual) {
    rung <- which(values == actual[["log_normalizer"]])
    expect_length(rung, 1L)
    rung
  }

  actual <- .test_cert_native_factor(block)
  expect_identical(actual[["relative_mcse"]], 0)
  # four rows of constant weight one: the normalizer is exactly one
  expect_equal(actual[["log_normalizer"]], 0, tolerance = 1e-12)
  # the first rung whose exterior mass fits the allowance (order 7, the first rung
  # with two successive changes)
  rung <- native_rung(actual)
  expect_identical(orders[rung], 7L)
  expect_identical(orders[which(tails <= .005 / 1.005)[1L]], 7L)
  # the change is that of the rung and the one before it (rounding noise: another
  # rung would show another value, which the absolute tolerance resolves)
  expect_equal(actual[["relative_change"]], abs(expm1(values[rung - 1L] - values[rung])),
               tolerance = 1e-20)
  # the allowance decides at its exact boundary tail = tolerance / (1 + tolerance),
  # which the plain exterior gate (tail <= tolerance) would pass
  threshold <- tails[rung] / (1 - tails[rung])
  for (side in c(-1, 1)) {
    boundary <- .test_cert_native_factor(block, tolerance = threshold * (1 + side * 1e-6))
    expect_identical(boundary[["relative_mcse"]], 0, info = as.character(side))
    expect_identical(orders[native_rung(boundary)], if (side > 0) 7L else 9L, info = as.character(side))
  }
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

test_that("loose requests do not accept unresolved product quadrature walls", {

  # These opposing walls were accepted with 5.09% error at tolerance .02 and
  # 99.998% error at .05. References integrate the one-dimensional Gaussian
  # expectation independently of the native Gauss-Hermite ladder.
  group <- function(count, loading, mean, sd) {
    list(count = count, loading = loading, sd = sd, mean = mean, sei = 1)
  }
  blocks <- list(
    .test_cert_block(list(
      group(21L, -.331695895502344, 1.41240275628494, .105668749149233),
      group(11L, 1.65803034252021, -7.03005380219133, .105668749149233)
    ), c(1, .0897774216356852)),
    .test_cert_block(list(
      group(10L, 1.37901563644409, -4.65666543654929, .0908403173702142),
      group(30L, -1.62369661589619, 9.74595427356444, .0908403173702142),
      group(18L, -.833751429128460, 5.95677370447769, .0908403173702142),
      group(19L, -.374783916166052, 3.75651554236848, .0908403173702142),
      group(31L, -1.83758793340530, 10.7586220994678, .0908403173702142),
      group(20L, -1.85229506243486, 10.8354220438896, .0908403173702142)
    ), c(1, .001090973524444), c(1.96, -Inf))
  )

  for (block in blocks) {
    reference <- .test_cert_reference_parallel(block)
    fast <- !is.null(.test_cert_emulate_ladder(block,
      orders = SELNORM_CLUSTER_QUADRATURE_ORDERS, changes = 1L)) ||
      .test_cert_tier(block, finest = 1023L, previous = 511L, cluster = TRUE)$accepted
    for (evaluate in list(.test_cert_native_cluster)) {
      strict <- evaluate(block, tolerance = .005)
      for (tolerance in c(.005, .02, .05)) {
        value <- evaluate(block, tolerance = tolerance)
        expect_silent(.selection_joint_factor_check_diagnostics(
          value, list(relative_tolerance = tolerance), "test"))
        expect_lt(.test_cert_relative_error(value[["log_normalizer"]], reference), tolerance)
        if (fast) {
          expect_identical(value[["log_normalizer"]], strict[["log_normalizer"]])
        } else {
          expect_true(.test_cert_deterministic(value))
          expect_identical(value$relative_mcse, 0)
          expect_lte(.test_cert_relative_error(value$log_normalizer, reference),
            value$relative_change + 1e-10)
        }
      }
    }
  }

  # The same calibration applies to tensor and forest product ladders, not
  # just the scalar cluster entry point.
  for (block in .test_cert_small_blocks()[-1L]) {
    strict <- .test_cert_native_factor(block)
    for (tolerance in c(.02, .05)) {
      loose <- .test_cert_native_factor(block, tolerance = tolerance)
      expect_identical(attr(loose, "deterministic"), attr(strict, "deterministic"))
      if (isTRUE(as.logical(attr(strict, "deterministic")))) {
        expect_identical(loose[["log_normalizer"]], strict[["log_normalizer"]])
      }
    }
  }

  # With deliberately inadequate QMC budgets, the second wall now reaches the
  # normal diagnostic refusal rather than reporting a wrong deterministic value.
  for (tolerance in c(.02, .05)) {
    value <- .test_cert_native_cluster(blocks[[2L]], tolerance = tolerance,
                                      initial = 16L, points = 32L, scrambles = 2L)
    expect_false(isTRUE(as.logical(attr(value, "deterministic"))))
    expect_error(.selection_joint_factor_check_diagnostics(
      value, list(relative_tolerance = tolerance), "cluster"),
      "Increase 'max_points_per_scramble' or 'scrambles'", fixed = TRUE)
  }
})

test_that("the best weight rule keeps its geometric acceptance, decisions and values", {

  # Seeded random best-rule blocks (one- and two-sided; rank-one clusters and
  # rank-two and rank-three factors) with the outputs of the build before the
  # certificate: the best rule keeps the precheck, the exact ladder start and
  # the geometric acceptance, so decisions and values are unchanged. The values
  # are compared at 1e-12 relative so that the test also holds on other
  # compilers and instruction sets (the rounding of the long double sums
  # differs there). That the best path reproduces the build before the
  # certificate bitwise was checked as an A/B comparison of the two builds during
  # development, not in the suite. The randomized fallback of a block stops later
  # than that build did when its one-standard-error stage fails the scaled rule
  # (2.36 times the standard error); the stored output is then the stage where
  # the old rule stopped (.test_cert_expect_stored_fallback).
  reference <- readRDS(testthat::test_path("fixtures", "selection-certificate-best-reference.rds"))
  factor_rules  <- .selection_joint_factor_quadrature_rules()
  cluster_rules <- .selection_joint_cluster_quadrature_rules(SELNORM_CLUSTER_QUADRATURE_ORDERS)
  routes <- character()
  outcomes <- character()
  for (i in seq_along(reference[["blocks"]])) {
    b <- reference[["blocks"]][[i]]
    actual <- .test_cert_native_best(b, factor_rules, cluster_rules)
    expected <- reference[["results"]][[i]]
    info <- paste("block", i)
    expect_identical(actual[["relative_mcse"]] == 0, expected[["relative_mcse"]] == 0, info = info)
    routes <- c(routes, if (actual[["relative_mcse"]] == 0) "deterministic" else "fallback")
    fields <- c("log_density", "relative_mcse", "relative_change", "log_normalizer")
    if (actual[["relative_mcse"]] == 0) {
      for (field in fields) {
        expect_equal(as.numeric(actual[[field]]), expected[[field]], tolerance = 1e-12, info = paste(info, field))
      }
    } else {
      outcomes <- c(outcomes, .test_cert_expect_stored_fallback(actual, expected,
        function(points) .test_cert_native_best(b, factor_rules, cluster_rules, initial = points, points = points),
        stages = c(64L, 128L, 256L, 512L), fields = fields, info = info))
    }
  }
  expect_true(all(c("deterministic", "fallback") %in% routes))
  # the stopping multiplier moves some fallback blocks to a later stage, not all of them
  expect_true(all(c("unchanged", "later stage") %in% outcomes))
})
