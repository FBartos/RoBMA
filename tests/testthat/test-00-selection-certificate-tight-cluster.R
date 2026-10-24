context("Selection cluster refinement at tight tolerances")
skip_on_cran()

test_that("cluster refinement agreement does not authorize tighter tolerances", {

  # A wall found by varying its location and width against the cluster ladder.
  # Both refinement comparisons are below 5e-10, but their common value has
  # relative error 2.12e-6. The count-adjusted sharpness gate also passes.
  # Composite Gauss-Legendre and adaptive integration independently gave
  # -3.2360574428815934 and -3.2360574428815898 for its log normalizer.
  counts <- c(217L, 361L, 371L, 30L, 561L, 174L)
  means <- c(-1.53267527128502, -1.50209298115994, -1.53252703426667,
             -1.48899978923644, -1.51248154381973, -1.44679752866181)
  groups <- Map(function(count, mean) {
    .test_cert_parallel_group(count, 1, .0715777471688762, mean)
  }, counts, means)
  block <- .test_cert_block(groups, c(1, .122769973506572))
  reference <- .test_cert_reference_parallel(block)
  expect_equal(reference, -3.2360574428815934, tolerance = 1e-12)

  # Group the coincident rows to keep the independent R evaluation small.
  evaluate <- function(order, shift = 0) {

    rule <- .gauss_hermite_nodes(order)
    nodes <- rule[["nodes"]]
    values <- rule[["log_weights"]] - shift * nodes - shift^2 / 2
    for (group in groups) {
      values <- values + group[["count"]] * .test_cert_log_probability(
        group[["mean"]] + nodes + shift, group[["sd"]], group[["sei"]],
        block[["omega"]], block[["z_lower"]])
    }
    .test_cert_logsumexp(values)
  }
  average <- function(order) {

    probe <- .test_cert_probe(block, order, cluster = TRUE, refine = FALSE)
    values <- vapply(0:3, function(phase) {
      evaluate(order, phase * probe[["gap"]][1L] / 4)
    }, numeric(1L))
    list(
      m2 = .test_cert_logsumexp(values[c(1L, 3L)]) - log(2),
      m4 = .test_cert_logsumexp(values) - log(4),
      gate = probe[["resolves_tier"]]
    )
  }
  fine <- average(1023L)
  previous <- average(511L)
  shift_change <- abs(expm1(fine[["m4"]] - fine[["m2"]]))
  order_change <- abs(expm1(fine[["m4"]] - previous[["m4"]]))
  tolerance <- 1e-6
  expect_true(fine[["gate"]])
  expect_lt(shift_change, 1e-8)
  expect_lt(order_change, 1e-8)
  expect_lt(.test_cert_tail_bound(block, 1023L, fine[["m4"]]), tolerance)
  expect_gt(.test_cert_relative_error(fine[["m4"]], reference), 2 * tolerance)

  # Every rung fails its successive-change check even at the default tolerance.
  values <- vapply(SELNORM_CLUSTER_QUADRATURE_ORDERS, evaluate, numeric(1L))
  expect_true(all(abs(expm1(-diff(values))) > .005))
  default <- .test_cert_native_cluster(block)
  expect_true(.test_cert_deterministic(default))
  expect_equal(default[["log_normalizer"]], fine[["m4"]], tolerance = 1e-11)
  expect_lt(.test_cert_relative_error(default[["log_normalizer"]], reference), .005)

  # At the tighter tolerance the calibrated tier is bypassed. A small fallback
  # budget may be inaccurate, but it must retain its randomized diagnostics.
  actual <- .test_cert_native_cluster(block, tolerance = tolerance,
    initial = 64L, points = 512L, scrambles = 4L)
  expect_false(.test_cert_deterministic(actual))
  expect_gt(actual[["relative_mcse"]], 0)
  expect_true(.test_cert_relative_error(actual[["log_normalizer"]], reference) <= tolerance ||
    max(.selection_qmc_z() * actual[["relative_mcse"]], actual[["relative_change"]]) > tolerance)
})
