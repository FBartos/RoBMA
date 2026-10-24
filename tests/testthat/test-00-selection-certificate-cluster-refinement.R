context("Rank-one cluster certificate refinement")
skip_on_cran()

# Independent four-shift averages of the scalar Gauss-Hermite rule. The gap is
# taken at each unshifted rung's own tilted mean; the native cluster probe does
# not supply these averages.
.test_cert_cluster_refinement <- function(block, orders = c(511L, 1023L)) {

  average <- function(order) {
    base <- .test_cert_tensor_log_integral(block, order, moments = TRUE)
    gap <- .test_cert_node_gap(order, base[["mean"]])
    values <- vapply(0:3, function(m) {
      .test_cert_tensor_log_integral(block, order, shift = m * gap / 4)
    }, numeric(1L))
    c(two = .test_cert_logsumexp(values[c(1L, 3L)]) - log(2),
      four = .test_cert_logsumexp(values) - log(4))
  }
  prior <- average(orders[[1L]])
  fine <- average(orders[[2L]])
  list(value = fine[["four"]],
    shift_change = abs(expm1(fine[["four"]] - fine[["two"]])),
    order_change = abs(expm1(fine[["four"]] - prior[["four"]])))
}

.test_cert_cluster_refinement_cases <- function() {

  list(
    # These two blocks previously exercised the cluster's QMC proposals. The
    # new tier resolves both after all seven original ladder rungs reject.
    rare_tail = .test_cert_block(list(.test_cert_parallel_group(4L,
      1.8707955473336113, 0.042233265807645956, -6.8303662251056725,
      0.22565989842880826)), c(1, .2)),
    wall = .test_cert_block(list(.test_cert_parallel_group(100L,
      0.65477525608481191, 0.026923363058843731, -2.5196777242881367,
      0.20064200965816639)), c(1, .05), c(1.645, -Inf)),
    band = .test_cert_block(list(.test_cert_parallel_group(128L,
      1, 0.0562285354392654, -0.617803297704086)), c(1, .05, 1), c(1, -1, -Inf)),
    zero_weight = .test_cert_block(list(.test_cert_parallel_group(2L,
      -1, 0.0290061025794998, -1.2016530171968)), c(1, 0))
  )
}

test_that("the full cluster ladder refines its final two rungs against independent integrals", {

  cases <- .test_cert_cluster_refinement_cases()
  for (name in names(cases)) {
    block <- cases[[name]]
    expect_null(.test_cert_emulate_ladder(block,
      orders = SELNORM_CLUSTER_QUADRATURE_ORDERS, changes = 1L), info = name)
    expect_true(.test_cert_resolves(block, 1023L, .test_cert_tier_width_scale(block)), info = name)
    refined <- .test_cert_cluster_refinement(block)
    change <- max(refined[["shift_change"]], refined[["order_change"]])
    expect_lte(change, .005 / 2, label = name)
    expect_lte(.test_cert_tail_bound(block, 1023L, refined[["value"]]), .005, label = name)

    actual <- .test_cert_native_cluster(block)
    expect_true(.test_cert_deterministic(actual), info = name)
    expect_identical(as.numeric(actual[["relative_mcse"]]), 0, info = name)
    expect_lt(abs(actual[["log_normalizer"]] - refined[["value"]]), 1e-11, label = name)
    expect_lt(abs(actual[["relative_change"]] - change), 1e-11, label = name)
    # Composite Gauss-Legendre panels are at most one eighth of the smallest
    # row transition width, independently of every Hermite rung and shift.
    reference <- .test_cert_reference_parallel(block)
    expect_lt(.test_cert_relative_error(actual[["log_normalizer"]], reference), .005, label = name)
    expect_identical(.test_cert_native_cluster(block, seed = 2L), actual, info = name)
  }
})

test_that("cluster refinement retains the calibrated tolerance and preceding-rung requirement", {

  block <- .test_cert_cluster_refinement_cases()[["rare_tail"]]
  expect_true(.test_cert_deterministic(.test_cert_native_cluster(block)))
  strict <- .test_cert_native_cluster(block, tolerance = .005 * (1 - 1e-8))
  expect_false(.test_cert_deterministic(strict))
  expect_gt(strict[["relative_mcse"]], 0)

  # A one-node preceding rung has no gap for its shifted average. This valid
  # two-rung schedule must retain QMC rather than access a missing node gap.
  no_previous_gap <- .test_cert_native_cluster(block,
    rules = .selection_joint_cluster_quadrature_rules(c(1L, 1023L)))
  # QMC scrambles can coincide and report zero MCSE; the explicit route flag
  # identifies fallback. This schedule retains its baseline QMC output.
  expect_false(.test_cert_deterministic(no_previous_gap))
  expect_error(.test_cert_native_cluster(block,
    rules = .selection_joint_cluster_quadrature_rules(1023L)),
    "Cluster selection controls are invalid", fixed = TRUE)
})

test_that("a cluster whose two-shift check fails keeps its randomized diagnostics", {

  block <- .test_cert_block(list(.test_cert_parallel_group(2L,
    1, 0.0157304352195527, -1.31269414862618)), c(1, 0))
  expect_null(.test_cert_emulate_ladder(block,
    orders = SELNORM_CLUSTER_QUADRATURE_ORDERS, changes = 1L))
  expect_true(.test_cert_resolves(block, 1023L, .test_cert_tier_width_scale(block)))
  refined <- .test_cert_cluster_refinement(block)
  expect_gt(refined[["shift_change"]], .005 / 2)
  expect_lt(refined[["order_change"]], .005 / 2)
  actual <- .test_cert_native_cluster(block, points = 4096L)
  expect_false(.test_cert_deterministic(actual))
  expect_gt(actual[["relative_mcse"]], 0)
  expect_lt(.test_cert_relative_error(actual[["log_normalizer"]],
    .test_cert_reference_parallel(block)), 3 * actual[["relative_mcse"]])
})
