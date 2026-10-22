context("Selection certificate shifted rules")
skip_on_cran()

# The shifted rungs of the shift-amplitude certificate, through the native
# development probe: the exact change of variables, the reuse of unshifted
# messages, the node-gap port, the tilted means and the decisions.

.test_cert_group <- function(count, loading, sd, mean, sei = 1) {

  list(count = count, loading = loading, sd = sd, mean = mean, sei = sei)
}

.test_cert_structures <- function(omega = c(1, .9), sd_scale = 1) {

  list(
    # shared axis with two leaves, rows on the shared axis only and a free row
    forest = .test_cert_block(list(
      .test_cert_group(3L, c(.6, .5, 0), sd_scale * 1.1, .2),
      .test_cert_group(2L, c(.6, 0, .4), sd_scale * 1.3, -.1),
      .test_cert_group(2L, c(.7, 0, 0), sd_scale * 1.2, .1),
      .test_cert_group(1L, c(0, 0, 0), sd_scale * 1, 0)), omega),
    # two independent two-axis groups
    multi_root = .test_cert_block(list(
      .test_cert_group(3L, c(.6, .5, 0, 0), sd_scale * 1.1, .2),
      .test_cert_group(2L, c(0, 0, .5, .4), sd_scale * 1.3, -.1)), omega),
    # nested supports of depth two
    chain = .test_cert_block(list(
      .test_cert_group(2L, c(.6, .4, .3), sd_scale * 1.1, .2),
      .test_cert_group(2L, c(.5, .4, 0), sd_scale * 1.2, -.1),
      .test_cert_group(2L, c(.7, 0, 0), sd_scale * 1.0, .1)), omega),
    # pairwise partial overlaps: not a forest, the tensor rule
    tensor = .test_cert_block(list(
      .test_cert_group(2L, c(.6, .4, 0), sd_scale * 1.1, .2),
      .test_cert_group(2L, c(0, .5, .4), sd_scale * 1.2, -.1),
      .test_cert_group(2L, c(.5, 0, .6), sd_scale * 1.0, .1)), omega)
  )
}

.test_cert_by_column <- function(probe, field, rank) {

  value <- rep(NA_real_, rank)
  value[probe[["column"]]] <- probe[[field]]
  value
}

test_that("shifted rules estimate the same integral as the rung (exact change of variables)", {

  shifts <- c(-.5, .3, 1)
  blocks <- .test_cert_structures()
  orders <- c(forest = 31L, multi_root = 31L, chain = 21L, tensor = 21L)
  for (name in names(blocks)) {
    block <- blocks[[name]]
    probe <- .test_cert_probe(block, orders[[name]],
      shift_axis = rep(seq_len(block[["rank"]]), each = length(shifts)),
      shift_delta = rep(shifts, block[["rank"]]))
    expect_true(is.finite(probe[["log_normalizer"]]), info = name)
    for (field in c("reuse_value", "full_value")) {
      expect_lt(max(abs(expm1(probe[[field]] - probe[["log_normalizer"]]))), 1e-12, label = paste(name, field))
    }
  }
})

test_that("a shifted rung of a coarse rule moves, and reusing messages equals recomputing them", {

  blocks <- .test_cert_structures(c(1, .2), sd_scale = .3)
  for (name in names(blocks)) {
    block <- blocks[[name]]
    probe <- .test_cert_probe(block, 9L,
      shift_axis = rep(seq_len(block[["rank"]]), each = 2L),
      shift_delta = rep(c(.4, -.25), block[["rank"]]))
    # not vacuous: the coarse rule is far from the shifted copies' common limit
    expect_gt(max(abs(expm1(probe[["full_value"]] - probe[["log_normalizer"]]))), 1e-3, label = name)
    expect_lt(max(abs(expm1(probe[["reuse_value"]] - probe[["full_value"]]))), 1e-14, label = name)
  }
})

test_that("the rung value of the certificate evaluators equals the unmodified evaluators", {

  blocks <- .test_cert_structures(c(1, .2))
  for (name in names(blocks)) for (order in c(7L, 15L, 21L)) {
    probe <- .test_cert_probe(blocks[[name]], order)
    expect_identical(probe[["log_normalizer"]], probe[["reference_log_normalizer"]],
                     info = paste(name, order))
  }
  expected_path <- c(forest = "nested", multi_root = "nested", chain = "nested", tensor = "tensor_direct")
  for (name in names(expected_path)) {
    expect_identical(.test_cert_probe(blocks[[name]], 9L)[["path"]],
                     .test_cert_paths[[expected_path[[name]]]], info = name)
  }
})

test_that("tilted means, shifted values and amplitudes follow the plain-R tensor certificate", {

  blocks <- .test_cert_structures(c(1, .2))
  for (name in names(blocks)) for (order in c(9L, 11L)) {
    block  <- blocks[[name]]
    probe  <- .test_cert_probe(block, order)
    oracle <- .test_cert_amplitudes(block, order)
    info   <- paste(name, order)
    expect_lt(abs(probe[["log_normalizer"]] - oracle[["log_integral"]]), 1e-10, label = info)
    for (field in c("mean", "gap", "half", "quarter", "amplitude")) {
      native <- .test_cert_by_column(probe, field, block[["rank"]])
      # rows on no axis keep a single axis out of the groups only when unused
      expect_lt(max(abs(native - oracle[[field]]), na.rm = TRUE), 1e-9, label = paste(info, field))
    }
    expect_equal(sort(probe[["column"]]), seq_len(block[["rank"]]), info = info)
    expect_lt(abs(probe[["amplitude_sum"]] - sum(oracle[["amplitude"]])), 1e-9, label = info)
  }
})

test_that("permuting the loading columns permutes the axes and keeps the rung", {

  block    <- .test_cert_structures(c(1, .2))[["forest"]]
  permute  <- c(3L, 1L, 2L)
  permuted <- block
  permuted[["loading"]] <- block[["loading"]][, permute]
  original <- .test_cert_probe(block, 11L)
  moved    <- .test_cert_probe(permuted, 11L)
  expect_lt(abs(original[["log_normalizer"]] - moved[["log_normalizer"]]), 1e-12)
  for (field in c("mean", "gap", "half", "quarter", "amplitude")) {
    # column j of the permuted block is column permute[j] of the original
    expect_lt(max(abs(.test_cert_by_column(moved, field, 3L) -
      .test_cert_by_column(original, field, 3L)[permute])), 1e-10, label = field)
  }
})

test_that("the node gap is the findInterval gap clamped to the interior intervals", {

  gap <- function(nodes, mean) .Call("RoBMA_selnorm_certificate_gap", nodes, mean, PACKAGE = "RoBMA")
  for (order in c(3L, 4L, 15L, 31L)) {
    nodes <- .gauss_hermite_nodes(order)[["nodes"]]
    means <- c(-Inf, nodes[1L] - 1, nodes[1L], (nodes[1L] + nodes[2L]) / 2, nodes,
               (nodes[-1L] + nodes[-order]) / 2, nodes[order] - 1e-12, nodes[order] + 1, Inf)
    expected <- vapply(means, function(m) .test_cert_node_gap(order, m), numeric(1L))
    expect_identical(gap(nodes, means), expected, info = order)
    # below the first node: the first gap; above the last: the last; on a node: to its right
    expect_identical(gap(nodes, nodes[1L] - 5), nodes[2L] - nodes[1L])
    expect_identical(gap(nodes, nodes[order] + 5), nodes[order] - nodes[order - 1L])
    expect_identical(gap(nodes, nodes[order]), nodes[order] - nodes[order - 1L])
    expect_identical(gap(nodes, nodes[2L]), nodes[3L] - nodes[2L])
  }
  # irregular nodes: the actual gap, not a local spacing formula
  nodes <- c(-3, -1, 0, 2, 2.5, 7)
  expect_identical(gap(nodes, c(-1, -.5, 1.9999, 2, 2.2, 6)), c(1, 1, 2, .5, .5, 4.5))
})

test_that("the half-shift early exit decides like the full certificate", {

  set.seed(61)
  early <- accepted <- 0L
  for (i in seq_len(36L)) {
    rows   <- sample(2:14, 1L)
    groups <- lapply(seq_len(rows), function(row) .test_cert_group(1L, stats::rnorm(2L, 0, .7),
      stats::runif(1L, .3, 1), stats::rnorm(1L, 0, 1.5)))
    block  <- .test_cert_block(groups, c(1, stats::runif(1L, 0, .8)))
    for (order in c(7L, 11L, 15L, 21L, 31L)) {
      probe <- .test_cert_probe(block, order, tolerance = c(.005, .02)[1L + i %% 2L])
      expect_identical(probe[["accepted"]], probe[["accepted_without_early_exit"]],
                       info = paste(i, order))
      early    <- early + probe[["early_reject"]]
      accepted <- accepted + probe[["accepted"]]
    }
  }
  # both outcomes occur
  expect_gt(early, 0L)
  expect_gt(accepted, 0L)
})

test_that("a block on the log-scale path is certified from that path", {

  block  <- .test_cert_structures(c(1, .2))[["forest"]]
  direct <- .test_cert_probe(block, 11L)
  logged <- .test_cert_probe(block, 11L, telescope = FALSE)
  expect_identical(direct[["path"]], .test_cert_paths[["nested"]])
  expect_identical(logged[["path"]], .test_cert_paths[["tensor_log"]])
  expect_identical(logged[["log_normalizer"]], logged[["reference_log_normalizer"]])
  oracle <- .test_cert_amplitudes(block, 11L)
  expect_lt(abs(logged[["log_normalizer"]] - oracle[["log_integral"]]), 1e-9)
  for (field in c("mean", "gap", "half", "quarter", "amplitude")) {
    expect_lt(max(abs(.test_cert_by_column(logged, field, 3L) - oracle[[field]])), 1e-8, label = field)
  }
  shifted <- .test_cert_probe(block, 11L, telescope = FALSE,
    shift_axis = 1:3, shift_delta = c(.3, -.2, .5))
  expect_true(all(is.finite(shifted[["full_value"]])))
  expect_identical(shifted[["full_value"]], shifted[["reuse_value"]])
})

test_that("the rank-one cluster rung gives its mean, shifts and amplitude without precomputed weights", {

  block <- .test_cert_block(list(.test_cert_group(6L, 1.3, .8, .5)), c(1, .3))
  for (order in c(15L, 31L)) {
    plain       <- .test_cert_probe(block, order, cluster = TRUE, shift_axis = c(1L, 1L), shift_delta = c(.3, -.4))
    precomputed <- .test_cert_probe(block, order, cluster = TRUE, shift_axis = c(1L, 1L),
      shift_delta = c(.3, -.4), precomputed = TRUE)
    oracle <- .test_cert_amplitudes(block, order)
    expect_identical(plain[["path"]], .test_cert_paths[["cluster"]])
    expect_identical(plain[["log_normalizer"]], plain[["reference_log_normalizer"]])
    # the shifted copies build their own weights and ignore a precomputed rule
    expect_identical(precomputed[["log_normalizer"]], plain[["log_normalizer"]])
    expect_identical(precomputed[["full_value"]], plain[["full_value"]])
    expect_identical(precomputed[["amplitude"]], plain[["amplitude"]])
    expect_lt(abs(plain[["log_normalizer"]] - oracle[["log_integral"]]), 1e-10)
    for (field in c("mean", "gap", "half", "quarter", "amplitude")) {
      expect_lt(max(abs(plain[[field]] - oracle[[field]])), 1e-9, label = paste(order, field))
    }
  }

  smooth <- .test_cert_block(list(.test_cert_group(3L, .6, 1.2, .2)), c(1, .9))
  probe  <- .test_cert_probe(smooth, 127L, cluster = TRUE, precomputed = TRUE,
    shift_axis = c(1L, 1L, 1L), shift_delta = c(-.5, .3, 1))
  expect_lt(max(abs(expm1(probe[["full_value"]] - probe[["log_normalizer"]]))), 1e-12)
})

test_that("blocks without a transition carry no amplitude to certify", {

  flat  <- .test_cert_block(list(.test_cert_group(4L, c(.6, .5), 1.1, .2)), c(1, 1))
  probe <- .test_cert_probe(flat, 21L)
  expect_false(probe[["transitions"]])
  step  <- .test_cert_block(list(.test_cert_group(4L, c(.6, .5), 1.1, .2)), c(1, .5))
  expect_true(.test_cert_probe(step, 21L)[["transitions"]])
})
