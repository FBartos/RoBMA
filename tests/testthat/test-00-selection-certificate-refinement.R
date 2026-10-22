context("Selection certificate refinement tier")
skip_on_cran()

# The refinement tier of the factor kernel: after the certificate ladder rejects
# every affordable rung, the average of the finest rung over a lattice of shifts
# (m / M) h per axis, M = 2 and M = 4, with h the node gap at the rung's own
# tilted mean. It is accepted when the M = 4 average agrees with the M = 2
# average and with the M = 4 average of the previous affordable rung, each
# within half the tolerance, and the finest rung passes the tail gate. Tensors
# of rank at most two and forests of depth at most one only.

.test_cert_ref_group <- function(count, loading, sd, mean, sei = 1) {

  list(count = count, loading = loading, sd = sd, mean = mean, sei = sei)
}

.test_cert_ref_blocks <- function(omega = c(1, .2), sd_scale = .3) {

  list(
    # one latent axis
    rank_one = .test_cert_block(list(
      .test_cert_ref_group(4L, .8, sd_scale * 1.1, .3)), omega),
    # two independent roots, each with rows on its own axis (a forest of depth zero)
    two_roots = .test_cert_block(list(
      .test_cert_ref_group(3L, c(.8, 0), sd_scale * 1.1, .2),
      .test_cert_ref_group(2L, c(0, .7), sd_scale * 1.3, -.1)), omega),
    # a shared axis with two leaves, rows on the shared axis only and a free row
    two_leaves = .test_cert_block(list(
      .test_cert_ref_group(3L, c(.6, .5, 0), sd_scale * 1.1, .2),
      .test_cert_ref_group(2L, c(.6, 0, .4), sd_scale * 1.3, -.1),
      .test_cert_ref_group(2L, c(.7, 0, 0), sd_scale * 1.2, .1),
      .test_cert_ref_group(1L, c(0, 0, 0), sd_scale * 1, 0)), omega),
    # two independent trees, each a root with one leaf
    two_trees = .test_cert_block(list(
      .test_cert_ref_group(3L, c(.6, .5, 0, 0), sd_scale * 1.1, .2),
      .test_cert_ref_group(2L, c(0, 0, .5, .4), sd_scale * 1.3, -.1)), omega),
    # a group on both axes with rows on either axis alone: partial overlaps, no
    # forest, the rank-two tensor rule
    crossing = .test_cert_block(list(
      .test_cert_ref_group(3L, c(.7, .3), sd_scale * 1.1, .2),
      .test_cert_ref_group(2L, c(.5, 0), sd_scale * 1.2, -.1),
      .test_cert_ref_group(2L, c(0, .6), sd_scale * 1.3, .1)), omega)
  )
}

# Mean over the grid of shifts (m / grid) * gap, m = 0..grid-1 on every axis, of
# the plain-R tensor rule (log scale). `gap` is by loading column.
.test_cert_shift_average <- function(block, order, gap, grid) {

  lattice <- as.matrix(expand.grid(rep(list(seq_len(grid) - 1L), block[["rank"]])))
  value <- vapply(seq_len(nrow(lattice)), function(k) {
    .test_cert_tensor_log_integral(block, order, shift = lattice[k, ] / grid * gap)
  }, numeric(1L))
  .test_cert_logsumexp(value) - log(length(value))
}

.test_cert_ref_by_column <- function(probe, field, rank) {

  value <- rep(NA_real_, rank)
  value[probe[["column"]]] <- probe[[field]]
  value
}

test_that("the refinement averages equal the average over the full lattice of shifts", {

  # Given the shift of a root, the leaves of its tree are independent: the
  # product of the leaf integrals, each averaged over its own shifts, is the
  # average of the products over the whole lattice of root and leaf shifts.
  blocks <- .test_cert_ref_blocks()
  orders <- c(rank_one = 9L, two_roots = 9L, two_leaves = 11L, two_trees = 7L, crossing = 15L)
  for (name in names(blocks)) {
    block <- blocks[[name]]
    probe <- .test_cert_probe(block, orders[[name]])
    expect_true(probe[["refined"]], info = name)
    gap <- .test_cert_ref_by_column(probe, "gap", block[["rank"]])
    for (grid in c(2L, 4L)) {
      expected <- .test_cert_shift_average(block, orders[[name]], gap, grid)
      actual   <- probe[[paste0("refined_", grid)]]
      expect_lt(abs(expm1(actual - expected)), 1e-13, label = paste(name, "M =", grid))
    }
    # not vacuous: the coarse rung is far from its refined average
    expect_gt(abs(expm1(probe[["refined_4"]] - probe[["log_normalizer"]])), 1e-5, label = name)
  }
})

test_that("only forests of depth one and tensors of rank two are refined", {

  blocks <- .test_cert_ref_blocks()
  for (name in names(blocks)) {
    expect_true(.test_cert_probe(blocks[[name]], 7L)[["refined"]], info = name)
  }
  # nested supports of depth two, and a rank-three tensor (partial overlaps)
  chain <- .test_cert_block(list(
    .test_cert_ref_group(2L, c(.6, .4, .3), .33, .2),
    .test_cert_ref_group(2L, c(.5, .4, 0), .36, -.1),
    .test_cert_ref_group(2L, c(.7, 0, 0), .3, .1)), c(1, .2))
  tensor <- .test_cert_block(list(
    .test_cert_ref_group(2L, c(.6, .4, 0), .33, .2),
    .test_cert_ref_group(2L, c(0, .5, .4), .36, -.1),
    .test_cert_ref_group(2L, c(.5, 0, .6), .3, .1)), c(1, .2))
  expect_identical(.test_cert_probe(chain, 7L)[["path"]], .test_cert_paths[["nested"]])
  expect_false(.test_cert_probe(chain, 7L)[["refined"]])
  expect_identical(.test_cert_probe(tensor, 7L)[["path"]], .test_cert_paths[["tensor_direct"]])
  expect_false(.test_cert_probe(tensor, 7L)[["refined"]])
})

test_that("the kernel returns the refined average of the finest rung when the ladder rejects", {

  blocks <- list(
    fed80 = .test_cert_block(list(
      .test_cert_parallel_group(80L, c(.496441, .542314), .246872, -2.125758)), c(1, .01)),
    crossing = .test_cert_block(list(
      list(count = 92L, loading = c(1.59892823963330177, .46248906601712753),
           sd = 1.17935060431781880, mean = -.12799436126613484, sei = 1),
      list(count = 32L, loading = c(-.37056447402550319, .48747887335740198),
           sd = .23232799104534663, mean = -2.32439955742265258, sei = 1)),
      c(1, .03633489870625161)),
    forest371 = .test_cert_block(list(
      list(count = 115L, loading = c(.001, 1, 0), sd = .249373398318475, mean = 1.93906024750322, sei = 1),
      list(count = 256L, loading = c(.001, 0, -1), sd = .103488784220421, mean = 1.91516049392521, sei = 1)),
      c(1, .01)))
  references <- list(fed80 = .test_cert_reference_parallel, crossing = .test_cert_reference_two_groups,
    forest371 = .test_cert_reference_forest)
  for (name in names(blocks)) {
    block  <- blocks[[name]]
    actual <- .test_cert_native_factor(block)
    tier   <- .test_cert_tier(block)
    info   <- name
    # the certificate ladder rejects every rung, so the value is the tier's
    expect_true(is.na(.test_cert_probe_ladder(block)), info = info)
    if (block[["rank"]] == 2L) expect_null(.test_cert_emulate_ladder(block), info = info)
    expect_true(tier[["accepted"]], info = info)
    expect_true(.test_cert_deterministic(actual), info = info)
    expect_lt(abs(actual[["log_normalizer"]] - tier[["value"]]), 1e-12, label = info)
    expect_lt(abs(actual[["relative_change"]] - max(tier[["shift_change"]], tier[["order_change"]])),
              1e-12, label = info)
    expect_lt(max(tier[["shift_change"]], tier[["order_change"]]), .005 / 2, label = info)
    expect_lt(.test_cert_relative_error(actual[["log_normalizer"]], references[[name]](block)),
              1e-4, label = info)
  }
})

test_that("a block whose refinement checks fail goes on to the randomized fallback", {

  # Blocks of the research battery (seed 41) that the refinement tier does not
  # certify: the four-shift average still differs from the two-shift average by
  # more than half the tolerance (shift check), or from the average of the
  # previous rung (order check).
  parallel <- function(count, loading, sd, mean, sei, omega, z_lower) {
    .test_cert_block(list(.test_cert_parallel_group(count, loading, sd, mean, sei)), omega, z_lower)
  }
  blocks <- list(
    # shift check passes (1.5e-3), order check fails (1.9e-2)
    order_check = list(shift = FALSE, order = TRUE, block = parallel(17L,
      c(-0.76072098928229415, 0.016162300638867866), 0.031622080054956925,
      -0.23809113421998562, 1.169331828146867, c(1, 0.050000000000000003), c(1.96, -Inf))),
    # both fail (7.3e-3 and 3.7e-3)
    both_checks = list(shift = TRUE, order = TRUE, block = parallel(12L,
      c(0.014692097199723289, -0.8616553723389665), 0.13653164601012455,
      -5.6000270018551461, 0.34873594478420489, c(1, 0.050000000000000003), c(0, -Inf))),
    # nearly parallel rows: the shift check fails (7.8e-3)
    shift_check = list(shift = TRUE, order = FALSE, block = .test_cert_block(list(
      list(count = 1L, loading = c(-0.010868789953203364, 0.67241886559550057),
           sd = 0.037873483667476435, mean = 0.96681311864952701, sei = 1),
      list(count = 1L, loading = c(-0.012359303641840926, 0.67239312088782954),
           sd = 0.037873483667476435, mean = 0.96681311864952701, sei = 1),
      list(count = 1L, loading = c(-0.010512092253634437, 0.67242453652453649),
           sd = 0.037873483667476435, mean = 0.96681311864952701, sei = 1)),
      c(1, 0.29999999999999999), c(1.96, -Inf)))
  )
  for (name in names(blocks)) {
    case  <- blocks[[name]]
    tier  <- .test_cert_tier(case[["block"]])
    expect_true(is.na(.test_cert_probe_ladder(case[["block"]])), info = name)
    expect_false(tier[["accepted"]], info = name)
    expect_lte(tier[["tail"]], .005, label = name)
    expect_identical(tier[["shift_change"]] > .005 / 2, case[["shift"]], info = name)
    expect_identical(tier[["order_change"]] > .005 / 2, case[["order"]], info = name)
    expect_false(.test_cert_deterministic(.test_cert_native_factor(case[["block"]])), info = name)
  }
})
