context("Selection certificate refinement tier")
skip_on_cran()

# The refinement tier of the factor kernel: after the certificate ladder rejects
# every affordable rung, the average of the finest rung over a lattice of shifts
# (m / M) h per axis, M = 2 and M = 4, with h the node gap at the rung's own
# tilted mean. It is accepted when the M = 4 average agrees with the M = 2
# average and with the M = 4 average of the previous affordable rung, each
# within half the tolerance, and the finest rung passes the tail gate. Tensors
# of rank at most two and forests of depth at most one only. The rung must also
# pass the sharpness gate of the tier, whose step widths narrow with the number
# of rows that share an edge (a wall of coincident rows hides errors of the
# refined lattice from both checks, see the wall tests below).

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
  blocks <- .test_cert_small_blocks()
  orders <- c(rank_one = 9L, two_roots = 9L, two_leaves = 11L, two_trees = 7L, crossing = 15L)
  for (name in names(blocks)) {
    block <- blocks[[name]]
    probe <- .test_cert_probe(block, orders[[name]])
    expect_true(probe[["refined"]], info = name)
    gap <- .test_cert_ref_by_column(probe, "gap", block[["rank"]])
    for (grid in c(2L, 4L)) {
      expected <- .test_cert_shift_average(block, orders[[name]], gap, grid)
      actual   <- probe[[paste0("refined_", grid)]]
      expect_lt(abs(expm1(actual - expected)), 1e-14, label = paste(name, "M =", grid))
    }
    # not vacuous: the coarse rung is far from its refined average
    expect_gt(abs(expm1(probe[["refined_4"]] - probe[["log_normalizer"]])), 1e-5, label = name)
  }
})

test_that("only forests of depth one and tensors of rank two are refined", {

  blocks <- .test_cert_small_blocks()
  for (name in names(blocks)) {
    expect_true(.test_cert_probe(blocks[[name]], 7L)[["refined"]], info = name)
  }
  # nested supports of depth two, and a rank-three tensor (partial overlaps)
  chain <- .test_cert_block(list(
    .test_cert_small_group(2L, c(.6, .4, .3), .33, .2),
    .test_cert_small_group(2L, c(.5, .4, 0), .36, -.1),
    .test_cert_small_group(2L, c(.7, 0, 0), .3, .1)), c(1, .2))
  tensor <- .test_cert_block(list(
    .test_cert_small_group(2L, c(.6, .4, 0), .33, .2),
    .test_cert_small_group(2L, c(0, .5, .4), .36, -.1),
    .test_cert_small_group(2L, c(.5, 0, .6), .3, .1)), c(1, .2))
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
      c(1, .03633489870625161)))
  references <- list(fed80 = .test_cert_reference_parallel, crossing = .test_cert_reference_two_groups)
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

test_that("a forest whose rows narrow the tier gate beyond its edges reaches the fallback", {

  # The 371-row forest of the research (a shared axis with loading .001 and two
  # leaf axes of opposite orientation, 115 and 256 coincident rows). Its tier
  # checks pass and its refined value is within 1e-4 of the reference, but the
  # gate asks the 371 rows on the shared axis for edges of 0.43 of the node gap
  # (the coupled rows are compared with the widest central gap, 0.33, against
  # 2.33 x the leaf SD 0.10) and refuses the rung: the block is no longer
  # certified by the tier. It reaches the randomized fallback, which is within
  # the tolerance of the reference.
  block <- .test_cert_block(list(
    list(count = 115L, loading = c(.001, 1, 0), sd = .249373398318475, mean = 1.93906024750322, sei = 1),
    list(count = 256L, loading = c(.001, 0, -1), sd = .103488784220421, mean = 1.91516049392521, sei = 1)),
    c(1, .01))
  reference <- .test_cert_reference_forest(block)
  tier      <- .test_cert_tier(block)
  expect_true(is.na(.test_cert_probe_ladder(block)))
  expect_true(tier[["checks"]])
  expect_lt(.test_cert_relative_error(tier[["value"]], reference), 1e-4)
  expect_true(.test_cert_resolves(block, 95L, .test_cert_width_scale_tier))
  expect_false(.test_cert_resolves(block, 95L, .test_cert_tier_width_scale(block)))
  expect_false(tier[["accepted"]])
  actual <- .test_cert_native_factor(block)
  expect_false(.test_cert_deterministic(actual))
  expect_gt(actual[["relative_mcse"]], 0)
  expect_lt(.test_cert_relative_error(actual[["log_normalizer"]], reference), .005)
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


# Walls of many coincident rows: the product of N steps at (nearly) one place is
# narrower than one step, and the harmonics of the refined lattice that the
# tier's two checks cannot see then leave a four-shift average wrong by more than
# the tolerance although both checks pass. Each block is a wall on the first
# latent axis, groups of coincident rows with edges within one node gap of the
# order-95 rung, plus one benign row on the second axis (a genuine rank-two
# block whose integral separates into two one-dimensional integrals). The native
# tier returned a deterministic value 13-44% over the tolerance for all of them
# before the gate narrowed with the rows of a wall.
.test_wall_block <- function(t0, f, phase, counts, w2, cut = 0) {

  gap    <- .test_cert_node_gap(95L, t0)
  groups <- Map(function(count, shift) {
    .test_cert_parallel_group(count, c(1, 0), f * gap, cut - (t0 + shift * gap))
  }, counts, phase)
  groups[[length(groups) + 1L]] <- .test_cert_parallel_group(1L, c(0, .5), 1, 0)
  .test_cert_block(groups, c(1, w2), c(cut, -Inf))
}

.test_wall_blocks <- function() {

  list(
    # 1,885 rows, weight ratio 0.16; the native tier value was 5.66e-3 off
    rows_1885 = .test_wall_block(1.6579623735980198, 0.25013219923606989,
      c(0, 0.11463717561324255, 0.2140615014357225, 0.34593264733960727,
        0.39613751196415065, 0.84293796738559346, 0.86655454325784653),
      c(12L, 414L, 461L, 355L, 1L, 232L, 409L), 0.16408394810806834),
    # 382 rows, weight ratio 0.02 (5.67e-3)
    rows_382 = .test_wall_block(1.703107919746371, 0.19253477477881378,
      c(0.31637816326468388, 0.64740204987636607, 0.66698628733865917,
        0.69550095548126811, 0.70583813432863307, 0.81818244639781434,
        0.84752143871610663, 1),
      c(8L, 8L, 85L, 3L, 51L, 62L, 78L, 86L), 0.02),
    # 592 rows in ten groups, weight ratio 0.02 (6.17e-3)
    rows_592 = .test_wall_block(1.6880720027710769, 0.19411067606025131,
      c(0.29349067457951605, 0.32465500372921557, 0.38307427920195997,
        0.64642551759727152, 0.74421513990069443, 0.76199385894799754,
        0.89207790611240567, 0.90579905921208725, 0.91079307921694319, 1),
      c(19L, 6L, 22L, 2L, 90L, 28L, 75L, 120L, 120L, 109L), 0.02),
    # 2,285 rows behind the cut 1.96 (5.75e-3)
    rows_2285 = .test_wall_block(1.676784259488975, 0.20867182559218952,
      c(0.15751853348933853, 0.26851360715113781, 0.3291741060971104,
        0.54908596548395849, 0.80868242339161378, 0.90206173135818646,
        0.95240261906150347),
      c(219L, 20L, 1000L, 31L, 166L, 573L, 275L), 0.38568548174855827, cut = 1.96),
    # 521 rows with the wall far in the tail, at z = 3.97 (7.22e-3)
    tail_wall = .test_wall_block(3.9709134940021951, 0.28650959105080109,
      c(0.2889351625976383, 0.49619592853865413, 0.52486973812968452,
        0.72541792092979407, 0.81832225466228026, 0.85065277806872586,
        0.97577916999207648),
      c(8L, 243L, 13L, 131L, 93L, 26L, 6L), 0.5072226005824273, cut = 1.96)
  )
}

test_that("walls of coincident rows are not certified by the tier", {

  tolerance <- .005
  blocks    <- .test_wall_blocks()
  for (name in names(blocks)) {
    block     <- blocks[[name]]
    reference <- .test_cert_reference_parallel(block)
    tier      <- .test_cert_tier(block)
    # the tier's checks pass for the finest rung and its refined value is wrong
    expect_true(tier[["checks"]], info = name)
    expect_gt(.test_cert_relative_error(tier[["value"]], reference), tolerance, label = name)
    # the gate with the plain tier scale resolves the edges; the scale narrowed
    # by the rows of the wall does not, and it is the gate that refuses the rung
    expect_true(.test_cert_resolves(block, 95L, .test_cert_width_scale_tier), info = name)
    expect_false(.test_cert_resolves(block, 95L, .test_cert_tier_width_scale(block)), info = name)
    expect_false(tier[["resolves"]], info = name)
    expect_false(tier[["accepted"]], info = name)
    expect_true(is.na(.test_cert_probe_ladder(block)), info = name)
    # native: no deterministic value over the tolerance; the block reaches the
    # randomized fallback with its Monte Carlo error
    actual <- .test_cert_native_factor(block)
    if (.test_cert_deterministic(actual)) {
      expect_lt(.test_cert_relative_error(actual[["log_normalizer"]], reference),
                tolerance, label = paste(name, "relative error"))
    }
    expect_false(.test_cert_deterministic(actual), label = paste(name, "deterministic"))
    expect_gt(actual[["relative_mcse"]], 0, label = paste(name, "relative MCSE"))
  }
})
