context("Selection certificate sharpness gate")
skip_on_cran()

# The shift-amplitude certificate reads the aliasing error of a rung from its
# first harmonic. An edge of the selection weight much narrower than the node gap
# makes that error sawtooth-like, with harmonics that the three sampled phases
# cannot see, so a rung is certified only when the transition rule resolves its
# transitions with every step width relaxed by 1 / 0.15 (edge width at least
# 0.15 of the node gap); the refinement tier asks for 0.125 at its rung, more
# when many rows share an edge (see the tier tests). These tests pin the width
# scales of the native rule against its R mirror, its exact boundary, the
# decisions of the ladders built on it, and the reason for the gate.

.test_sharp_tolerance <- .005

.test_sharp_probe <- function(block, order, cluster = FALSE) {

  .test_cert_probe(block, order, cluster = cluster, refine = FALSE)
}

# Rows of a random block of coincident groups on `rank` axes. The widths cover
# edges from far below to far above the node gaps of the rungs.
.test_sharp_random_block <- function(rank, groups) {

  two_cuts <- stats::runif(1L) < .5
  omega    <- if (two_cuts) c(1, .6, stats::runif(1L, .05, .5)) else c(1, stats::runif(1L, .1, .9))
  cuts     <- if (two_cuts) c(1.96, 0, -Inf) else c(0, -Inf)
  .test_cert_block(lapply(seq_len(groups), function(g) {
    coupled <- rank > 1L && stats::runif(1L) < .5
    loading <- if (coupled) stats::rnorm(rank, 0, .8) else {
      value <- numeric(rank); value[sample(rank, 1L)] <- stats::rnorm(1L, 0, .8); value
    }
    list(count = sample(1:6, 1L), loading = loading, sd = exp(stats::runif(1L, log(.01), log(1))),
         mean = stats::rnorm(1L, 0, 2), sei = stats::runif(1L, .6, 1.4))
  }), omega, cuts)
}

test_that("the native gate equals the R rule at the scaled step widths", {

  set.seed(20261022)
  decisions <- list(ladder = logical(), tier = logical())
  for (i in seq_len(24L)) {
    rank  <- sample(1:3, 1L)
    block <- .test_sharp_random_block(rank, sample(1:3, 1L))
    # a one-axis block is evaluated by the cluster kernel's rule orders
    cluster <- rank == 1L
    orders  <- if (cluster) c(15L, 31L, 63L, 127L) else c(5L, 9L, 15L, 31L)
    for (order in orders) {
      probe <- .test_sharp_probe(block, order, cluster)
      info  <- paste("block", i, "order", order)
      expect_identical(probe[["resolves_ladder"]],
                       .test_cert_resolves(block, order, .test_cert_width_scale_ladder), info = info)
      expect_equal(probe[["tier_width_scale"]], .test_cert_tier_width_scale(block),
                   tolerance = 1e-14, info = info)
      expect_identical(probe[["resolves_tier"]],
                       .test_cert_resolves(block, order, .test_cert_tier_width_scale(block)),
                       info = info)
      decisions[["ladder"]] <- c(decisions[["ladder"]], probe[["resolves_ladder"]])
      decisions[["tier"]]   <- c(decisions[["tier"]], probe[["resolves_tier"]])
    }
  }
  # the geometries separate sharp from resolved rungs for both scales
  expect_true(any(decisions[["ladder"]]) && any(!decisions[["ladder"]]))
  expect_true(any(decisions[["tier"]]) && any(!decisions[["tier"]]))
})

test_that("the gate flips exactly where the edge width reaches its fraction of the node gap", {

  # One rank-one group: its transition sits in the middle of an interior node
  # interval, so the rule compares that interval's gap with scale x residual SD
  # (uncoupled axis). The fractions 0.15 (ladder) and 0.125 / sqrt(max(1, 2 ln N))
  # (tier, N rows sharing the edge) of the gap are the thresholds of the SD, for
  # any loading.
  order  <- 31L
  nodes  <- .gauss_hermite_nodes(order)[["nodes"]]
  k      <- 20L
  gap    <- nodes[k + 1L] - nodes[k]
  center <- (nodes[k] + nodes[k + 1L]) / 2
  wall   <- function(count, loading, sd) .test_cert_block(list(list(
    count = count, loading = loading, sd = sd, mean = -loading * center, sei = 1)), c(1, .5))
  for (loading in c(1, 2.5, -.7)) {
    for (count in c(1L, 3L, 40L)) {
      # the tier scale narrows with the rows that share the edge
      scales <- c(ladder = .test_cert_width_scale_ladder,
                  tier = .test_cert_tier_width_scale(wall(count, loading, 1)))
      for (gate in names(scales)) {
        for (side in c(-1, 1)) {
          # the SD at which the edge width sd / |loading| is the fraction 1 / scale of the gap
          sd    <- gap * abs(loading) / scales[[gate]] * (1 + side * 1e-6)
          block <- wall(count, loading, sd)
          probe <- .test_sharp_probe(block, order, cluster = TRUE)
          info  <- paste("loading", loading, "rows", count, gate, "side", side)
          expect_identical(probe[[paste0("resolves_", gate)]], side > 0, info = info)
          expect_identical(.test_cert_resolves(block, order, scales[[gate]]), side > 0, info = info)
        }
      }
    }
  }

  # A row that loads on two axes couples them: the rule then compares the largest
  # node gap over |z| <= 3 with the widths of the rows of both axes.
  nodes_nine <- .gauss_hermite_nodes(9L)[["nodes"]]
  coupled_gap <- max(diff(nodes_nine)[utils::head(nodes_nine, -1L) <= 3 & utils::tail(nodes_nine, -1L) >= -3])
  loading <- c(.6, .8)
  for (side in c(-1, 1)) {
    sd <- coupled_gap * max(abs(loading)) / .test_cert_width_scale_ladder * (1 + side * 1e-6)
    block <- .test_cert_block(list(list(count = 3L, loading = loading, sd = sd, mean = .2, sei = 1)),
                              c(1, .5))
    expect_identical(.test_sharp_probe(block, 9L)[["resolves_ladder"]], side > 0,
                     info = as.character(side))
  }
})

test_that("the tier gate asks more of a wall of many rows than of a single row", {

  # Edge width as a fraction of the gap of one interior node interval of the
  # order-31 rule. The ladder needs 0.15 of the gap whatever the number of rows.
  # The tier needs 0.125 for a single row, where it admits what the ladder
  # refuses, and 0.125 sqrt(2 ln N) for N rows sharing the edge (0.27 for 10 rows,
  # 0.46 for 1000): the product of N coincident steps is a wall narrower than one
  # step, and it refuses edges that the ladder admits.
  order  <- 31L
  nodes  <- .gauss_hermite_nodes(order)[["nodes"]]
  k      <- 20L
  gap    <- nodes[k + 1L] - nodes[k]
  center <- (nodes[k] + nodes[k + 1L]) / 2
  cases <- data.frame(
    rows     = c(1L, 1L, 1L, 10L, 10L, 1000L, 1000L),
    fraction = c(.10, .14, .20, .20, .30, .30, .50),
    ladder   = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
    tier     = c(FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE)
  )
  for (i in seq_len(nrow(cases))) {
    block <- .test_cert_block(list(list(count = cases[["rows"]][i], loading = 1,
      sd = cases[["fraction"]][i] * gap, mean = -center, sei = 1)), c(1, .5))
    probe <- .test_sharp_probe(block, order, cluster = TRUE)
    info  <- paste(cases[["rows"]][i], "rows, edge", cases[["fraction"]][i], "of the gap")
    expect_equal(probe[["tier_width_scale"]],
                 8 / sqrt(max(1, 2 * log(cases[["rows"]][i]))), tolerance = 1e-14, info = info)
    expect_identical(probe[["resolves_ladder"]], cases[["ladder"]][i], info = info)
    expect_identical(probe[["resolves_tier"]], cases[["tier"]][i], info = info)
  }
})

test_that("a block without transitions is always resolved", {

  block <- .test_cert_block(list(list(count = 3L, loading = c(.6, .5), sd = .001, mean = .2, sei = 1)),
                            c(1, 1))
  probe <- .test_sharp_probe(block, 5L)
  expect_false(probe[["transitions"]])
  expect_true(probe[["resolves_ladder"]])
  expect_true(probe[["resolves_tier"]])
})

# The nine sharp blocks, in order: six staircase blocks (four equal-height steps
# a quarter node gap apart, edge width 0.018 to 0.102 of the gap) and three
# failing forests of an adversarial search (a sharp leaf and a wide one, edge
# width 0.042 of the gap at rung 15). The certificate alone accepted each with a
# relative error of 5.0e-3 to 2.04e-2.
.test_sharp_staircase <- function(width, means, factor) {

  groups <- Map(function(count, mean) {
    list(count = count, loading = if (factor) c(1, 0) else 1, sd = width, mean = mean, sei = 1)
  }, c(23L, 13L, 9L, 7L), means)
  if (factor) groups[[5L]] <- list(count = 1L, loading = c(0, .5), sd = 1, mean = 0, sei = 1)
  .test_cert_block(groups, c(1, .97))
}

.test_sharp_forest <- function(leaves, omega) {

  .test_cert_block(lapply(seq_along(leaves), function(l) {
    leaf    <- leaves[[l]]
    loading <- numeric(1L + length(leaves))
    loading[c(1L, 1L + l)] <- c(leaf[[2L]], leaf[[1L]])
    list(count = leaf[[5L]], loading = loading, sd = leaf[[4L]], mean = leaf[[3L]], sei = 1)
  }), omega)
}

.test_sharp_cases <- function() {

  steps <- list(
    f01 = c(0.14999999999999999, 0.010013103039749255, -0.12997379392050148, -0.26996069088075225),
    f02 = c(0.14999999999999999, 0.030896278472876404, -0.088207443054247187, -0.2073111645813708),
    f035 = c(0.17499999999999999, 0.035013103039749249, -0.10497379392050149, -0.24496069088075223),
    c125 = c(0.125, 0.026430946114204684, -0.072138107771590632, -0.17070716165738598),
    c200 = c(-0.20000000000000001, -0.26955752269625372, -0.33911504539250736, -0.40867256808876107))
  list(
    # factor blocks: the staircase on one axis and a benign row on a second
    stair_factor_w01 = list(kind = "factor", order = 31L, pinned = -0.55120497817541914,
      block = .test_sharp_staircase(.01, steps[["f01"]], TRUE)),
    stair_factor_w02 = list(kind = "factor", order = 43L, pinned = -0.53557274901208896,
      block = .test_sharp_staircase(.02, steps[["f02"]], TRUE)),
    stair_factor_w035 = list(kind = "factor", order = 31L, pinned = -0.53929641127242878,
      block = .test_sharp_staircase(.035, steps[["f035"]], TRUE)),
    # rank-one clusters; the last has its steps at the spacing of order 127
    stair_cluster_w02 = list(kind = "cluster", order = 63L, pinned = -0.51781202989719155,
      block = .test_sharp_staircase(.02, steps[["c125"]], FALSE)),
    stair_cluster_w01 = list(kind = "cluster", order = 63L, pinned = -0.51726157140565476,
      block = .test_sharp_staircase(.01, steps[["c125"]], FALSE)),
    stair_cluster_w01_127 = list(kind = "cluster", order = 127L, pinned = -0.67676298891803555,
      block = .test_sharp_staircase(.01, steps[["c200"]], FALSE)),
    # depth-one forests: a root with two leaves (failing forests of the search)
    forest_w3 = list(kind = "forest", order = 15L, pinned = -3.7567291859187097,
      block = .test_sharp_forest(list(
        list(1.3974906387562571, 0.001, 3.5847573875242569, 0.050000000000000003, 154L),
        list(-0.30073570385413184, -0.001, 0.29715387627494971, 0.5, 48L)), c(1, 0.3))),
    forest_w5_a = list(kind = "forest", order = 15L, pinned = -1.8225403247458907,
      block = .test_sharp_forest(list(
        list(1.4833930865976188, 0.0013631474275185578, 3.7945195155167086, 0.050000000000000003, 252L),
        list(-0.31578158945293461, -0.0013631474275185578, 0.33126352892326133,
             0.43900707378758846, 26L)), c(1, 0.5))),
    forest_w5_b = list(kind = "forest", order = 15L, pinned = -1.6029106549630621,
      block = .test_sharp_forest(list(
        list(1.4286127189941824, 0.0017909370485516987, 3.6693119508996848,
             0.051778904646699535, 142L),
        list(-0.40507943508703881, -0.0017909370485516987, 0.39634593039578653,
             0.47969819025573263, 27L)), c(1, 0.5)))
  )
}

test_that("sharp-edged blocks that the certificate alone accepted wrongly are not certified", {

  # Native outcome: every block reaches the randomized fallback. The edges are
  # narrower than 0.15 of the gap at every rung of their ladder (the finest rung
  # of the cluster and factor ladders still has gaps above 0.1 and 0.3), and the
  # refinement tier needs at least 0.125 of the gap (more for many rows) at the
  # factor rung 95. A deterministic value, if one were returned, would have to
  # be right.
  cases <- .test_sharp_cases()
  for (name in names(cases)) {
    case      <- cases[[name]]
    block     <- case[["block"]]
    reference <- switch(case[["kind"]],
      forest = .test_cert_reference_forest(block), .test_cert_reference_parallel(block))
    # a changed reference must be visible
    expect_lt(abs(reference - case[["pinned"]]), 1e-8, label = paste(name, "reference"))
    actual <- if (case[["kind"]] == "cluster") .test_cert_native_cluster(block) else
      .test_cert_native_factor(block)
    if (.test_cert_deterministic(actual)) {
      expect_lt(.test_cert_relative_error(actual[["log_normalizer"]], reference),
                .test_sharp_tolerance, label = paste(name, "relative error"))
    }
    expect_false(.test_cert_deterministic(actual), label = paste(name, "deterministic"))
  }
})

test_that("the gate is what rejects the rungs that the certificate alone accepted", {

  # At the rung that the certificate alone accepted (the order recorded per
  # block) the amplitude and the combined allowance, with the factor 1.3 that
  # the ungated rule used, pass and the rung is wrong by more than half the
  # tolerance; only the sharpness gate refuses it (the allowance factor 2 refuses
  # the first forest as well). Its edges are between 0.018 and 0.102 of the gap
  # there.
  cases <- .test_sharp_cases()
  for (name in names(cases)) {
    case      <- cases[[name]]
    block     <- case[["block"]]
    cluster   <- case[["kind"]] == "cluster"
    probe     <- .test_sharp_probe(block, case[["order"]], cluster)
    error     <- .test_cert_relative_error(probe[["log_normalizer"]], case[["pinned"]])
    ungated <- probe[["finite"]] && !probe[["early_reject"]] &&
      probe[["amplitude_sum"]] <= .test_sharp_tolerance / 2 &&
      probe[["tail_bound"]] + 1.3 * probe[["amplitude_sum"]] <=
        .test_sharp_tolerance / (1 + .test_sharp_tolerance)
    expect_true(ungated, label = paste(name, "ungated certificate"))
    expect_gt(error, .test_sharp_tolerance / 2, label = paste(name, "error of the rung"))
    expect_false(probe[["resolves_ladder"]], label = paste(name, "ladder gate"))
    expect_false(probe[["resolves_tier"]], label = paste(name, "tier gate"))
  }
})

test_that("the native ladder follows the R emulation with the gate on sharp random blocks", {

  # Two-group blocks of narrow edges (SD 0.03 to 0.2): rungs that the certificate
  # accepts are refused for sharpness until the gaps are fine enough, so the
  # accepted order is set by the gate. Blocks that the ladder rejects reach the
  # refinement tier or the fallback, as the tier checks of the probe decide.
  set.seed(77)
  routes <- character()
  for (i in seq_len(10L)) {
    group <- function() list(count = sample(2:7, 1L), loading = stats::rnorm(2L, 0, .9),
      sd = stats::runif(1L, .03, .2), mean = stats::rnorm(1L, 0, 1.2), sei = 1)
    block  <- .test_cert_block(list(group(), group()), c(1, stats::runif(1L, .2, .8)))
    actual <- .test_cert_native_factor(block)
    expected <- .test_cert_emulate_ladder(block)
    info <- paste("block", i)
    if (!is.null(expected)) {
      routes <- c(routes, "ladder")
      expect_identical(actual[["relative_mcse"]], 0, info = info)
      expect_equal(actual[["log_normalizer"]], expected[["value"]], tolerance = 1e-10, info = info)
      expect_equal(actual[["relative_change"]], expected[["change"]], tolerance = 1e-8, info = info)
      expect_lt(.test_cert_relative_error(actual[["log_normalizer"]],
        .test_cert_reference_two_groups(block)), .test_sharp_tolerance, label = info)
    } else {
      tier <- .test_cert_tier(block)
      routes <- c(routes, if (tier[["accepted"]]) "tier" else "fallback")
      expect_identical(.test_cert_deterministic(actual), tier[["accepted"]], info = info)
      if (tier[["accepted"]]) {
        expect_equal(actual[["log_normalizer"]], tier[["value"]], tolerance = 1e-10, info = info)
        expect_lt(.test_cert_relative_error(actual[["log_normalizer"]],
          .test_cert_reference_two_groups(block)), .test_sharp_tolerance, label = info)
      }
    }
  }
  expect_true(any(routes == "ladder"))
})
