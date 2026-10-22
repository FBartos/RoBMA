context("Selection double-precision node products")
skip_on_cran()

# The product ladders accumulate the per-node normalizer products of the nested
# forest rule (base and shifted copies) and of the rank-one cluster rule in
# double. The native development probe evaluates a rung on that path and on the
# unchanged long double path (`reference_log_normalizer`, and `reference_value`
# for the shifted copies): they agree to 1e-13 where both apply, and a rung
# whose products leave the double range, whose normalizer fails the checks of
# the batch routine, or whose rows cannot reach it by the up-front skip, is
# evaluated by the long double path itself.

.test_lean_tolerance <- 1e-13

.test_lean_group <- function(count, loading, sd, mean, sei = 1) {

  list(count = count, loading = loading, sd = sd, mean = mean, sei = sei)
}

# A random block of one of the structures the nested rule and the cluster rule
# evaluate; the numbers are drawn, the structure is fixed.
.test_lean_block <- function(kind, omega, z_lower) {

  draw <- function(low, high) stats::runif(1L, low, high)
  group <- function(count, loading) {
    .test_lean_group(count, loading, draw(.6, 1.4), draw(-.4, .6), draw(.3, 1))
  }
  switch(kind,
    cluster = .test_cert_block(list(
      group(sample(2:4, 1L), draw(.3, 1)),
      group(sample(1:3, 1L), draw(-.8, -.2)),
      group(sample(1:2, 1L), draw(.2, .6))), omega, z_lower),
    # a shared axis with two leaves, rows on the shared axis only and a free row
    forest = .test_cert_block(list(
      group(sample(2:4, 1L), c(draw(.3, .8), draw(.3, .8), 0)),
      group(sample(1:3, 1L), c(draw(.3, .8), 0, draw(.3, .8))),
      group(2L, c(draw(.3, .8), 0, 0)),
      group(1L, c(0, 0, 0))), omega, z_lower),
    # nested supports of depth two
    chain = .test_cert_block(list(
      group(2L, c(draw(.3, .8), draw(.3, .8), draw(.2, .6))),
      group(2L, c(draw(.3, .8), draw(.3, .8), 0)),
      group(2L, c(draw(.3, .8), 0, 0))), omega, z_lower))
}

.test_lean_relative <- function(x, reference) abs(expm1(x - reference))

test_that("double node products agree with the long double rules on random blocks", {

  set.seed(20261023)
  omegas <- list(c(1, .6), c(1, .7, .3), c(1, .5, .8, .2), c(1, 1))
  cuts   <- list(c(1.96, -Inf), c(1.96, 0, -Inf), c(1.96, .5, -.5, -Inf), c(1.96, -Inf))
  kinds  <- c("cluster", "forest", "chain")
  difference <- 0
  differing <- 0L
  for (i in seq_len(36L)) {
    kind    <- kinds[(i - 1L) %% 3L + 1L]
    pattern <- (i - 1L) %% 4L + 1L
    block   <- .test_lean_block(kind, omegas[[pattern]], cuts[[pattern]])
    cluster <- kind == "cluster"
    info    <- paste(kind, "omega", pattern, "draw", i)
    for (order in if (cluster) c(15L, 31L) else c(7L, 11L)) {
      axes  <- seq_len(block[["rank"]])
      probe <- .test_cert_probe(block, order, cluster = cluster,
        shift_axis = rep(axes, each = 2L), shift_delta = rep(c(.3, -.2), length(axes)),
        refine = FALSE)
      label <- paste(info, "order", order)
      expect_true(probe[["lean"]], info = label)
      relative <- .test_lean_relative(probe[["log_normalizer"]], probe[["reference_log_normalizer"]])
      expect_lt(relative, .test_lean_tolerance, label = label)
      difference <- max(difference, relative)
      differing  <- differing + (probe[["log_normalizer"]] != probe[["reference_log_normalizer"]])
      # the shifted copies, with the unshifted messages reused and recomputed
      reference <- probe[["reference_value"]]
      known <- is.finite(reference)
      if (!cluster) expect_true(any(known), info = label)
      for (field in c("reuse_value", "full_value")) {
        expect_lt(max(.test_lean_relative(probe[[field]][known], reference[known])), .test_lean_tolerance,
                  label = paste(label, field))
      }
    }
  }
  # not vacuous: the double products are a different arithmetic, not a relabelled copy
  expect_gt(differing, 0L)
  expect_gt(difference, 0)
})

test_that("a rung whose node products leave the double range is evaluated in long double, unchanged", {

  # Rows far in the selection tail: the product of their normalizers (about
  # 1e-57 each) is below 2^-960 at every node, so the double products are
  # abandoned at run time and the long double rule evaluates the rung.
  rows <- 36L
  omega <- c(1, 0)
  cuts  <- c(0, -Inf)
  blocks <- list(
    cluster = .test_cert_block(list(.test_lean_group(rows, .4, .5, -8)), omega, cuts),
    forest  = .test_cert_block(list(
      .test_lean_group(rows / 2L, c(.4, .3, 0), .5, -8),
      .test_lean_group(rows / 2L, c(.4, 0, .3), .5, -8)), omega, cuts))
  for (name in names(blocks)) {
    cluster <- name == "cluster"
    probe <- .test_cert_probe(blocks[[name]], if (cluster) 31L else 9L, cluster = cluster,
                              shift_axis = 1L, shift_delta = .3, refine = FALSE)
    expect_false(probe[["lean"]], info = name)
    expect_true(is.finite(probe[["log_normalizer"]]), info = name)
    expect_lt(probe[["log_normalizer"]], -960 * log(2), label = name)
    expect_identical(probe[["log_normalizer"]], probe[["reference_log_normalizer"]], info = name)
    if (!cluster) {
      expect_identical(probe[["reuse_value"]], probe[["reference_value"]], info = name)
    }
  }

  # A normalizer with massive cancellation (the weight falls from the last bin
  # to the first, and a row sits where the first bin has all the mass) fails
  # the checks of the batch routine, and so does the double path.
  cancelled <- .test_cert_block(list(
    .test_lean_group(2L, .4, .5, 6),
    .test_lean_group(2L, .3, 1.1, 0)), c(0, 1), c(0, -Inf))
  probe <- .test_cert_probe(cancelled, 31L, cluster = TRUE, refine = FALSE)
  expect_false(probe[["lean"]])
  expect_identical(probe[["log_normalizer"]], probe[["reference_log_normalizer"]])
})

test_that("the double path is skipped up front where K log2(smallest positive omega) < -960", {

  # Rows at the lowest weight of the partition: every normalizer is the last
  # weight 2^-24 or a little more, so 40 rows stay in range (2^-960 exactly at a
  # node without selection mass) and 41 rows cannot.
  weight <- 2^-24
  for (rows in c(40L, 41L)) {
    first <- rows %/% 2L
    blocks <- list(
      cluster = .test_cert_block(list(.test_lean_group(rows, .3, .3, -9)), c(1, weight), c(0, -Inf)),
      forest  = .test_cert_block(list(
        .test_lean_group(first, c(.3, .2), .3, -9),
        .test_lean_group(rows - first, c(.3, 0), .3, -9)), c(1, weight), c(0, -Inf)))
    for (name in names(blocks)) {
      cluster <- name == "cluster"
      probe <- .test_cert_probe(blocks[[name]], if (cluster) 31L else 9L, cluster = cluster, refine = FALSE)
      info <- paste(name, rows, "rows")
      expect_identical(probe[["lean"]], rows == 40L, info = info)
      expect_lt(.test_lean_relative(probe[["log_normalizer"]], probe[["reference_log_normalizer"]]),
                .test_lean_tolerance, label = info)
      if (rows > 40L) {
        expect_identical(probe[["log_normalizer"]], probe[["reference_log_normalizer"]], info = info)
      }
    }
  }
})

test_that("tensor rules keep the long double arithmetic, and so do blocks without telescoping", {

  tensor <- .test_cert_block(list(
    .test_lean_group(2L, c(.6, .4, 0), 1.1, .2),
    .test_lean_group(2L, c(0, .5, .4), 1.2, -.1),
    .test_lean_group(2L, c(.5, 0, .6), 1.0, .1)), c(1, .5), c(0, -Inf))
  probe <- .test_cert_probe(tensor, 9L, refine = FALSE)
  expect_identical(probe[["path"]], .test_cert_paths[["tensor_direct"]])
  expect_false(probe[["lean"]])
  expect_identical(probe[["log_normalizer"]], probe[["reference_log_normalizer"]])

  forest <- .test_lean_block("forest", c(1, .6), c(1.96, -Inf))
  probe <- .test_cert_probe(forest, 9L, telescope = FALSE, refine = FALSE)
  expect_false(probe[["lean"]])
  expect_identical(probe[["log_normalizer"]], probe[["reference_log_normalizer"]])
})
