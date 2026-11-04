context("Selection refinement counts per axis")
skip_on_cran()

test_that("a small leaf keeps its own row count in the tier gate", {

  block <- .test_cert_block(list(
    list(count = 80L, loading = c(.2, .2, 0), sd = 1, mean = 0, sei = 1),
    list(count = 3L, loading = c(.01, 0, 1), sd = .08, mean = 0, sei = 1)),
    c(1, .1))
  probe <- .test_cert_probe(block, 95L, refine = FALSE)
  expect_equal(.test_cert_tier_width_scales(block),
    8 / sqrt(2 * log(c(83, 80, 3))), tolerance = 1e-14)
  expect_false(.test_cert_resolves(block, 95L, .test_cert_tier_width_scale(block)))
  expect_true(.test_cert_tier_resolves(block, 95L))
  expect_true(probe[["resolves_tier"]])
  expect_equal(probe[["tier_width_scale"]], .test_cert_tier_width_scale(block), tolerance = 1e-14)

  # Reordering and reversing axes leaves their row counts and geometry intact.
  block[["loading"]] <- block[["loading"]][, c(3L, 1L, 2L)] * rep(c(-1, 1, -1), each = block[["rows"]])
  expect_true(.test_cert_probe(block, 95L, refine = FALSE)[["resolves_tier"]])
})

test_that("shared rows retain the full count on every axis", {

  order <- 95L
  nodes <- .gauss_hermite_nodes(order)[["nodes"]]
  gap <- max(diff(nodes)[utils::head(nodes, -1L) <= 3 & utils::tail(nodes, -1L) >= -3])
  scale <- 8 / sqrt(2 * log(80))
  for (side in c(-1, 1)) {
    block <- .test_cert_block(list(
      list(count = 80L, loading = c(1, 1), sd = gap / scale * (1 + side * 1e-6),
           mean = 0, sei = 1)), c(1, .1))
    probe <- .test_cert_probe(block, order, refine = FALSE)
    expect_equal(.test_cert_tier_width_scales(block), rep(scale, 2), tolerance = 1e-14)
    expect_identical(probe[["resolves_tier"]], side > 0)
    expect_identical(.test_cert_tier_resolves(block, order), side > 0)
    expect_identical(.test_cert_resolves(block, order, scale), side > 0)
  }
})

test_that("a forest with unequal leaf counts can return its accurate refined value", {

  block <- .test_cert_block(list(
    list(count = 24L, loading = c(.396212124289880, .598265364819431, 0, 0),
         sd = .166745965827061, mean = 2.34033728058624, sei = 1),
    list(count = 15L, loading = c(.360349379222890, 0, .434351429645156, 0),
         sd = .0449433707043086, mean = 2.25913594757555, sei = 1),
    list(count = 2L, loading = c(-1.25006768340786, 0, 0, .131449513438867),
         sd = .244317296614989, mean = 2.62623363831159, sei = 1)),
    c(1, .5), c(1.96, -Inf))
  reference <- .test_cert_reference_forest(block)
  tier <- .test_cert_tier(block)
  expect_true(is.na(.test_cert_probe_ladder(block)))
  expect_false(.test_cert_resolves(block, 95L, .test_cert_tier_width_scale(block)))
  expect_true(.test_cert_tier_resolves(block, 95L))
  expect_true(tier[["accepted"]])
  expect_lt(.test_cert_relative_error(tier[["value"]], reference), 1e-6)
  actual <- .test_cert_native_factor(block)
  expect_true(.test_cert_deterministic(actual))
  expect_equal(actual[["relative_mcse"]], 0)
  expect_equal(as.numeric(actual[["log_normalizer"]]), tier[["value"]], tolerance = 1e-13)
  expect_lt(.test_cert_relative_error(actual[["log_normalizer"]], reference), .005)
})

test_that("the tier mirror rejects a rule without node gaps like the native gate", {

  # A single node resolves no transition and has no central gap; the native rule
  # returns false for order < 2 whenever the weights change, and accepts a block
  # without transitions whatever the rule.
  for (loading in list(c(1, 1), c(1, 0))) {
    info <- paste("loading", paste(loading, collapse = ","))
    block <- .test_cert_block(list(
      list(count = 5L, loading = loading, sd = 1, mean = 0, sei = 1)), c(1, .1))
    expect_false(.test_cert_tier_resolves(block, 1L), info = info)
    expect_true(.test_cert_tier_resolves(block, 3L), info = info)
    constant <- .test_cert_block(list(
      list(count = 5L, loading = loading, sd = 1, mean = 0, sei = 1)), c(1, 1))
    expect_true(.test_cert_tier_resolves(constant, 1L), info = info)
  }
})
