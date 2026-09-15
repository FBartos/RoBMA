context("(00) qCMDE pilot acceptance gate")

# The qCMDE estimator refines its normalization grid until the normalizer
# settles. That refinement is the expensive part of a density line, and it
# cannot change the acceptance gate's verdict when the line's row weights are
# already far below the gate: the grid refines a per-row constant, while the
# gate is driven by the effective sample size those weights produce. These
# tests pin the predicate that stops the refinement.

test_that("the pilot stop needs both a settled sequence and a wide margin", {

  rows    <- 500L
  minimum <- .iwmde_density_min_ess(rows)
  margin  <- .iwmde_qcmde_pilot_gate_margin()
  expect_gt(margin, 1)
  expect_gt(minimum, 0)

  hopeless <- minimum / margin / 2

  # Settled and far below the gate: the refinement cannot rescue the line.
  expect_true(.iwmde_qcmde_pilot_gate_hopeless(
    c(hopeless, hopeless * 1.05), rows))

  # Far below the gate but still moving: the sequence has not settled, so the
  # refinement is allowed to continue.
  expect_false(.iwmde_qcmde_pilot_gate_hopeless(
    c(hopeless, hopeless * margin * 1.2), rows))

  # Settled but near the gate: a refinement could still carry it across.
  expect_false(.iwmde_qcmde_pilot_gate_hopeless(
    c(minimum * 0.95, minimum), rows))

  # Comfortably accepted lines are never stopped.
  expect_false(.iwmde_qcmde_pilot_gate_hopeless(
    c(minimum * 3, minimum * 3.2), rows))

  # One grid, or a non-finite value, is never enough to stop.
  expect_false(.iwmde_qcmde_pilot_gate_hopeless(hopeless, rows))
  expect_false(.iwmde_qcmde_pilot_gate_hopeless(c(hopeless, NA_real_), rows))
  expect_false(.iwmde_qcmde_pilot_gate_hopeless(numeric(0), rows))
})


test_that("the pilot bulk ESS reads the central mass of the pilot curve", {

  # Three rows with very different weights on a grid whose curve is unimodal.
  display_grid <- seq(-4, 4, length.out = 65L)
  log_q <- vapply(c(0, -0.4, 0.3), function(shift) {
    stats::dnorm(display_grid, shift, 1, log = TRUE)
  }, numeric(length(display_grid)))

  balanced <- .iwmde_qcmde_pilot_bulk_ess(
    display_grid   = display_grid,
    log_q_display  = log_q,
    log_normalizer = rep(0, 3L),
    active_mass    = 1,
    denominator    = 3L
  )
  expect_true(is.finite(balanced))
  expect_gt(balanced, 1)
  expect_lte(balanced, 3)

  # One row given an overwhelmingly larger weight must drive the effective
  # sample size down towards one.
  dominated <- .iwmde_qcmde_pilot_bulk_ess(
    display_grid   = display_grid,
    log_q_display  = log_q,
    log_normalizer = c(-12, 0, 0),
    active_mass    = 1,
    denominator    = 3L
  )
  expect_true(is.finite(dominated))
  expect_lt(dominated, balanced)
  expect_lt(dominated, 1.1)

  # Inputs the estimator cannot summarize report no pilot value rather than a
  # number the stop would act on.
  expect_true(is.na(.iwmde_qcmde_pilot_bulk_ess(
    display_grid = display_grid, log_q_display = log_q,
    log_normalizer = c(0, NA_real_, 0), active_mass = 1, denominator = 3L)))
  expect_true(is.na(.iwmde_qcmde_pilot_bulk_ess(
    display_grid = display_grid[1:2], log_q_display = log_q,
    log_normalizer = rep(0, 3L), active_mass = 1, denominator = 3L)))
})
