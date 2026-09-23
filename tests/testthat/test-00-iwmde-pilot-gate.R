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


test_that("a stopped line reports the gate and the value it stopped on", {

  # The verdict and the values behind it were recorded on the density curve and
  # then dropped before any consumer saw them. They belong in the public
  # diagnostics, so a user can tell a line whose refinement was cut short from
  # one that settled on its own.
  template <- .iwmde_empty_public_density_diagnostics()
  expect_true(all(c("pilot_gate_stopped", "pilot_bulk_ess") %in% names(template)))
  expect_identical(typeof(template[["pilot_gate_stopped"]]), "logical")
  expect_identical(typeof(template[["pilot_bulk_ess"]]), "double")

  # The gate compares the largest pilot value against the minimum it needs, so
  # that is the one a stopped line has to report.
  expect_identical(.iwmde_pilot_gate_bulk_ess(c(25.9, 24.1)), 25.9)
  expect_identical(.iwmde_pilot_gate_bulk_ess(c(25.9, NA_real_, Inf)), 25.9)
  expect_true(is.na(.iwmde_pilot_gate_bulk_ess(numeric(0))))
  expect_true(is.na(.iwmde_pilot_gate_bulk_ess(NULL)))

  # An estimator that keeps no pilot sequence reports no value rather than a
  # number that would read as a measurement.
  expect_true(is.na(.iwmde_public_logical(NULL)))
  expect_true(is.na(.iwmde_public_numeric(.iwmde_pilot_gate_bulk_ess(NULL))))
})
