context("(00) qCMDE pilot acceptance gate")

# The qCMDE pilot gate stopped the old refinement sequence of widening grids
# once a density line's bulk effective sample size had settled far below the
# acceptance gate. The nested normalization grids evaluate exactly one selected
# grid and its validation grid, so there is no refinement for the gate to cut
# short; the public columns keep reporting that no line was stopped.

test_that("nested qCMDE grids report no pilot stop", {

  normalization_values <- seq(-6, 6, length.out = 25L)
  testthat::local_mocked_bindings(
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {
      matrix(stats::dnorm(values, log = TRUE), nrow = length(values),
             ncol = length(row_states))
    },
    .package = "RoBMA"
  )
  density <- .iwmde_density_grid(
    context            = list(),
    parameter          = "mu",
    display_grid       = seq(-3, 3, length.out = 31L),
    normalization_grid = list(
      x            = normalization_values,
      z            = normalization_values,
      log_jacobian = rep(0, length(normalization_values))
    ),
    transform          = .iwmde_parameter_transform(c(-Inf, Inf)),
    row_states         = rep(list(list()), 5L),
    active_mass        = 1,
    replacement        = list(type = "scalar")
  )
  expect_false(density[["pilot_gate_stopped"]])
  expect_length(density[["pilot_bulk_ess"]], 0L)
  expect_true(is.na(.iwmde_public_numeric(
    .iwmde_pilot_gate_bulk_ess(density[["pilot_bulk_ess"]])
  )))
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
