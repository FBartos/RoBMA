context("(00) qCMDE pilot acceptance gate")

# The qCMDE pilot gate stopped the old refinement sequence of widening grids
# once a density line's bulk effective sample size had settled far below the
# acceptance gate. The nested normalization grids evaluate exactly one selected
# grid and its validation grid, so there is no refinement for the gate to cut
# short; the public columns keep reporting that no line was stopped.

test_that("batch ESS counts rounding-level variation as constant and keeps real variation", {

  n <- 240L
  batch <- function(x) .iwmde_batch_mcse(matrix(x, nrow = 1L))
  withr::local_seed(20261109L)
  noise <- stats::rnorm(n)

  # exactly constant contributions, and constant up to rounding at any scale: the
  # Kish value, where the noise alone used to estimate an inefficiency of 2 to 7
  expect_lt(abs(batch(rep(.7, n))[["ess"]] / n - 1), 1e-12)
  for (scale in c(1e-100, 1e-10, .7, 1e100)) {
    x <- scale * (1 + 3e-16 * noise)
    expect_gt(stats::var(x), 0)
    expect_lt(abs(batch(x)[["ess"]] / n - 1), 1e-12, label = paste("rounding noise at scale", scale))
  }
  # the largest level measured across platforms
  expect_lt(abs(batch(1 + 1.35e-15 * noise)[["ess"]] / n - 1), 1e-12)

  # real variation keeps its estimated inefficiency down to small scales, and the
  # floor follows the size of the contributions
  chain <- as.numeric(stats::arima.sim(list(ar = .9), n))
  for (level in c(.1, 1e-6, 1e-12, 1e-14)) {
    for (scale in c(1, 1e-100)) {
      x <- scale * (1 + level * chain)
      result <- batch(x)
      kish <- sum(x)^2 / sum(x^2)
      expect_lt(abs(result[["ess"]] / (kish / max(1, n * result[["mcse"]]^2 / stats::var(x))) - 1), 1e-12,
                label = paste("level", level, "scale", scale))
      expect_lt(result[["ess"]], n / 2, label = paste("level", level, "scale", scale))
    }
  }

  # the MCSE, which the densities' uncertainty reports, does not change: a mixed matrix
  # gives each row its own decision
  mixed <- rbind(1 + 3e-16 * noise, 1 + 1e-12 * chain)
  result <- .iwmde_batch_mcse(mixed)
  expect_lt(abs(result[["ess"]][[1L]] / n - 1), 1e-12)
  expect_lt(result[["ess"]][[2L]], n / 2)
  expect_gt(result[["mcse"]][[1L]], 0)
})


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
