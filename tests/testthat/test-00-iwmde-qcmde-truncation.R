context("(00) qCMDE truncation reporting")

# Rows whose Gaussian likelihood kernel meets a normal prior are truncated
# normals: their conditional law, the mass outside the normalization range,
# and the ordinate all follow from the normal distribution directly.
.qcmde_truncation_rows <- function() {

  list(
    center = c(-.4, .1, .35),
    scale  = c(.15, .08, .2),
    prior  = BayesTools::prior("normal", list(mean = 0, sd = 1))
  )
}


.qcmde_truncation_density <- function(rows, probability, value = .05) {

  current <- rows[["center"]] + .02
  testthat::local_mocked_bindings(
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {

      out <- vapply(seq_along(rows[["center"]]), function(row) {
        -.5 * ((values - rows[["center"]][[row]]) / rows[["scale"]][[row]])^2 +
          stats::dnorm(values, 0, 1, log = TRUE)
      }, numeric(length(values)))
      out <- matrix(out, nrow = length(values))
      attr(out, "gaussian_kernel") <- list(
        current     = current,
        linear      = (rows[["center"]] - current) / rows[["scale"]]^2,
        quadratic   = 1 / rows[["scale"]]^2,
        prior_route = rep("focal", length(current))
      )
      out
    },
    .package = "RoBMA"
  )
  states <- lapply(seq_along(rows[["center"]]), function(row) {
    list(row_index = row, focal_prior = rows[["prior"]],
         use_focal_prior_delta = TRUE)
  })
  base <- seq(-.2, .2, length.out = 101L)

  .iwmde_density_grid(
    context            = list(),
    parameter          = "mu",
    display_grid       = value,
    normalization_grid = list(x = base, z = base, log_jacobian = rep(0, 101L)),
    transform          = .iwmde_parameter_transform(c(-Inf, Inf)),
    row_states         = states,
    active_mass        = 1,
    replacement        = list(type = "scalar"),
    normalization_prob = probability
  )
}


.qcmde_truncation_laws <- function(rows) {

  precision <- 1 / rows[["scale"]]^2 + 1
  list(
    mean = rows[["center"]] / rows[["scale"]]^2 / precision,
    sd   = 1 / sqrt(precision)
  )
}


test_that("exact rows report their truncation and bound the ordinate error", {

  rows <- .qcmde_truncation_rows()
  laws <- .qcmde_truncation_laws(rows)
  for (probability in c(.99, .999, 1 - 1e-8)) {
    density <- .qcmde_truncation_density(rows, probability)
    range   <- density[["normalization_range"]]

    # The range is the union of the rows' central intervals; the draws range
    # is not used when every row is exact.
    quantile <- stats::qnorm((1 - probability) / 2, lower.tail = FALSE)
    expect_equal(range, c(min(laws$mean - quantile * laws$sd),
                          max(laws$mean + quantile * laws$sd)),
                 tolerance = 1e-10)

    outside <- stats::pnorm(range[[1L]], laws$mean, laws$sd) +
      stats::pnorm(range[[2L]], laws$mean, laws$sd, lower.tail = FALSE)
    expect_identical(density[["normalization_truncation_status"]], "exact")
    expect_equal(density[["normalization_truncation"]], max(outside),
                 tolerance = 1e-8)
    expect_lte(density[["normalization_truncation"]], 1 - probability + 1e-12)
    expect_equal(density[["truncation_ordinate_bound"]],
                 max(outside) / (1 - max(outside)), tolerance = 1e-8)

    # The ordinate error is the truncation's overstatement plus the
    # discretization the nested grids measure.
    truth <- mean(stats::dnorm(.05, laws$mean, laws$sd))
    error <- density[["y"]] / truth - 1
    expect_lte(error, density[["truncation_ordinate_bound"]] +
      density[["ordinate_relative_change"]])
    expect_gte(error, -density[["ordinate_relative_change"]])
  }

  # The coverage is what sets the result: a lower 'normalization_prob' moves
  # the ordinate by about the truncation it reports.
  loose <- .qcmde_truncation_density(rows, .99)
  tight <- .qcmde_truncation_density(rows, 1 - 1e-8)
  expect_gt(loose[["y"]] / tight[["y"]] - 1, 1e-4)
  expect_lt(loose[["y"]] / tight[["y"]] - 1,
            loose[["truncation_ordinate_bound"]])
})


test_that("the qCMDE ordinate error bound adds the grid change and truncation", {

  diagnostics <- function(change, bound, status = "estimate") {
    list(
      estimator                       = "q_grid_cmde",
      ordinate                        = .4,
      ordinate_relative_change        = change,
      normalization_truncation        = bound / (1 + bound),
      normalization_truncation_status = status,
      truncation_ordinate_bound       = bound
    )
  }

  expect_equal(
    .iwmde_diagnostics_qcmde_ordinate_relative_change(diagnostics(.01, .02)),
    .03
  )
  expect_identical(
    .iwmde_diagnostics_bf_failure_reason(diagnostics(.001, .06)),
    paste0(
      "qCMDE posterior ordinate error bound is 6.1% (grid change 0.1%; ",
      "tail truncation 6%, estimated). Try setting 'normalization_prob' ",
      "closer to 1 in the 'density_control' argument"
    )
  )
  expect_identical(
    .iwmde_diagnostics_bf_failure_reason(diagnostics(.05, .01, "exact")),
    paste0(
      "qCMDE posterior ordinate error bound is 6% (grid change 5%; ",
      "tail truncation 1%, exact). Try increasing 'normalization_points' in ",
      "the 'density_control' argument"
    )
  )
  expect_null(.iwmde_diagnostics_bf_failure_reason(diagnostics(.01, .01)))
  expect_identical(
    .iwmde_diagnostics_bf_warning(diagnostics(.01, .02, "bound")),
    paste0(
      "qCMDE posterior ordinate error bound is 3% (grid change 1%; ",
      "tail truncation 2%, upper bound) (warning threshold ",
      .iwmde_percent(.025), "; BF rejection threshold ", .iwmde_percent(.05),
      ")."
    )
  )

  # A row that does not decay at an end of the range has no finite bound; a
  # reported but unknown truncation leaves the check unavailable.
  expect_identical(
    .iwmde_diagnostics_bf_failure_reason(diagnostics(.001, Inf)),
    paste0(
      "qCMDE posterior ordinate error is unbounded: a conditional density ",
      "does not decrease at an end of the normalization range. Inspect ",
      "density_diagnostics() and try another 'density_method'"
    )
  )
  expect_identical(
    .iwmde_diagnostics_bf_failure_reason(diagnostics(.001, NA_real_)),
    "normalization diagnostics are unavailable"
  )
})


test_that("density diagnostics report the truncation and its status", {

  entry <- BayesTools::posterior_ordinate_attribute(
    value          = 0,
    ordinate       = .4,
    method         = "q_grid_cmde",
    density_method = "qCMDE",
    diagnostics    = list(
      evaluation_value                = 0,
      relative_mcse                   = .01,
      finite_terms                    = 200,
      ess                             = 150,
      max_weight_share                = .05,
      active_mass                     = 1,
      normalization_relative_error    = 0,
      ordinate_relative_change        = .002,
      normalization_truncation        = .0099,
      normalization_truncation_status = "bound",
      truncation_ordinate_bound       = .01,
      estimator                       = "q_grid_cmde",
      weight_method                   = "conditional_grid"
    )
  )
  row <- .iwmde_public_density_diagnostic_row(entry)

  expect_identical(row[["stability_metric"]], "ordinate_error_bound")
  expect_equal(row[["stability_relative_error"]], .012)
  expect_equal(row[["ordinate_relative_change"]], .002)
  expect_equal(row[["normalization_truncation"]], .0099)
  expect_identical(row[["normalization_truncation_status"]], "bound")
  expect_equal(row[["truncation_ordinate_bound"]], .01)
  expect_identical(row[["status"]], "ok")
  expect_identical(names(row),
                   names(.iwmde_empty_public_density_diagnostics()))

  # Records from before the truncation was reported gain typed missing values.
  older <- row[setdiff(names(row), c("normalization_truncation",
    "normalization_truncation_status", "truncation_ordinate_bound"))]
  class(older) <- c("RoBMA_density_diagnostics", "data.frame")
  validated <- .density_diagnostics_validate(older)
  expect_identical(names(validated),
                   names(.iwmde_empty_public_density_diagnostics()))
  expect_identical(validated[["normalization_truncation_status"]],
                   NA_character_)
  expect_identical(validated[["truncation_ordinate_bound"]], NA_real_)
})


test_that("qCMDE requires a normalization probability below 1", {

  expect_error(
    .density_control_normalize("qCMDE", list(normalization_prob = 1)),
    "'density_control$normalization_prob' must be lower than 1 for qCMDE.",
    fixed = TRUE
  )
  expect_equal(
    .density_control_normalize("IWMDE",
                               list(normalization_prob = 1))[["normalization_prob"]],
    1
  )
})
