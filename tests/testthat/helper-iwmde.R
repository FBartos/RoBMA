# ============================================================================ #
# helper-iwmde.R
# ============================================================================ #

.expect_iwmde_ok <- function(out, parameters) {

  expect_named(out, parameters)

  for (parameter in parameters) {
    diagnostic <- out[[parameter]]

    expect_equal(diagnostic[["status"]], "ok")
    expect_length(diagnostic[["iwmde"]][["x"]], 20L)
    expect_length(diagnostic[["iwmde"]][["y"]], 20L)
    expect_true(all(is.finite(diagnostic[["iwmde"]][["y"]])))
    expect_true(all(diagnostic[["iwmde"]][["y"]] >= 0))
    expect_true(is.finite(diagnostic[["diagnostics"]][["integral"]]))
    expect_gt(diagnostic[["diagnostics"]][["integral"]], 0)
    expect_gt(diagnostic[["diagnostics"]][["display_mass_fraction"]], 0)
    expect_true(is.finite(diagnostic[["diagnostics"]][["display_mass_fraction"]]))
    expect_lt(diagnostic[["diagnostics"]][["display_mass_fraction"]], 10)
    expect_equal(
      diagnostic[["diagnostics"]][["normalization_integral"]],
      diagnostic[["active_mass"]],
      tolerance = 1e-6
    )
    expect_true(is.finite(diagnostic[["diagnostics"]][["min_ess"]]))
    expect_gt(diagnostic[["diagnostics"]][["min_ess"]], 0)
    expect_true(
      diagnostic[["diagnostics"]][["max_weight_share"]] >= 0 &&
        diagnostic[["diagnostics"]][["max_weight_share"]] <= 1
    )
    expect_true(is.finite(diagnostic[["diagnostics"]][["normalization_points"]]))
    expect_gt(diagnostic[["diagnostics"]][["normalization_points"]], 0)
    expect_true(is.finite(diagnostic[["diagnostics"]][["n_batches"]]))
    expect_equal(diagnostic[["diagnostics"]][["estimator"]], "q_grid_cmde")
    expect_equal(
      diagnostic[["active_mass"]] +
        sum(diagnostic[["point_masses"]][["mass"]]),
      1,
      tolerance = 1e-8
    )
  }
}


.skip_if_missing_raw_fits <- function(names) {

  skip_if_missing_fits(names)
}

.expect_iwmde_batch_equals_scalar <- function(fit_name, parameter,
                                             parameter_spec = NULL,
                                             n_rows = 3L,
                                             n_values = 3L) {

  .skip_if_missing_raw_fits(fit_name)

  context <- .iwmde_context(load_fit(fit_name, validate = FALSE))
  spec    <- .iwmde_parameter_spec(context, parameter, parameter_spec)
  values  <- .iwmde_parameter_values(context, parameter, spec)
  finite  <- is.finite(values)
  component <- .iwmde_parameter_components(context, parameter, spec)
  rows      <- which(component[["active"]] & finite)
  rows      <- head(rows, n_rows)

  expect_gt(length(rows), 0L)

  row_states <- .iwmde_row_states(context, rows, parameter, spec)
  keep       <- vapply(row_states, function(state) {
    is.finite(state[["baseline_log_q"]])
  }, logical(1))
  row_states <- row_states[keep]
  rows       <- rows[keep]

  expect_gt(length(rows), 0L)

  active_values <- values[component[["active"]] & finite]
  grid_values   <- as.numeric(stats::quantile(
    active_values,
    probs = seq(.25, .75, length.out = n_values),
    names = FALSE,
    type  = 8
  ))
  grid_values <- unique(grid_values[is.finite(grid_values)])

  expect_gt(length(grid_values), 0L)

  replacement <- .iwmde_replacement_spec(context, parameter, spec)
  batch       <- .iwmde_log_q_grid(
    context     = context,
    parameter   = parameter,
    values      = grid_values,
    row_states  = row_states,
    replacement = replacement
  )
  scalar      <- .iwmde_log_q_grid_scalar(
    context     = context,
    parameter   = parameter,
    values      = grid_values,
    row_states  = row_states,
    replacement = replacement
  )

  expect_equal(dim(batch), dim(scalar))
  expect_equal(is.finite(batch), is.finite(scalar))
  finite_terms <- is.finite(batch) & is.finite(scalar)
  if (any(finite_terms)) {
    expect_equal(batch[finite_terms], scalar[finite_terms], tolerance = 1e-8)
  }
}


.expect_iwmde_predictor_fast_equals_scalar <- function(fit_name, parameter,
                                                       parameter_spec = NULL,
                                                       n_rows = 3L,
                                                       n_values = 3L) {

  .skip_if_missing_raw_fits(fit_name)

  context <- .iwmde_context(load_fit(fit_name, validate = FALSE))
  spec    <- .iwmde_parameter_spec(context, parameter, parameter_spec)
  values  <- .iwmde_parameter_values(context, parameter, spec)
  finite  <- is.finite(values)
  component <- .iwmde_parameter_components(context, parameter, spec)
  rows      <- which(component[["active"]] & finite)
  rows      <- head(rows, n_rows)

  expect_gt(length(rows), 0L)

  row_states <- .iwmde_row_states(context, rows, parameter, spec)
  keep       <- vapply(row_states, function(state) {
    is.finite(state[["baseline_log_q"]])
  }, logical(1))
  row_states <- row_states[keep]

  expect_gt(length(row_states), 0L)

  active_values <- values[component[["active"]] & finite]
  grid_values   <- as.numeric(stats::quantile(
    active_values,
    probs = seq(.25, .75, length.out = n_values),
    names = FALSE,
    type  = 8
  ))
  grid_values <- unique(grid_values[is.finite(grid_values)])

  expect_gt(length(grid_values), 0L)

  replacement <- .iwmde_replacement_spec(context, parameter, spec)
  fast        <- .iwmde_log_q_grid_predictor_batch(
    context     = context,
    parameter   = parameter,
    values      = grid_values,
    row_states  = row_states,
    replacement = replacement
  )
  scalar      <- .iwmde_log_q_grid_scalar(
    context     = context,
    parameter   = parameter,
    values      = grid_values,
    row_states  = row_states,
    replacement = replacement
  )

  expect_true(is.matrix(fast))
  expect_equal(dim(fast), dim(scalar))
  expect_equal(is.finite(fast), is.finite(scalar))
  finite_terms <- is.finite(fast) & is.finite(scalar)
  if (any(finite_terms)) {
    expect_equal(fast[finite_terms], scalar[finite_terms], tolerance = 1e-8)
  }
}


.expect_iwmde_normal_quadratic_equals_scalar <- function(fit_name, parameter,
                                                         parameter_spec = NULL,
                                                         n_rows = 3L,
                                                         n_values = 3L) {

  .skip_if_missing_raw_fits(fit_name)

  context <- .iwmde_context(load_fit(fit_name, validate = FALSE))
  spec    <- .iwmde_parameter_spec(context, parameter, parameter_spec)
  values  <- .iwmde_parameter_values(context, parameter, spec)
  finite  <- is.finite(values)
  component <- .iwmde_parameter_components(context, parameter, spec)
  rows      <- which(component[["active"]] & finite)
  rows      <- head(rows, n_rows)

  expect_gt(length(rows), 0L)

  row_states <- .iwmde_row_states(context, rows, parameter, spec)
  keep       <- vapply(row_states, function(state) {
    is.finite(state[["baseline_log_q"]])
  }, logical(1))
  row_states <- row_states[keep]

  expect_gt(length(row_states), 0L)

  active_values <- values[component[["active"]] & finite]
  grid_values   <- as.numeric(stats::quantile(
    active_values,
    probs = seq(.25, .75, length.out = n_values),
    names = FALSE,
    type  = 8
  ))
  grid_values <- unique(grid_values[is.finite(grid_values)])

  expect_gt(length(grid_values), 0L)

  replacement  <- .iwmde_replacement_spec(context, parameter, spec)
  active_setup <- row_states[[1L]][["active_setup"]]
  unit         <- if (.is_data_multilevel(context[["data"]])) {
    "cluster"
  } else {
    "estimate"
  }
  setup <- .iwmde_predictor_setup(
    context      = context,
    row_states   = row_states,
    active_setup = active_setup,
    unit         = unit
  )
  basis <- .iwmde_predictor_update_basis(
    context     = context,
    parameter   = parameter,
    row_states  = row_states,
    replacement = replacement,
    setup       = setup
  )
  basis <- .iwmde_predictor_resolve_formula_basis(
    context      = context,
    active_setup = active_setup,
    setup        = setup,
    basis        = basis,
    row_states   = row_states
  )

  fast <- .iwmde_log_q_grid_normal_location_group(
    context     = context,
    parameter   = parameter,
    values      = grid_values,
    row_states  = row_states,
    replacement = replacement,
    setup       = setup,
    basis       = basis
  )
  scalar <- .iwmde_log_q_grid_scalar(
    context     = context,
    parameter   = parameter,
    values      = grid_values,
    row_states  = row_states,
    replacement = replacement
  )

  expect_true(is.matrix(fast))
  expect_equal(dim(fast), dim(scalar))
  expect_equal(is.finite(fast), is.finite(scalar))
  finite_terms <- is.finite(fast) & is.finite(scalar)
  if (any(finite_terms)) {
    expect_equal(fast[finite_terms], scalar[finite_terms], tolerance = 1e-8)
  }
}
