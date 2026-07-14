# ============================================================================ #
# test-02-iwmde-glmm-local.R
# ============================================================================ #

context("IWMDE GLMM local-state q-grid equivalence")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-iwmde.R"))


.expect_glmm_local_q_grid_parity <- function(fit_name, outcome_type,
                                             nuisance_parameter) {

  skip_if_missing_fits(fit_name)

  fit     <- load_fit(fit_name)
  context <- .iwmde_context(fit)

  expect_false(
    inherits(fit, "RoBMA"),
    info = paste(fit_name, "must exercise an ordinary GLMM fit")
  )
  expect_equal(
    .data_outcome_type(context[["data"]]),
    outcome_type,
    info = paste(fit_name, "must exercise the expected GLMM family")
  )

  K <- nrow(context[["data"]][["outcome"]])

  for (parameter in c("mu", "tau")) {
    spec      <- .iwmde_parameter_spec(context, parameter, NULL)
    values    <- .iwmde_parameter_values(context, parameter, spec)
    component <- .iwmde_parameter_components(context, parameter, spec)
    rows      <- which(component[["active"]] & is.finite(values))
    rows      <- head(rows, 3L)

    expect_true(
      length(rows) > 0L,
      info = paste(fit_name, parameter, "must have active posterior rows")
    )

    row_states <- .iwmde_row_states(context, rows, parameter, spec)
    keep       <- vapply(row_states, function(state) {
      is.finite(state[["baseline_log_q"]])
    }, logical(1))
    row_states <- row_states[keep]

    expect_true(
      length(row_states) > 0L,
      info = paste(fit_name, parameter, "must have finite local row states")
    )

    for (state in row_states) {
      local_parameters <- c(nuisance_parameter, "theta")

      expect_equal(
        state[["likelihood_mode"]],
        "conditional",
        info = paste(fit_name, parameter, "must use conditional likelihood")
      )
      expect_equal(
        state[["state_scope"]],
        "local",
        info = paste(fit_name, parameter, "must retain local state")
      )
      expect_true(
        all(local_parameters %in% names(state[["parameters"]])),
        info = paste(fit_name, parameter, "must retain local parameters")
      )
      expect_true(
        all(local_parameters %in% names(state[["prior_list"]])),
        info = paste(fit_name, parameter, "must retain local priors")
      )

      for (local_parameter in local_parameters) {
        columns <- paste0(local_parameter, "[", seq_len(K), "]")

        expect_equal(
          state[["parameters"]][[local_parameter]],
          as.numeric(state[["row"]][columns]),
          tolerance = 0,
          info = paste(fit_name, parameter, "must retain sampled", local_parameter)
        )
      }
    }

    active_values <- values[component[["active"]] & is.finite(values)]
    grid_values   <- as.numeric(stats::quantile(
      active_values,
      probs = c(.25, .50, .75),
      names = FALSE,
      type  = 8
    ))
    grid_values <- unique(grid_values[is.finite(grid_values)])

    expect_true(
      length(grid_values) > 0L,
      info = paste(fit_name, parameter, "must provide finite q-grid values")
    )

    replacement <- .iwmde_replacement_spec(context, parameter, spec)
    batch <- .iwmde_log_q_grid_glmm_conditional_batch(
      context     = context,
      parameter   = parameter,
      values      = grid_values,
      row_states  = row_states,
      replacement = replacement
    )
    dispatched <- .iwmde_log_q_grid(
      context     = context,
      parameter   = parameter,
      values      = grid_values,
      row_states  = row_states,
      replacement = replacement
    )
    scalar <- .iwmde_log_q_grid_scalar(
      context     = context,
      parameter   = parameter,
      values      = grid_values,
      row_states  = row_states,
      replacement = replacement
    )

    expect_equal(
      dim(batch),
      c(length(grid_values), length(row_states)),
      info = paste(fit_name, parameter, "conditional batch dimensions")
    )
    expect_equal(
      is.finite(batch),
      is.finite(scalar),
      info = paste(fit_name, parameter, "batch and scalar finite support")
    )

    finite <- is.finite(batch) & is.finite(scalar)
    expect_true(
      any(finite),
      info = paste(fit_name, parameter, "must have finite q-grid terms")
    )
    expect_equal(
      batch[finite],
      scalar[finite],
      tolerance = 1e-8,
      info = paste(fit_name, parameter, "conditional batch versus scalar")
    )
    expect_equal(
      dispatched,
      batch,
      tolerance = 1e-8,
      info = paste(fit_name, parameter, "q-grid dispatch versus GLMM batch")
    )
  }
}


test_that("binomial GLMM q-grid retains local states and matches scalar", {

  .expect_glmm_local_q_grid_parity(
    fit_name           = "bcg_glmm",
    outcome_type       = "bin",
    nuisance_parameter = "pi"
  )
})


test_that("Poisson GLMM q-grid retains local states and matches scalar", {

  .expect_glmm_local_q_grid_parity(
    fit_name           = "nielweise2008_glmm",
    outcome_type       = "pois",
    nuisance_parameter = "phi"
  )
})
