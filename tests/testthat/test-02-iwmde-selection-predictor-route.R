# ============================================================================ #
# test-02-iwmde-selection-predictor-route.R
# ============================================================================ #

context("IWMDE selection predictor route")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-iwmde.R"))


# Every cached fixture that may carry a joint selection model; the guard of
# .iwmde_log_q_grid() decides which of them actually take the predictor route.
.selection_route_fit_names <- function() {

  unique(c(
    list_fits(feature = "selection"),
    list_fits(class = c("bselmodel", "RoBMA", "bselmodel.mv", "RoBMA.mv"))
  ))
}


# One location target, one moderator slope where present, one heterogeneity
# target and every bias coordinate the posterior carries.
.selection_route_parameters <- function(context) {

  columns    <- colnames(context[["posterior_samples"]])
  location   <- intersect(c("mu", "mu_intercept"), columns)
  moderators <- setdiff(grep("^mu_", columns, value = TRUE), "mu_intercept")
  scale      <- intersect(c("tau", "log_tau_intercept"), columns)
  bias       <- intersect(c("PET", "PEESE"), columns)
  weights    <- grep("^omega\\[", columns, value = TRUE)

  candidates <- c(
    utils::head(location, 1L),
    utils::head(moderators, 1L),
    utils::head(scale, 1L),
    bias,
    utils::tail(weights, 1L)
  )

  return(unique(candidates))
}


# The dispatch guard of .iwmde_log_q_grid(): a joint selection model whose
# execution plan has at least one block and only singleton blocks.
.selection_route_qualifies <- function(context) {

  if (!.is_data_joint_selection(context[["data"]])) {
    return(FALSE)
  }
  plan <- .data_selection_execution_plan(context[["data"]])

  return(length(plan[["block_methods"]]) > 0L &&
           all(plan[["block_methods"]] == "singleton"))
}


# Rows and grid values for one target, or NULL when the target is unavailable.
.selection_route_inputs <- function(context, parameter, n_rows = 6L,
                                    n_values = 4L) {

  spec <- tryCatch(
    .iwmde_parameter_spec(context, parameter, NULL),
    error = function(e) NULL
  )
  if (is.null(spec) || !identical(spec[["status"]], "ok")) {
    return(NULL)
  }
  values <- tryCatch(
    .iwmde_parameter_values(context, parameter, spec),
    error = function(e) NULL
  )
  component <- tryCatch(
    .iwmde_parameter_components(context, parameter, spec),
    error = function(e) NULL
  )
  if (is.null(values) || is.null(component)) {
    return(NULL)
  }

  active <- component[["active"]] & is.finite(values)
  if (sum(active) < 2L) {
    return(NULL)
  }
  rows <- utils::head(which(active), n_rows)

  row_states <- tryCatch(
    .iwmde_row_states(context, rows, parameter, spec),
    error = function(e) NULL
  )
  if (is.null(row_states)) {
    return(NULL)
  }
  finite_baseline <- vapply(row_states, function(state) {
    is.finite(state[["baseline_log_q"]])
  }, logical(1))
  row_states <- row_states[finite_baseline]
  if (length(row_states) == 0L) {
    return(NULL)
  }

  grid_values <- as.numeric(stats::quantile(
    values[active],
    probs = seq(.2, .8, length.out = n_values),
    names = FALSE,
    type  = 8
  ))
  grid_values <- unique(grid_values[is.finite(grid_values)])
  if (length(grid_values) == 0L) {
    return(NULL)
  }

  return(list(
    row_states  = row_states,
    values      = grid_values,
    replacement = .iwmde_replacement_spec(context, parameter, spec)
  ))
}


.expect_log_q_agreement <- function(observed, reference, label) {

  expect_equal(dim(observed), dim(reference), info = label)
  expect_identical(
    observed == -Inf,
    reference == -Inf,
    info = paste0(label, ": -Inf pattern")
  )
  expect_identical(is.finite(observed), is.finite(reference), info = label)

  finite <- is.finite(reference)
  expect_true(any(finite), info = paste0(label, ": no finite log density"))
  # Relative agreement, absolute where the reference is below one in magnitude.
  deviation <- max(abs(observed[finite] - reference[finite]) /
                     pmax(abs(reference[finite]), 1))
  expect_lt(deviation, 1e-10, label = paste0(label, ": max scaled deviation"))
}


test_that("all-singleton selection models take the predictor route", {

  fit_names <- .selection_route_fit_names()
  skip_if(length(fit_names) == 0L, "No cached selection fixtures are available.")

  checked  <- character()
  declined <- character()

  for (fit_name in fit_names) {
    context <- .iwmde_context(load_fit(fit_name, validate = FALSE))
    if (!.selection_route_qualifies(context)) {
      next
    }

    for (parameter in .selection_route_parameters(context)) {
      inputs <- .selection_route_inputs(context, parameter)
      if (is.null(inputs)) {
        next
      }
      label <- paste0(fit_name, " / ", parameter)

      predictor <- .iwmde_log_q_grid_predictor_batch(
        context     = context,
        parameter   = parameter,
        values      = inputs[["values"]],
        row_states  = inputs[["row_states"]],
        replacement = inputs[["replacement"]]
      )
      if (!is.matrix(predictor) ||
          nrow(predictor) != length(inputs[["values"]]) ||
          ncol(predictor) != length(inputs[["row_states"]])) {
        # The batch declined; the generic route keeps serving this target.
        declined <- c(declined, label)
        next
      }

      generic <- local({
        testthat::local_mocked_bindings(
          .iwmde_log_q_grid_predictor_batch = function(...) NULL,
          .package = "RoBMA"
        )
        .iwmde_log_q_grid(
          context     = context,
          parameter   = parameter,
          values      = inputs[["values"]],
          row_states  = inputs[["row_states"]],
          replacement = inputs[["replacement"]]
        )
      })
      .expect_log_q_agreement(predictor, generic, label)

      dispatched <- .iwmde_log_q_grid(
        context     = context,
        parameter   = parameter,
        values      = inputs[["values"]],
        row_states  = inputs[["row_states"]],
        replacement = inputs[["replacement"]]
      )
      expect_equal(
        as.numeric(dispatched),
        as.numeric(predictor),
        tolerance = 0,
        info      = paste0(label, ": dispatch did not take the predictor route")
      )

      checked <- c(checked, label)
    }
  }

  # The fixtures that reach the route are data, not an assumption: record them
  # so a catalog change that removes them fails here instead of passing empty.
  expect_gt(length(checked), 0L)
  expect_true(any(grepl("3PSM", checked)), info = paste(checked, collapse = ", "))
  expect_true(any(grepl("RoBMA", checked)), info = paste(checked, collapse = ", "))
  if (length(declined) > 0L) {
    message(
      "Predictor batch declined (generic route retained): ",
      paste(declined, collapse = ", ")
    )
  }
})
