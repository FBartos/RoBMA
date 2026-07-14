# ============================================================================ #
# IWMDE Parameter Diagnostic Orchestration
# ============================================================================ #

.iwmde_parameter_ordinate_diagnostic <- function(context, parameter, values,
                                                 max_samples,
                                                 normalization_points,
                                                 normalization_prob,
                                                 method = c("q_grid_cmde", "iwmde"),
                                                 parameter_spec = NULL,
                                                 diagnostic_cache = NULL) {

  method <- .iwmde_normalize_method(method)
  plan <- .iwmde_plan(
    context         = context,
    parameter       = parameter,
    density_method  = if (identical(method, "q_grid_cmde")) "qCMDE" else "IWMDE",
    density_control = list(
      n_points             = max(20L, length(values)),
      max_samples          = max_samples,
      normalization_points = normalization_points,
      normalization_prob   = normalization_prob,
      display_grid         = "ordinate"
    ),
    outputs         = "ordinate",
    values          = values,
    parameter_spec  = parameter_spec,
    metadata        = NULL
  )

  return(.iwmde_execute_plan_diagnostic(
    context          = context,
    plan             = plan,
    output           = "ordinate",
    execution_cache  = new.env(parent = emptyenv()),
    diagnostic_cache = diagnostic_cache
  ))
}


.iwmde_ordinate_interior_values <- function(values, support, xlim) {

  out   <- values
  width <- diff(xlim)
  if (!is.finite(width) || width <= 0) {
    width <- 1
  }
  eps <- sqrt(.Machine$double.eps) * max(1, width)

  if (is.finite(support[1])) {
    out[out <= support[1]] <- support[1] + eps
  }
  if (is.finite(support[2])) {
    out[out >= support[2]] <- support[2] - eps
  }

  return(out)
}


.iwmde_normalize_display_grid <- function(display_grid) {

  return(match.arg(display_grid, c("adaptive", "uniform")))
}


.iwmde_parameter_diagnostic <- function(context, parameter, n_points,
                                        max_samples, normalization_points,
                                        normalization_prob,
                                        method = c("q_grid_cmde", "iwmde"),
                                        display_grid_method = "adaptive",
                                        include_values = NULL,
                                        parameter_spec = NULL,
                                        diagnostic_cache = NULL) {

  method <- .iwmde_normalize_method(method)
  plan <- .iwmde_plan(
    context         = context,
    parameter       = parameter,
    density_method  = if (identical(method, "q_grid_cmde")) "qCMDE" else "IWMDE",
    density_control = list(
      n_points             = n_points,
      max_samples          = max_samples,
      normalization_points = normalization_points,
      normalization_prob   = normalization_prob,
      display_grid         = display_grid_method
    ),
    outputs         = "density",
    values          = include_values,
    parameter_spec  = parameter_spec,
    metadata        = NULL
  )

  return(.iwmde_execute_plan_diagnostic(
    context          = context,
    plan             = plan,
    output           = "density",
    execution_cache  = new.env(parent = emptyenv()),
    diagnostic_cache = diagnostic_cache
  ))
}


.iwmde_relabel_diagnostic <- function(diagnostic, parameter) {

  diagnostic[["parameter"]] <- parameter

  return(diagnostic)
}
