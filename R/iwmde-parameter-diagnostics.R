# ============================================================================ #
# IWMDE Parameter Diagnostic Orchestration
# ============================================================================ #

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


.iwmde_relabel_diagnostic <- function(diagnostic, parameter) {

  diagnostic[["parameter"]] <- parameter

  return(diagnostic)
}
