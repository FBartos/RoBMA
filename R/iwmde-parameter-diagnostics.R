# ============================================================================ #
# IWMDE Parameter Diagnostic Orchestration
# ============================================================================ #

.iwmde_normalize_display_grid <- function(display_grid) {

  return(match.arg(display_grid, c("adaptive", "uniform")))
}


.iwmde_relabel_diagnostic <- function(diagnostic, parameter) {

  diagnostic[["parameter"]] <- parameter

  return(diagnostic)
}
