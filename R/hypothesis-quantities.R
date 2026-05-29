#' @title Available Hypothesis Quantities
#'
#' @description Lists parameter names and aliases accepted by
#' \code{hypothesis()} for a fitted object.
#'
#' @param object a fitted \code{brma} or \code{marginal_means.brma} object.
#' @param ... unused.
#'
#' @return A data frame of available quantities, components, and aliases.
#'
#' @export
hypothesis_quantities <- function(object, ...) {

  UseMethod("hypothesis_quantities")
}


#' @rdname hypothesis_quantities
#' @export
hypothesis_quantities.brma <- function(object, ...) {

  .warn_unused_dots(
    dots    = list(...),
    allowed = character(),
    caller  = "hypothesis_quantities()"
  )
  catalog <- .brma_parameter_catalog(object)
  out <- catalog[, c("alias", "parameter", "component", "term"), drop = FALSE]
  out[["bracket"]] <- paste0(out[["parameter"]], "[level]")
  rownames(out) <- NULL
  return(out)
}


#' @rdname hypothesis_quantities
#' @export
hypothesis_quantities.marginal_means.brma <- function(object, ...) {

  .warn_unused_dots(
    dots    = list(...),
    allowed = character(),
    caller  = "hypothesis_quantities()"
  )
  term_map <- object[["term_map"]]
  out <- data.frame(
    alias      = term_map[["term"]],
    parameter  = term_map[["parameter"]],
    component  = "marginal_means",
    term       = term_map[["term"]],
    bracket    = paste0(term_map[["parameter"]], "[level]"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(out) <- NULL
  return(out)
}
