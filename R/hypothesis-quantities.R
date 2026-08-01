#' @title Available Hypothesis Quantities
#'
#' @description Lists parameter names and aliases accepted by
#' \code{hypothesis()} for a fitted object.
#'
#' @param object a fitted \code{brma} or \code{marginal_means.brma} object.
#' @param ... unused.
#'
#' @return A data frame of available quantities, components, aliases, and
#' method-level test eligibility notes.
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
  catalog <- catalog[catalog[["component"]] != "bias", , drop = FALSE]
  out <- catalog[, c("alias", "parameter", "component", "term"), drop = FALSE]
  out[["bracket"]] <- paste0(out[["parameter"]], "[level]")
  out <- .hypothesis_quantities_add_eligibility(
    out,
    point_test_methods = "KDE, qCMDE, IWMDE"
  )
  is_random <- out[["component"]] == "random"
  if (any(is_random)) {
    bundle <- .brma_random_parameter_bundle(object)
    specs  <- bundle[["specs"]]
    index <- match(out[["parameter"]][is_random], specs[["parameter"]])
    reasons <- vapply(index, function(i) {
      .brma_random_parameter_point_test_reason(
        as.list(specs[i, , drop = FALSE])
      )
    }, character(1))
    out[["bracket"]][is_random]                <- NA_character_
    out[["point_test"]][is_random]             <- !nzchar(reasons)
    out[["point_test_methods"]][is_random]     <- ifelse(
      nzchar(reasons), "", "KDE, normal"
    )
    out[["direction_test_methods"]][is_random] <- "KDE, normal"
    out[["reason"]][is_random]                 <- reasons

    fixed <- vapply(index, function(i) {
      values <- bundle[["samples"]][, i]
      finite <- values[is.finite(values)]
      length(finite) > 0L && all(finite == finite[[1L]])
    }, logical(1))
    random_rows <- which(is_random)
    if (any(fixed)) {
      fixed_rows <- random_rows[fixed]
      out[["point_test"]][fixed_rows]             <- FALSE
      out[["direction_test"]][fixed_rows]         <- FALSE
      out[["point_test_methods"]][fixed_rows]     <- ""
      out[["direction_test_methods"]][fixed_rows] <- ""
      out[["reason"]][fixed_rows] <- paste0(
        ifelse(nzchar(out[["reason"]][fixed_rows]),
               paste0(out[["reason"]][fixed_rows], " "), ""),
        "The quantity is fixed by the fitted model; posterior hypothesis tests are undefined."
      )
    }
  }
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
  out <- .hypothesis_quantities_add_eligibility(
    out,
    point_test_methods = "KDE, qCMDE, IWMDE"
  )
  rownames(out) <- NULL
  return(out)
}


.hypothesis_quantities_add_eligibility <- function(out, point_test_methods) {

  out[["point_test"]]             <- TRUE
  out[["direction_test"]]         <- TRUE
  out[["point_test_methods"]]     <- point_test_methods
  out[["direction_test_methods"]] <- "KDE, qCMDE, IWMDE"
  out[["reason"]]                 <- ""

  return(out)
}
