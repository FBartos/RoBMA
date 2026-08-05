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
  metadata <- if (!is.null(object[["fit"]])) {
    .brma_parameter_catalog_metadata(object)[["entries"]]
  } else {
    NULL
  }
  keys <- if ("quantity_id" %in% names(catalog)) {
    catalog[["quantity_id"]]
  } else {
    paste(catalog[["parameter"]], catalog[["component"]], sep = "\r")
  }
  unique_rows <- match(unique(keys), keys)
  routes <- lapply(unique_rows, function(i) {

    entry <- if (is.null(metadata) ||
                 !"quantity_id" %in% names(catalog)) {
      NULL
    } else {
      entry_i <- match(catalog[["quantity_id"]][[i]], metadata[["quantity_id"]])
      if (is.na(entry_i)) NULL else as.list(metadata[entry_i, , drop = FALSE])
    }
    .hypothesis_quantities_brma_route(
      object = object,
      row    = catalog[i, , drop = FALSE],
      entry  = entry
    )
  })
  routes <- routes[match(keys, unique(keys))]
  out[["bracket"]] <- vapply(seq_along(routes), function(i) {

    if (routes[[i]][["bracket"]]) {
      paste0(out[["parameter"]][[i]], "[level]")
    } else {
      NA_character_
    }
  }, character(1))
  out[["point_test"]] <- vapply(routes, `[[`, logical(1), "point_test")
  out[["direction_test"]] <- vapply(
    routes, `[[`, logical(1), "direction_test"
  )
  out[["point_test_methods"]] <- vapply(
    routes, `[[`, character(1), "point_test_methods"
  )
  out[["direction_test_methods"]] <- vapply(
    routes, `[[`, character(1), "direction_test_methods"
  )
  out[["reason"]] <- vapply(routes, `[[`, character(1), "reason")
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


.hypothesis_quantities_brma_route <- function(object, row, entry = NULL) {

  point_methods <- c("KDE", "qCMDE", "IWMDE")
  if (inherits(object, "brma.glmm")) {
    point_methods <- setdiff(point_methods, "IWMDE")
  }
  out <- list(
    bracket               = !is.null(entry) &&
      identical(entry[["role"]], "formula_coefficient_group"),
    point_test             = TRUE,
    direction_test         = TRUE,
    point_test_methods     = paste(point_methods, collapse = ", "),
    direction_test_methods = "KDE, normal",
    reason                 = ""
  )

  fixed_value <- if (is.null(entry)) NULL else entry[["fixed_value"]]
  fixed <- !is.null(entry) && (
    identical(entry[["status"]], "fixed") ||
      (length(fixed_value) == 1L && is.finite(fixed_value))
  )
  if (fixed) {
    out[["point_test"]]             <- FALSE
    out[["direction_test"]]         <- FALSE
    out[["point_test_methods"]]     <- ""
    out[["direction_test_methods"]] <- ""
    out[["reason"]] <- paste0(
      "The quantity is fixed by the fitted model; posterior hypothesis ",
      "tests are undefined."
    )
    return(out)
  }
  formula_parameter <- if (is.null(entry)) NULL else entry[["formula_parameter"]]
  if (is.null(entry) || identical(row[["component"]], "random") ||
      identical(entry[["role"]], "formula_coefficient_group") ||
      length(formula_parameter) != 1L || is.na(formula_parameter) ||
      !nzchar(formula_parameter)) {
    return(out)
  }

  selected <- list(
    parameter = row[["parameter"]],
    component = row[["component"]],
    entry     = entry
  )
  target <- .hypothesis_brma_formula_coefficient_target(
    object                    = object,
    selected                  = selected,
    standardized_coefficients = FALSE
  )
  if (is.null(target)) {
    return(out)
  }
  route <- .hypothesis_brma_formula_transform_route(target)
  if (route[["type"]] %in% c("identity", "affine")) {
    return(out)
  }
  if (identical(route[["type"]], "exp_affine")) {
    out[["point_test_methods"]]     <- "KDE"
    out[["direction_test_methods"]] <- "KDE"
    return(out)
  }

  out[["point_test"]]             <- FALSE
  out[["direction_test"]]         <- FALSE
  out[["point_test_methods"]]     <- ""
  out[["direction_test_methods"]] <- ""
  out[["reason"]]                 <- route[["reason"]]
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
    bracket    = NA_character_,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out <- .hypothesis_quantities_add_eligibility(
    out,
    point_test_methods = if (inherits(object[["source_object"]], "brma.glmm")) {
      "KDE, qCMDE"
    } else {
      "KDE, qCMDE, IWMDE"
    }
  )
  out[["direction_test_methods"]] <- "KDE"
  for (i in seq_len(nrow(out))) {
    samples <- object[["inference"]][["conditional"]][[out[["parameter"]][[i]]]]
    if (is.list(samples)) {
      out[["bracket"]][[i]] <- paste0(out[["parameter"]][[i]], "[level]")
    } else {
      values <- as.numeric(samples)
      values <- values[is.finite(values)]
      fixed <- length(values) > 0L && all(values == values[[1L]])
      if (fixed) {
        out[["point_test"]][[i]]             <- FALSE
        out[["direction_test"]][[i]]         <- FALSE
        out[["point_test_methods"]][[i]]     <- ""
        out[["direction_test_methods"]][[i]] <- ""
        out[["reason"]][[i]] <- paste0(
          "The quantity is fixed by the fitted model; posterior hypothesis ",
          "tests are undefined."
        )
      }
    }
  }
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
