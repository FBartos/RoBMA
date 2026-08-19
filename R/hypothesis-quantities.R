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
    likelihood_aware <- .hypothesis_quantities_iwmde_capability(object)
    random_parameters <- unique(out[["parameter"]][is_random])
    random_routes <- lapply(random_parameters, function(parameter) {

      i <- match(parameter, specs[["parameter"]])
      spec         <- as.list(specs[i, , drop = FALSE])
      source_prior <- .brma_random_parameter_source_prior(object, spec)
      direct_reason <- .brma_random_parameter_point_test_reason(
        spec         = spec,
        prior        = bundle[["priors"]][[parameter]],
        source_prior = source_prior
      )
      likelihood_reason <- .brma_random_parameter_point_test_reason(
        spec         = spec,
        prior        = bundle[["priors"]][[parameter]],
        source_prior = source_prior,
        derived      = TRUE
      )
      point_methods <- .hypothesis_quantities_random_point_methods(
        parameter          = parameter,
        object             = object,
        methods            = likelihood_aware[["methods"]],
        direct_allowed     = !nzchar(direct_reason),
        likelihood_allowed = !nzchar(likelihood_reason)
      )
      list(
        point_methods = point_methods,
        reason        = if (nzchar(point_methods)) "" else direct_reason
      )
    })
    names(random_routes) <- random_parameters
    out[["bracket"]][is_random]                <- NA_character_
    point_methods <- vapply(
      random_routes[out[["parameter"]][is_random]],
      `[[`,
      character(1),
      "point_methods"
    )
    out[["point_test"]][is_random]             <- nzchar(point_methods)
    out[["point_test_methods"]][is_random]     <- point_methods
    out[["direction_test_methods"]][is_random] <- "KDE, normal"
    out[["reason"]][is_random] <- vapply(
      random_routes[out[["parameter"]][is_random]],
      `[[`,
      character(1),
      "reason"
    )

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


.hypothesis_quantities_random_point_methods <- function(
    parameter, object, methods, direct_allowed = TRUE,
    likelihood_allowed = TRUE) {

  target <- tryCatch(
    .brma_random_parameter_density_target(object, parameter),
    error = function(error) NULL
  )
  supported <- is.list(target) &&
    is.character(target[["parameter"]]) &&
    length(target[["parameter"]]) == 1L &&
    !is.na(target[["parameter"]]) && nzchar(target[["parameter"]])

  paste(c(
    if (isTRUE(unname(direct_allowed))) "KDE" else character(),
    if (supported && isTRUE(unname(likelihood_allowed))) {
      methods
    } else {
      character()
    }
  ), collapse = ", ")
}


.hypothesis_quantities_brma_route <- function(object, row, entry = NULL) {

  likelihood_aware <- .hypothesis_quantities_iwmde_capability(object)
  point_methods     <- c("KDE", likelihood_aware[["methods"]])
  out <- list(
    bracket               = !is.null(entry) &&
      identical(entry[["role"]], "formula_coefficient_group"),
    point_test             = TRUE,
    direction_test         = TRUE,
    point_test_methods     = paste(point_methods, collapse = ", "),
    direction_test_methods = "KDE, normal",
    reason                 = likelihood_aware[["reason"]]
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
    object   = object,
    selected = selected
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
    out[["reason"]] <- paste0(
      "KDE requires runtime structural certification of an atom-free, ",
      "unconditional scalar target."
    )
    return(out)
  }

  out[["point_test"]]             <- FALSE
  out[["direction_test"]]         <- FALSE
  out[["point_test_methods"]]     <- ""
  out[["direction_test_methods"]] <- ""
  out[["reason"]]                 <- route[["reason"]]
  return(out)
}


.hypothesis_quantities_iwmde_capability <- function(object) {

  methods      <- c("qCMDE", "IWMDE")
  capabilities <- lapply(methods, function(density_method) {

    .iwmde_capability(
      object         = object,
      density_method = density_method
    )
  })
  available <- vapply(capabilities, `[[`, logical(1), "available")
  reasons   <- unique(vapply(
    capabilities[!available],
    `[[`,
    character(1),
    "reason"
  ))

  return(list(
    methods = methods[available],
    reason  = paste(reasons, collapse = " ")
  ))
}


#' @rdname hypothesis_quantities
#' @export
hypothesis_quantities.marginal_means.brma <- function(object, ...) {

  .warn_unused_dots(
    dots    = list(...),
    allowed = character(),
    caller  = "hypothesis_quantities()"
  )
  term_map    <- object[["term_map"]]
  conditional <- object[["inference"]][["conditional"]]
  rows        <- list()
  samples     <- list()
  row_i       <- 0L
  for (i in seq_len(nrow(term_map))) {

    parameter_samples <- conditional[[term_map[["parameter"]][[i]]]]
    if (is.list(parameter_samples)) {
      level_names <- names(parameter_samples)
      if (is.null(level_names) || length(level_names) != length(parameter_samples) ||
          any(!nzchar(level_names))) {
        stop(
          "Grouped marginal means must have non-empty level names.",
          call. = FALSE
        )
      }
      for (level_i in seq_along(parameter_samples)) {
        row_i <- row_i + 1L
        rows[[row_i]] <- .hypothesis_quantities_marginal_row(
          term_map[i, , drop = FALSE],
          level = level_names[[level_i]]
        )
        samples[[row_i]] <- parameter_samples[[level_i]]
      }
    } else {
      row_i <- row_i + 1L
      rows[[row_i]] <- .hypothesis_quantities_marginal_row(
        term_map[i, , drop = FALSE]
      )
      samples[[row_i]] <- parameter_samples
    }
  }
  out <- do.call(rbind, rows)
  likelihood_aware <- .hypothesis_quantities_iwmde_capability(
    object[["source_object"]]
  )
  out <- .hypothesis_quantities_add_eligibility(
    out,
    point_test_methods = paste(
      c("KDE", likelihood_aware[["methods"]]),
      collapse = ", "
    )
  )
  out[["reason"]] <- likelihood_aware[["reason"]]
  out[["direction_test_methods"]] <- "KDE"
  for (i in seq_len(nrow(out))) {
    values <- as.numeric(samples[[i]])
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
  rownames(out) <- NULL
  return(out)
}


.hypothesis_quantities_marginal_row <- function(term_row, level = NULL) {

  parameter <- term_row[["parameter"]][[1L]]
  return(data.frame(
    alias      = term_row[["term"]][[1L]],
    parameter  = parameter,
    component  = "marginal_means",
    term       = term_row[["term"]][[1L]],
    bracket    = if (is.null(level)) {
      NA_character_
    } else {
      paste0(parameter, "[", level, "]")
    },
    stringsAsFactors = FALSE,
    check.names = FALSE
  ))
}


.hypothesis_quantities_add_eligibility <- function(out, point_test_methods) {

  out[["point_test"]]             <- TRUE
  out[["direction_test"]]         <- TRUE
  out[["point_test_methods"]]     <- point_test_methods
  out[["direction_test_methods"]] <- "KDE, qCMDE, IWMDE"
  out[["reason"]]                 <- ""

  return(out)
}
