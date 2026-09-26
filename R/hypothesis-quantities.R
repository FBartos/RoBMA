#' @title Available Hypothesis Quantities
#'
#' @description Lists parameter names and aliases accepted by
#' \code{hypothesis()} for a fitted object, with the tests that
#' \code{hypothesis()} runs on each of them.
#'
#' @param object a fitted \code{brma} or \code{marginal_means.brma} object.
#' @param ... unused.
#'
#' @details The eligibility columns render the plans that \code{hypothesis()}
#' executes (see the details of [hypothesis()]) for the statements
#' \code{<q> = <null>} (point tests), \code{<q> > <null>} (directional and
#' interval tests), and, for factor terms, \code{<level> = <other level>}
#' (level contrasts) of every quantity \code{q}. The null value is 0; when 0
#' is a prior point mass, a support boundary without a regular prior
#' ordinate, or outside the support (for example, a gated random-effect
#' standard deviation or a variance proportion), the point test is rendered
#' at an interior value of the support instead (the midpoint of a bounded
#' support, one unit inside a half-bounded support, and 1 otherwise): point
#' hypotheses at such values follow the per-value rules of
#' \code{hypothesis()}. A factor term lists the point-test methods available
#' for all of its levels that the contrast does not fix, and the
#' contrast-test methods available for all pairs of its levels. Region tests
#' do not use a density method. The normal approximation is a rough check
#' and is not listed.
#'
#' @return A data frame with the columns \code{alias} (a name that
#' \code{hypothesis()} resolves to the row's quantity with the row's
#' \code{component}; a name shared by several quantities of a component, such
#' as the \code{intercept} of several scale formulas, is not listed),
#' \code{parameter},
#' \code{component}, \code{term}, \code{bracket} (the level selector of factor
#' terms), \code{point_test}, \code{direction_test}, \code{contrast_test}
#' (\code{NA} for quantities without levels), \code{point_test_methods},
#' \code{contrast_test_methods}, and \code{reason} (the refusals of the
#' rendered statements).
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
  out  <- catalog[, c("alias", "parameter", "component", "term"), drop = FALSE]
  keys <- if ("quantity_id" %in% names(catalog)) {
    catalog[["quantity_id"]]
  } else {
    paste(catalog[["parameter"]], catalog[["component"]], sep = "\r")
  }
  unique_rows <- match(unique(keys), keys)

  if (is.null(object[["fit"]]) || length(object[["fit"]]) == 0L) {
    # hypothesis() requires a fitted object.
    rows <- lapply(unique_rows, function(i) {
      .hypothesis_quantities_row(
        reason = "'hypothesis' requires a fitted brma object."
      )
    })
  } else {
    metadata <- .brma_parameter_catalog_metadata(object)
    entries  <- metadata[["entries"]]
    cache    <- .hypothesis_plan_cache()
    rows <- lapply(unique_rows, function(i) {
      entry <- as.list(entries[
        match(catalog[["quantity_id"]][[i]], entries[["quantity_id"]]),
        setdiff(names(entries), "aliases"),
        drop = FALSE
      ])
      .hypothesis_quantities_render(
        object   = object,
        entry    = entry,
        metadata = metadata,
        cache    = cache
      )
    })
  }
  rows <- do.call(rbind, rows)[match(keys, unique(keys)), , drop = FALSE]
  out  <- cbind(out, rows)
  if (!is.null(object[["fit"]]) && length(object[["fit"]]) > 0L) {
    # Only aliases that name their own quantity within its component, as
    # hypothesis() resolves them; ambiguous aliases are not listed.
    resolves <- vapply(seq_len(nrow(catalog)), function(i) {
      .hypothesis_quantities_alias_resolves(
        catalog     = metadata[["catalog"]],
        alias       = catalog[["alias"]][[i]],
        component   = catalog[["component"]][[i]],
        quantity_id = catalog[["quantity_id"]][[i]]
      )
    }, logical(1))
    out <- out[resolves, , drop = FALSE]
  }
  out[["bracket"]] <- ifelse(
    out[["bracket"]],
    paste0(out[["parameter"]], "[level]"),
    NA_character_
  )
  out <- out[, c(
    "alias", "parameter", "component", "term", "bracket", "point_test",
    "direction_test", "contrast_test", "point_test_methods",
    "contrast_test_methods", "reason"
  ), drop = FALSE]
  rownames(out) <- NULL

  return(out)
}


# One row of eligibility columns.
.hypothesis_quantities_row <- function(bracket = FALSE, point_test = FALSE,
                                       direction_test = FALSE,
                                       contrast_test = NA,
                                       point_test_methods = "",
                                       contrast_test_methods = NA_character_,
                                       reason = "") {

  data.frame(
    bracket               = bracket,
    point_test            = point_test,
    direction_test        = direction_test,
    contrast_test         = contrast_test,
    point_test_methods    = point_test_methods,
    contrast_test_methods = contrast_test_methods,
    reason                = reason,
    stringsAsFactors      = FALSE
  )
}


# The rendered statements of a catalog entry and their plans. Statements
# name the quantity as .hypothesis_quantities_reference_root() does.
.hypothesis_quantities_plans <- function(object, entry, metadata, cache) {

  # Statements are planned as hypothesis() plans them.
  plan_of <- function(statement) {
    .hypothesis_plans(
      object     = object,
      hypothesis = statement,
      component  = entry[["component"]],
      metadata   = metadata,
      cache      = cache
    )[[1L]]
  }
  root <- .hypothesis_quantities_reference_root(metadata, entry)
  reference <- function(level = NULL) {
    if (is.null(level)) {
      return(paste0("`", root, "`"))
    }
    paste0("`", root, "[", level, "]`")
  }

  if (!identical(entry[["role"]], "formula_coefficient_group")) {
    point <- .hypothesis_quantities_point_plan(plan_of, reference())
    return(list(
      point    = list(point),
      region   = list(plan_of(paste0(
        reference(), " > ", .hypothesis_quantities_region_value(point)
      ))),
      contrast = NULL
    ))
  }

  levels <- .hypothesis_plan_term_levels(metadata, entry)
  tested <- levels[["level"]][!levels[["fixed"]]]
  point  <- lapply(tested, function(level) {
    .hypothesis_quantities_point_plan(plan_of, reference(level))
  })
  region <- lapply(seq_along(tested), function(i) {
    plan_of(paste0(
      reference(tested[[i]]), " > ",
      .hypothesis_quantities_region_value(point[[i]])
    ))
  })
  pairs <- if (nrow(levels) > 1L) utils::combn(levels[["level"]], 2L) else NULL
  contrast <- lapply(seq_len(NCOL(pairs)), function(i) {
    plan_of(paste0(reference(pairs[1L, i]), " = ", reference(pairs[2L, i])))
  })

  list(point = point, region = region, contrast = contrast)
}


# The name of a catalog entry's quantity in its rendered statements.
# Random-effect quantities and factor terms are named by their parameter (the
# level references of factor terms are resolved across components). Other
# quantities are named by their displayed alias within the entry's component,
# unless that alias also names another quantity of the component (e.g. the
# 'intercept' of several scale formulas): such quantities are named by their
# catalog selector.
.hypothesis_quantities_reference_root <- function(metadata, entry) {

  parameter <- entry[["parameter"]]
  if (identical(entry[["component"]], "random") ||
      identical(entry[["role"]], "formula_coefficient_group")) {
    return(parameter)
  }

  entries <- metadata[["entries"]]
  aliases <- entries[["aliases"]][[
    match(entry[["quantity_id"]], entries[["quantity_id"]])
  ]]
  aliases <- as.list(stats::setNames(rep(parameter, length(aliases)), aliases))
  label   <- .hypothesis_brma_alias_label(aliases, parameter)
  if (.hypothesis_quantities_alias_resolves(
    catalog     = metadata[["catalog"]],
    alias       = label,
    component   = entry[["component"]],
    quantity_id = entry[["quantity_id"]]
  )) {
    return(label)
  }

  quantities <- metadata[["catalog"]][["quantities"]]
  BayesTools::parameter_labels(
    quantities[
      match(entry[["quantity_id"]], quantities[["quantity_id"]]), ,
      drop = FALSE
    ],
    style = "selector"
  )
}


# Whether 'alias' names the quantity 'quantity_id' within 'component': the
# resolution hypothesis() applies to the statements of that component. An
# alias that also names another quantity of the component is ambiguous there
# (e.g. the 'intercept' of several scale formulas).
.hypothesis_quantities_alias_resolves <- function(catalog, alias, component,
                                                  quantity_id) {

  tryCatch(
    identical(
      BayesTools::parameter_catalog_resolve(
        catalog        = catalog,
        alias          = alias,
        component      = component,
        simplify_names = TRUE
      )[["quantity_id"]],
      quantity_id
    ),
    BayesTools_parameter_ambiguous = function(error) FALSE
  )
}


# The plan of '<q> = 0', or of '<q> = <interior value>' when the value 0
# itself is refused (a prior point mass, a support boundary without a
# regular prior ordinate, or a value outside the support): the statement
# rendered for a quantity's point tests.
.hypothesis_quantities_point_plan <- function(plan_of, reference) {

  plan    <- plan_of(paste0(reference, " = 0"))
  refusal <- .hypothesis_plan_status(plan, "KDE")
  if (!is.null(refusal) && isTRUE(refusal[["value_specific"]])) {
    value <- .hypothesis_quantities_interior_value(plan[["support"]])
    plan  <- plan_of(paste0(reference, " = ", format(value, digits = 15L)))
    plan[["null_value"]] <- value
    return(plan)
  }
  plan[["null_value"]] <- 0

  plan
}


# The value of the rendered region statement '<q> > <value>': the null value
# of the point statement, or an interior value when that null value is a
# support bound (a region beyond a bound has no complementary prior mass).
.hypothesis_quantities_region_value <- function(point) {

  value   <- point[["null_value"]]
  support <- point[["support"]]
  if (!is.null(support) && (value <= support[[1L]] || value >= support[[2L]])) {
    value <- .hypothesis_quantities_interior_value(support)
  }

  format(value, digits = 15L)
}


# A value inside the support bounds 'support' (c(lower, upper), or NULL for
# an unbounded support): the midpoint of a bounded support, one unit inside a
# half-bounded support, and 1 otherwise.
.hypothesis_quantities_interior_value <- function(support) {

  bounds <- if (is.null(support)) c(-Inf, Inf) else as.numeric(support)
  if (all(is.finite(bounds))) {
    return(mean(bounds))
  }
  if (is.finite(bounds[[1L]])) {
    return(bounds[[1L]] + 1)
  }
  if (is.finite(bounds[[2L]])) {
    return(bounds[[2L]] - 1)
  }

  1
}


.hypothesis_quantities_render <- function(object, entry, metadata, cache) {

  plans <- .hypothesis_quantities_plans(object, entry, metadata, cache)
  .hypothesis_quantities_render_plans(
    plans   = plans,
    bracket = identical(entry[["role"]], "formula_coefficient_group")
  )
}


# The eligibility columns of the plans of one quantity: the methods that
# evaluate every point (contrast) plan, whether every region plan runs, and
# the distinct refusals of the plans.
.hypothesis_quantities_render_plans <- function(plans, bracket) {

  advertised <- .hypothesis_plan_advertised_methods()
  methods_of <- function(plans) {
    if (length(plans) == 0L) {
      return(character())
    }
    advertised[vapply(advertised, function(method) {
      all(vapply(plans, function(plan) {
        is.null(.hypothesis_plan_status(plan, method))
      }, logical(1)))
    }, logical(1))]
  }
  point_methods    <- methods_of(plans[["point"]])
  contrast_methods <- if (is.null(plans[["contrast"]])) {
    NULL
  } else {
    methods_of(plans[["contrast"]])
  }
  direction_test <- length(plans[["region"]]) > 0L && all(vapply(
    plans[["region"]],
    function(plan) is.null(.hypothesis_plan_status(plan, "KDE")),
    logical(1)
  ))
  reasons <- unlist(lapply(
    c(plans[["point"]], plans[["region"]], plans[["contrast"]]),
    function(plan) {
      vapply(advertised, function(method) {
        refusal <- .hypothesis_plan_status(plan, method)
        if (is.null(refusal)) NA_character_ else refusal[["reason"]]
      }, character(1))
    }
  ), use.names = FALSE)
  reasons <- unique(reasons[!is.na(reasons)])

  .hypothesis_quantities_row(
    bracket               = bracket,
    point_test            = length(point_methods) > 0L,
    direction_test        = direction_test,
    contrast_test         = if (is.null(contrast_methods)) {
      NA
    } else {
      length(contrast_methods) > 0L
    },
    point_test_methods    = paste(point_methods, collapse = ", "),
    contrast_test_methods = if (is.null(contrast_methods)) {
      NA_character_
    } else {
      paste(contrast_methods, collapse = ", ")
    },
    reason                = paste(reasons, collapse = " ")
  )
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
  averaged <- object[["inference"]][["averaged"]]
  if (is.null(averaged)) {
    averaged <- object[["inference"]][["conditional"]]
  }
  rows  <- list()
  cache <- .hypothesis_plan_cache()
  for (i in seq_len(nrow(term_map))) {

    parameter <- term_map[["parameter"]][[i]]
    levels    <- if (is.list(averaged[[parameter]])) {
      level_names <- names(averaged[[parameter]])
      if (is.null(level_names) ||
          length(level_names) != length(averaged[[parameter]]) ||
          any(!nzchar(level_names))) {
        stop(
          "Grouped marginal means must have non-empty level names.",
          call. = FALSE
        )
      }
      level_names
    } else {
      list(NULL)
    }
    for (level in levels) {
      rows[[length(rows) + 1L]] <- cbind(
        .hypothesis_quantities_marginal_row(term_map[i, , drop = FALSE], level),
        .hypothesis_quantities_marginal_render(
          object    = object,
          parameter = parameter,
          level     = level,
          levels    = if (is.null(level)) NULL else unlist(levels),
          cache     = cache
        )
      )
    }
  }
  out <- do.call(rbind, rows)
  out <- out[, c(
    "alias", "parameter", "component", "term", "bracket", "point_test",
    "direction_test", "contrast_test", "point_test_methods",
    "contrast_test_methods", "reason"
  ), drop = FALSE]
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


# The plans of one marginal mean (a level of a factor term, or a scalar
# term): its point and region statements and its contrasts with the other
# levels.
.hypothesis_quantities_marginal_render <- function(object, parameter, level,
                                                   levels, cache) {

  plan_of <- function(statement) {
    .hypothesis_plan_marginal_means(
      object    = object,
      statement = statement,
      parameter = parameter,
      cache     = cache
    )
  }
  reference <- function(level) {
    if (is.null(level)) {
      return(paste0("`", parameter, "`"))
    }
    paste0("`", parameter, "[", level, "]`")
  }
  point  <- .hypothesis_quantities_point_plan(plan_of, reference(level))
  plans  <- list(
    point    = list(point),
    region   = list(plan_of(paste0(
      reference(level), " > ", .hypothesis_quantities_region_value(point)
    ))),
    contrast = if (!is.null(level) && length(levels) > 1L) {
      lapply(setdiff(levels, level), function(other) {
        plan_of(paste0(reference(level), " = ", reference(other)))
      })
    }
  )
  out <- .hypothesis_quantities_render_plans(plans, bracket = FALSE)
  out[["bracket"]] <- NULL

  out
}
