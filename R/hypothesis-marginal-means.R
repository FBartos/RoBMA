#' @rdname hypothesis
#'
#' @param type for \code{marginal_means.brma} objects, whether to use
#' model-averaged (\code{"averaged"}) or conditional (\code{"conditional"})
#' marginal means.
#' @param parameter optional marginal-means parameter/term used to disambiguate
#' the hypothesis expression.
#'
#' @export
hypothesis.marginal_means.brma <- function(object, hypothesis,
                                           parameter = NULL,
                                           type = NULL,
                                           logBF = FALSE, BF01 = FALSE,
                                           seed = NULL,
                                           density_method = c(
                                             "KDE", "normal", "qCMDE", "IWMDE"
                                           ),
                                           columns = "default",
                                           ...) {

  BayesTools::check_char(hypothesis, "hypothesis", check_length = 0, allow_NA = FALSE)
  BayesTools::check_bool(logBF, "logBF")
  BayesTools::check_bool(BF01, "BF01")
  BayesTools::check_real(seed, "seed", check_length = 1, allow_NULL = TRUE, allow_NA = FALSE)
  BayesTools::check_char(columns, "columns", check_length = 0, allow_NA = FALSE)
  .warn_unused_dots(
    dots    = list(...),
    allowed = character(),
    caller  = "hypothesis.marginal_means()"
  )

  type <- .hypothesis_marginal_means_type(object = object, type = type)
  selected <- .hypothesis_marginal_means_select_parameter(
    object     = object,
    hypothesis = hypothesis,
    parameter  = parameter
  )
  parameter <- selected[["parameter"]]
  hypothesis <- .hypothesis_brma_rewrite(
    hypothesis = hypothesis,
    aliases    = selected[["aliases"]],
    parameter  = parameter
  )

  density_method <- .density_method_normalize(
    density_method = density_method,
    allow_normal   = TRUE
  )
  if (density_method %in% c("qCMDE", "IWMDE")) {
    .hypothesis_marginal_means_require_precomputed(
      object         = object,
      parameter      = parameter,
      hypothesis     = hypothesis,
      type           = type,
      density_method = density_method
    )
    density_method <- "precomputed"
  }

  inference <- object[["inference"]]
  inference[["conditional"]] <- inference[[type]]
  class(inference) <- unique(c(class(inference), "marginal_inference"))

  BayesTools::hypothesis_BF(
    posterior      = inference,
    hypothesis     = hypothesis,
    parameter      = parameter,
    logBF          = logBF,
    BF01           = BF01,
    seed           = seed,
    columns        = columns,
    density_method = density_method
  )
}


.hypothesis_marginal_means_type <- function(object, type) {

  if (is.null(type)) {
    return("averaged")
  }

  BayesTools::check_char(type, "type", check_length = 1)
  type <- match.arg(type, c("averaged", "conditional"))
  if (is.null(object[["inference"]][[type]])) {
    stop("No marginal means are available for type = '", type, "'.",
         call. = FALSE)
  }

  return(type)
}


.hypothesis_marginal_means_select_parameter <- function(object, hypothesis,
                                                        parameter) {

  aliases <- .hypothesis_marginal_means_aliases(object)
  if (!is.null(parameter)) {
    selected <- .marginal_means_select_parameter(object, parameter)[["parameter"]]
    aliases <- .hypothesis_brma_aliases_for_parameter(aliases, selected)
    return(list(parameter = selected, aliases = aliases))
  }

  roots <- .hypothesis_brma_symbol_roots(hypothesis)
  roots <- unique(roots[nzchar(roots)])
  if (length(roots) == 0L) {
    stop("Hypothesis must reference a marginal-means parameter.",
         call. = FALSE)
  }

  matches <- lapply(roots, function(root) aliases[[root]])
  known <- unlist(matches[!vapply(matches, is.null, logical(1))],
                  use.names = FALSE)
  known <- unique(known)

  if (length(known) == 1L) {
    aliases <- .hypothesis_brma_aliases_for_parameter(aliases, known)
    return(list(parameter = known, aliases = aliases))
  }

  if (length(known) > 1L) {
    stop(
      "Hypothesis references multiple marginal-means parameters (",
      paste(known, collapse = ", "),
      "). Specify 'parameter'.",
      call. = FALSE
    )
  }

  stop(
    "Could not infer a marginal-means parameter from the hypothesis. ",
    "Available quantities are: ", .hypothesis_brma_available_aliases(aliases),
    ".",
    call. = FALSE
  )
}


.hypothesis_marginal_means_aliases <- function(object) {

  term_map <- object[["term_map"]]
  aliases  <- list()

  for (i in seq_len(nrow(term_map))) {
    parameter <- term_map[["parameter"]][[i]]
    term      <- term_map[["term"]][[i]]

    aliases[[parameter]] <- parameter
    if (!is.null(term) && nzchar(term)) {
      aliases[[term]] <- parameter
      aliases[[.formula_design_display_names(term)]] <- parameter
    }
  }

  return(aliases)
}


.hypothesis_marginal_means_require_precomputed <- function(
    object, parameter, hypothesis, type, density_method) {

  point_refs <- .hypothesis_brma_point_refs(hypothesis, parameter)
  if (nrow(point_refs) == 0L) {
    return(invisible(TRUE))
  }

  object_density_method <- object[["density_method"]]
  if (is.null(object_density_method) ||
      !identical(
        .density_method_normalize(object_density_method, allow_normal = TRUE),
        .density_method_normalize(density_method)
      )) {
    stop(
      "The marginal-means object does not contain ", density_method,
      " precomputed ordinates. Recompute it with ",
      "'marginal_means(..., density_method = \"", density_method, "\")'.",
      call. = FALSE
    )
  }

  samples <- object[["inference"]][[type]][[parameter]]
  for (i in seq_len(nrow(point_refs))) {
    ref <- point_refs[i, , drop = FALSE]
    sample <- samples
    label <- parameter
    if (!is.na(ref[["level"]])) {
      label <- paste0(parameter, "[", ref[["level"]], "]")
      if (!is.list(samples) || !ref[["level"]] %in% names(samples)) {
        stop("Hypothesis references unknown level '", ref[["level"]],
             "' for parameter '", parameter, "'.", call. = FALSE)
      }
      sample <- samples[[ref[["level"]]]]
    }

    ordinate <- attr(sample, "posterior_ordinate", exact = TRUE)
    if (!.hypothesis_brma_ordinate_has_value(ordinate, ref[["value"]])) {
      stop(
        "Precomputed ", density_method,
        " posterior ordinate is unavailable for '", label, " = ",
        ref[["value"]], "'. Use type = 'conditional' when testing ",
        "marginal-means point nulls, recompute 'marginal_means()' with ",
        "the requested null_hypothesis, or use density_method = 'KDE'.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}


.hypothesis_brma_ordinate_has_value <- function(ordinate, value) {

  if (is.null(ordinate)) {
    return(FALSE)
  }
  if (is.list(ordinate) &&
      !is.null(ordinate[["status"]]) &&
      !identical(ordinate[["status"]], "ok")) {
    return(FALSE)
  }

  source <- ordinate
  if (is.list(ordinate)) {
    if (!is.null(ordinate[["ordinate"]]) &&
        (is.list(ordinate[["ordinate"]]) ||
         is.data.frame(ordinate[["ordinate"]]))) {
      source <- ordinate[["ordinate"]]
    }
    if (!is.null(ordinate[["ordinates"]]) &&
        (is.list(ordinate[["ordinates"]]) ||
         is.data.frame(ordinate[["ordinates"]]))) {
      source <- ordinate[["ordinates"]]
    }
  }

  if (is.list(source) &&
      !is.data.frame(source) &&
      is.null(source[["x"]]) &&
      is.null(source[["value"]]) &&
      is.null(source[["null_hypothesis"]])) {
    return(any(vapply(
      source,
      .hypothesis_brma_ordinate_has_value,
      logical(1),
      value = value
    )))
  }

  if (is.data.frame(source)) {
    x_name <- intersect(c("x", "value", "null_hypothesis"), colnames(source))[1L]
    y_name <- intersect(c("y", "ordinate", "height", "posterior_height"),
                        colnames(source))[1L]
    if (is.na(x_name) || is.na(y_name)) {
      return(FALSE)
    }
    x <- source[[x_name]]
    y <- source[[y_name]]
  } else if (is.list(source)) {
    x <- NULL
    y <- NULL
    for (name in c("x", "value", "null_hypothesis")) {
      if (!is.null(source[[name]])) {
        x <- source[[name]]
        break
      }
    }
    for (name in c("y", "ordinate", "height", "posterior_height")) {
      if (!is.null(source[[name]])) {
        y <- source[[name]]
        break
      }
    }
    if (is.null(x) || is.null(y)) {
      return(FALSE)
    }
  } else {
    return(FALSE)
  }

  x <- suppressWarnings(as.numeric(x))
  y <- suppressWarnings(as.numeric(y))
  if (length(x) != length(y)) {
    return(FALSE)
  }

  tolerance <- sqrt(.Machine$double.eps) * max(1, abs(value))
  any(is.finite(x) & is.finite(y) & y > 0 & abs(x - value) <= tolerance)
}

