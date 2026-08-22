#' @rdname hypothesis
#'
#' @param parameter optional marginal-means parameter/term used to disambiguate
#' the hypothesis expression.
#' @details Marginal-means hypotheses are specified on the fitted
#' linear-predictor scale. Display transformations stored by
#' \code{marginal_means()} do not transform hypothesis constants. Single-model
#' hypotheses and model-averaged region hypotheses use averaged marginal draws.
#' Model-averaged point-null hypotheses use alternative-conditioned draws to
#' avoid null atoms in averaged marginals. For model-averaged objects,
#' hypotheses that mix point and region events, and point hypotheses spanning
#' multiple factor levels, are not supported. Single-model hypotheses use one
#' common averaged posterior and may combine these events.
#'
#' @export
hypothesis.marginal_means.brma <- function(object, hypothesis,
                                           parameter = NULL,
                                           logBF = FALSE, BF01 = FALSE,
                                           seed = NULL,
                                           density_method = NULL,
                                           density_control = NULL,
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

  hypothesis <- BayesTools::hypothesis_parse(hypothesis)
  selected <- .hypothesis_marginal_means_select_parameter(
    object     = object,
    hypothesis = hypothesis,
    parameter  = parameter
  )
  parameter <- selected[["parameter"]]
  parameter_label <- .hypothesis_brma_alias_label(
    aliases   = selected[["aliases"]],
    parameter = parameter
  )
  hypothesis <- .hypothesis_brma_rewrite(
    hypothesis = hypothesis,
    aliases    = selected[["aliases"]],
    parameter  = parameter
  )
  route <- .hypothesis_marginal_means_route(
    object     = object,
    hypothesis = hypothesis,
    parameter  = parameter
  )
  level_contrast <- .hypothesis_brma_level_contrast_candidate(
    hypothesis = hypothesis,
    parameter  = parameter
  )

  density_method           <- .marginal_means_density_method(
    object,
    density_method
  )
  requested_density_method <- density_method
  estimated_ordinate        <- density_method %in% c("qCMDE", "IWMDE") &&
    (level_contrast || nrow(.hypothesis_brma_point_refs(
      hypothesis     = hypothesis,
      parameter      = parameter,
      require_direct = !level_contrast
    )) > 0L)
  precomputed_density       <- estimated_ordinate
  if (precomputed_density) {
    density_control <- .hypothesis_marginal_means_density_control(
      object          = object,
      density_method  = density_method,
      density_control = density_control
    )
    if (is.null(density_control[["normalization_points"]])) {
      density_control[["normalization_points"]] <- max(
        50L,
        density_control[["n_points"]]
      )
    }
    if (!level_contrast) {
      object <- .hypothesis_marginal_means_attach_iwmde(
        object          = object,
        parameter       = parameter,
        parameter_label = parameter_label,
        hypothesis      = hypothesis,
        inference_type  = route[["inference_type"]],
        density_method  = density_method,
        density_control = density_control
      )
    }
    density_method <- "precomputed"
  } else if (!is.null(density_control)) {
    if (density_method %in% c("qCMDE", "IWMDE")) {
      .density_control_normalize(
        density_method  = density_method,
        density_control = density_control,
        purpose         = "ordinate"
      )
    } else {
      stop("'density_control' is only used when 'density_method' is ",
           "'qCMDE' or 'IWMDE'.", call. = FALSE)
    }
  }

  if (!estimated_ordinate && density_method %in% c("qCMDE", "IWMDE")) {
    density_method <- "KDE"
  }

  inference <- object[["inference"]]
  if (is.null(inference[[route[["inference_type"]]]])) {
    stop(
      "'marginal_means' object does not contain ",
      route[["inference_type"]], " marginal means.",
      call. = FALSE
    )
  }
  if (level_contrast) {
    return(.hypothesis_brma_level_contrast_BF(
      object          = object[["source_object"]],
      posterior       = inference[[route[["inference_type"]]]][[parameter]],
      hypothesis      = hypothesis,
      parameter       = parameter,
      density_method  = requested_density_method,
      density_control = density_control,
      logBF           = logBF,
      BF01            = BF01,
      seed            = seed,
      columns         = columns
    ))
  }
  inference[["conditional"]] <- inference[[route[["inference_type"]]]]
  class(inference) <- unique(c(class(inference), "marginal_inference"))

  out <- BayesTools::hypothesis_BF(
    posterior      = inference,
    hypothesis     = hypothesis,
    parameter      = parameter,
    logBF          = logBF,
    BF01           = BF01,
    seed           = seed,
    columns        = columns,
    density_method = density_method
  )

  if (precomputed_density) {
    out <- .hypothesis_brma_append_iwmde_warnings(
      table     = out,
      posterior = inference[["conditional"]][[parameter]]
    )
  }

  return(out)
}


.hypothesis_marginal_means_route <- function(object, hypothesis, parameter) {

  model_averaged <- isTRUE(object[["model_averaged"]]) ||
    inherits(object[["source_object"]], "RoBMA")
  statement_routes <- vapply(
    hypothesis[["statements"]],
    function(statement) {

      types <- c(
        statement[["left"]][["type"]],
        statement[["right"]][["type"]]
      )
      if (all(types == "region")) {
        return("region")
      }
      if (all(types %in% c("point", "not_point"))) {
        return("point")
      }

      if (model_averaged) {
        stop(
          "Model-averaged marginal-means hypotheses cannot mix point and ",
          "region events because they require different posterior ",
          "conditioning. Use a pure point-null or a pure region hypothesis.",
          call. = FALSE
        )
      }
      return("mixed")
    },
    character(1)
  )
  if (model_averaged && length(unique(statement_routes)) != 1L) {
    stop(
      "A model-averaged marginal-means hypothesis request cannot mix ",
      "point-null and region statements because they require different ",
      "posterior conditioning.",
      call. = FALSE
    )
  }

  route <- if (length(unique(statement_routes)) == 1L) {
    statement_routes[[1L]]
  } else {
    "mixed"
  }
  if (model_averaged && identical(route, "point")) {
    levels <- unique(.hypothesis_marginal_means_ast_levels(
      hypothesis = hypothesis,
      parameter  = parameter
    ))
    if (length(levels) > 1L) {
      stop(
        "Point hypotheses spanning multiple marginal-means levels are not ",
        "supported because the levels need not share one conditioning event.",
        call. = FALSE
      )
    }
  }

  inference_type <- if (model_averaged && identical(route, "point")) {
    "conditional"
  } else {
    "averaged"
  }

  return(list(route = route, inference_type = inference_type))
}


.hypothesis_marginal_means_ast_levels <- function(hypothesis, parameter) {

  find_levels <- function(node) {

    if (!is.list(node)) {
      return(character())
    }
    if (identical(node[["type"]], "level_reference") &&
        identical(node[["parameter"]], parameter)) {
      return(as.character(node[["level"]]))
    }

    return(unlist(lapply(node, find_levels), use.names = FALSE))
  }

  return(find_levels(hypothesis[["statements"]]))
}


.hypothesis_marginal_means_select_parameter <- function(object, hypothesis,
                                                         parameter) {

  alias_catalog <- .hypothesis_marginal_means_alias_catalog(object)
  if (!is.null(parameter)) {
    selected <- .marginal_means_select_parameter(object, parameter)[["parameter"]]
    aliases <- .hypothesis_marginal_means_aliases(alias_catalog, selected)
    return(list(parameter = selected, aliases = aliases))
  }

  roots <- .hypothesis_brma_symbol_roots(hypothesis)
  roots <- unique(roots[nzchar(roots)])
  if (length(roots) == 0L) {
    stop("Hypothesis must reference a marginal-means parameter.",
         call. = FALSE)
  }

  matches <- lapply(roots, function(root) {

    unique(alias_catalog[["parameter"]][alias_catalog[["alias"]] == root])
  })
  ambiguous <- lengths(matches) > 1L
  if (any(ambiguous)) {
    root <- roots[ambiguous][[1L]]
    stop(
      "Marginal-means alias '", root, "' is ambiguous (",
      paste(matches[[which(ambiguous)[[1L]]]], collapse = ", "),
      "). Specify 'parameter' using an internal parameter name.",
      call. = FALSE
    )
  }
  known <- unlist(matches[lengths(matches) == 1L], use.names = FALSE)
  known <- unique(known)

  if (length(known) == 1L) {
    aliases <- .hypothesis_marginal_means_aliases(alias_catalog, known)
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
    "Available quantities are: ",
    paste0("'", sort(unique(alias_catalog[["alias"]])), "'", collapse = ", "),
    ".",
    call. = FALSE
  )
}


.hypothesis_marginal_means_alias_catalog <- function(object) {

  term_map <- object[["term_map"]]
  rows     <- list()

  for (i in seq_len(nrow(term_map))) {
    parameter <- term_map[["parameter"]][[i]]
    term      <- term_map[["term"]][[i]]

    aliases <- parameter
    if (!is.null(term) && nzchar(term)) {
      aliases <- c(aliases, term, .formula_design_display_names(term))
    }
    rows[[length(rows) + 1L]] <- data.frame(
      alias      = unique(aliases),
      parameter  = parameter,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  intercept <- term_map[["parameter"]][term_map[["term"]] == "intercept"]
  if (length(intercept) == 1L) {
    rows[[length(rows) + 1L]] <- data.frame(
      alias      = "mu",
      parameter  = intercept,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  aliases <- unique(do.call(rbind, rows))
  rownames(aliases) <- NULL
  return(aliases)
}


.hypothesis_marginal_means_aliases <- function(alias_catalog, parameter) {

  rows <- alias_catalog[alias_catalog[["parameter"]] == parameter, , drop = FALSE]
  aliases <- as.list(rows[["parameter"]])
  names(aliases) <- rows[["alias"]]
  return(aliases)
}


.hypothesis_marginal_means_attach_iwmde <- function(
    object, parameter, parameter_label, hypothesis, density_method,
    density_control, inference_type = "conditional") {

  point_refs <- .hypothesis_brma_point_refs(hypothesis, parameter)
  if (nrow(point_refs) == 0L) {
    return(object)
  }
  samples <- object[["inference"]][[inference_type]][[parameter]]
  point_refs <- .hypothesis_marginal_means_resolve_singleton_levels(
    samples    = samples,
    point_refs = point_refs
  )

  .hypothesis_marginal_means_check_point_refs(
    object          = object,
    parameter       = parameter,
    parameter_label = parameter_label,
    point_refs      = point_refs,
    inference_type  = inference_type
  )

  source_object <- object[["source_object"]]
  if (is.null(source_object) ||
      !inherits(source_object, "brma") ||
      is.null(source_object[["fit"]])) {
    stop(
      "The marginal-means object does not contain the source fitted brma ",
      "object needed to compute ", density_method, " ordinates.",
      call. = FALSE
    )
  }
  .check_iwmde_available(source_object, "qCMDE/IWMDE marginal-means hypothesis()")
  .iwmde_check_point_ordinate_supported(
    source_object,
    density_method
  )

  resources      <- .iwmde_workspace_resources(
    source_object,
    density_control[["workspace"]]
  )
  context        <- resources[["context"]]
  estimate_cache <- resources[["estimate_cache"]]
  object[["inference"]][[inference_type]][[parameter]] <-
    .hypothesis_brma_keep_requested_ordinates(
      posterior  = samples,
      point_refs = point_refs
    )

  for (i in seq_len(nrow(point_refs))) {
    object <- .hypothesis_marginal_means_attach_iwmde_ref(
      object              = object,
      context             = context,
      estimate_cache      = estimate_cache,
      parameter           = parameter,
      ref                 = point_refs[i, , drop = FALSE],
      inference_type      = inference_type,
      density_method      = density_method,
      density_control     = density_control
    )
  }

  return(object)
}


.hypothesis_marginal_means_resolve_singleton_levels <- function(
    samples, point_refs) {

  missing_level <- is.na(point_refs[["level"]])
  if (!is.list(samples) || length(samples) != 1L || !any(missing_level)) {
    return(point_refs)
  }
  level <- names(samples)
  if (length(level) != 1L || is.na(level) || !nzchar(level)) {
    return(point_refs)
  }
  point_refs[["level"]][missing_level] <- level

  return(point_refs)
}


.hypothesis_marginal_means_check_point_refs <- function(
    object, parameter, parameter_label, point_refs, inference_type) {

  samples <- object[["inference"]][[inference_type]][[parameter]]
  for (i in seq_len(nrow(point_refs))) {
    ref <- point_refs[i, , drop = FALSE]
    if (is.na(ref[["level"]]) && is.list(samples)) {
      stop(
        "qCMDE/IWMDE point hypotheses for marginal-means factor ",
        "parameters must specify a level, e.g. '", parameter_label,
        "[level] = ", ref[["value"]], "'.",
        call. = FALSE
      )
    }
    if (!is.na(ref[["level"]]) &&
        (!is.list(samples) || !ref[["level"]] %in% names(samples))) {
      stop("Hypothesis references unknown level '", ref[["level"]],
           "' for parameter '", parameter, "'.", call. = FALSE)
    }
  }

  invisible(TRUE)
}


.hypothesis_marginal_means_density_control <- function(object, density_method,
                                                       density_control) {

  if (is.null(density_control) && is.list(object[["density_settings"]])) {
    settings <- if (is.list(object[["ordinate_settings"]])) {
      object[["ordinate_settings"]]
    } else {
      object[["density_settings"]]
    }
    keep <- c(
      "n_points",
      "samples",
      "target_relative_mcse",
      "display_grid"
    )
    if (density_method %in% c("qCMDE", "IWMDE")) {
      keep <- c(keep, "normalization_points", "normalization_prob")
    }
    density_control <- settings[intersect(keep, names(settings))]
  }

  density_control <- .iwmde_density_control_resolve(
    density_method  = density_method,
    density_control = density_control,
    purpose         = "ordinate"
  )

  return(density_control)
}


.hypothesis_marginal_means_attach_iwmde_ref <- function(
    object, context, estimate_cache, parameter, ref, inference_type,
    density_method, density_control) {

  samples <- object[["inference"]][[inference_type]][[parameter]]
  if (is.null(samples)) {
    stop(
      "'marginal_means' object does not contain ", inference_type,
      " samples for '", parameter, "'.",
      call. = FALSE
    )
  }

  level <- ref[["level"]]
  value <- ref[["value"]]
  sample <- .hypothesis_marginal_means_ref_sample(
    samples = samples,
    level   = level
  )
  label <- if (is.na(level)) {
    parameter
  } else {
    paste0(parameter, "[", level, "]")
  }
  existing <- attr(sample, "posterior_ordinate", exact = TRUE)

  specs <- .iwmde_marginal_means_specs(
    marginal_means_object = object,
    parameter             = parameter,
    type                  = inference_type,
    levels                = if (is.na(level)) NULL else level
  )
  if (length(specs) != 1L) {
    stop(
      "Could not construct a qCMDE/IWMDE target for '", label, "'.",
      call. = FALSE
    )
  }
  spec <- specs[[1L]]
  metadata <- .iwmde_posterior_metadata(
    samples   = sample,
    parameter = parameter,
    level     = spec[["level"]]
  )
  expected_provenance <- .iwmde_request_provenance(
    context         = context,
    parameter       = spec[["label"]],
    density_method  = density_method,
    density_control = density_control,
    attribute       = "ordinate",
    value           = value,
    parameter_spec  = spec,
    metadata        = metadata
  )
  if (.iwmde_posterior_ordinate_matches_request(
        posterior_ordinate = existing,
        value              = value,
        provenance         = expected_provenance)) {
    return(object)
  }
  existing <- .iwmde_posterior_ordinate_drop_value(existing, value)

  estimate <- .iwmde_estimate(
    context         = context,
    parameter       = spec[["label"]],
    density_method  = density_method,
    density_control = density_control,
    outputs         = "ordinate",
    values          = value,
    parameter_spec  = spec,
    metadata        = metadata,
    cache           = estimate_cache
  )
  diagnostic <- estimate[["diagnostics"]][["ordinate"]]
  ordinate   <- estimate[["posterior_ordinate"]]
  if (is.null(ordinate)) {
    .iwmde_stop_ordinate_unavailable(
      message = .hypothesis_brma_iwmde_ordinate_failure_message(
        density_method = density_method,
        target         = paste0(label, " = ", value),
        diagnostic     = diagnostic,
        reason         = .hypothesis_brma_diagnostic_reason(diagnostic)
      ),
      estimate = estimate
    )
  }

  attr(sample, "posterior_ordinate") <- BayesTools::posterior_ordinate_append(
    existing = existing,
    ordinate = ordinate
  )
  samples <- .hypothesis_marginal_means_set_ref_sample(
    samples = samples,
    level   = level,
    sample  = sample
  )
  object[["inference"]][[inference_type]][[parameter]] <- samples

  return(object)
}


.hypothesis_marginal_means_ref_sample <- function(samples, level) {

  if (is.na(level)) {
    return(samples)
  }

  return(samples[[level]])
}


.hypothesis_marginal_means_set_ref_sample <- function(samples, level, sample) {

  if (is.na(level)) {
    return(sample)
  }

  samples[[level]] <- sample
  return(samples)
}
