#' @rdname hypothesis
#'
#' @param parameter optional marginal-means parameter/term used to disambiguate
#' the hypothesis expression.
#'
#' @export
hypothesis.marginal_means.brma <- function(object, hypothesis,
                                           parameter = NULL,
                                           logBF = FALSE, BF01 = FALSE,
                                           seed = NULL,
                                           density_method = c(
                                             "KDE", "qCMDE", "IWMDE"
                                           ),
                                           density_control = NULL,
                                           columns = "default",
                                           ...) {

  BayesTools::check_char(hypothesis, "hypothesis", check_length = 0, allow_NA = FALSE)
  BayesTools::check_bool(logBF, "logBF")
  BayesTools::check_bool(BF01, "BF01")
  BayesTools::check_real(seed, "seed", check_length = 1, allow_NULL = TRUE, allow_NA = FALSE)
  BayesTools::check_char(columns, "columns", check_length = 0, allow_NA = FALSE)
  dots <- list(...)
  if (identical(dots[["type"]], "conditional")) {
    dots[["type"]] <- NULL
  }
  .warn_unused_dots(
    dots    = dots,
    allowed = character(),
    caller  = "hypothesis.marginal_means()"
  )

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

  density_method       <- .density_method_normalize(density_method)
  precomputed_density <- density_method %in% c("qCMDE", "IWMDE")
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
    object <- .hypothesis_marginal_means_attach_iwmde(
      object          = object,
      parameter       = parameter,
      parameter_label = parameter_label,
      hypothesis      = hypothesis,
      density_method  = density_method,
      density_control = density_control
    )
    density_method <- "precomputed"
  } else if (!is.null(density_control)) {
    stop("'density_control' is only used when 'density_method' is ",
         "'qCMDE' or 'IWMDE'.", call. = FALSE)
  }

  inference <- object[["inference"]]
  if (is.null(inference[["conditional"]])) {
    stop(
      "'marginal_means' object does not contain alternative-conditioned ",
      "marginal means.",
      call. = FALSE
    )
  }
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


.hypothesis_marginal_means_attach_iwmde <- function(
    object, parameter, parameter_label, hypothesis, density_method,
    density_control) {

  point_refs <- .hypothesis_brma_point_refs(hypothesis, parameter)
  if (nrow(point_refs) == 0L) {
    return(object)
  }

  .hypothesis_marginal_means_check_point_refs(
    object          = object,
    parameter       = parameter,
    parameter_label = parameter_label,
    point_refs      = point_refs
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

  context        <- .iwmde_context(source_object)
  estimate_cache <- .iwmde_estimate_cache()

  for (i in seq_len(nrow(point_refs))) {
    object <- .hypothesis_marginal_means_attach_iwmde_ref(
      object              = object,
      context             = context,
      estimate_cache      = estimate_cache,
      parameter           = parameter,
      ref                 = point_refs[i, , drop = FALSE],
      density_method      = density_method,
      density_control     = density_control
    )
  }

  return(object)
}


.hypothesis_marginal_means_check_point_refs <- function(
    object, parameter, parameter_label, point_refs) {

  samples <- object[["inference"]][["conditional"]][[parameter]]
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

  if (is.null(density_control) && is.list(object[["iwmde_settings"]])) {
    settings <- if (is.list(object[["iwmde_ordinate_settings"]])) {
      object[["iwmde_ordinate_settings"]]
    } else {
      object[["iwmde_settings"]]
    }
    keep <- c(
      "n_points",
      "max_samples",
      "initial_samples",
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
    object, context, estimate_cache, parameter, ref, density_method,
    density_control) {

  samples <- object[["inference"]][["conditional"]][[parameter]]
  if (is.null(samples)) {
    stop(
      "'marginal_means' object does not contain alternative-conditioned ",
      "samples for '", parameter, "'.",
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
    type                  = "conditional",
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
      message = paste0(
        density_method, " posterior ordinate is unavailable for '", label,
        " = ", value, "': ", .hypothesis_brma_diagnostic_reason(diagnostic)
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
  object[["inference"]][["conditional"]][[parameter]] <- samples

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
