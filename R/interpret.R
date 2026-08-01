#' @title Interpret brma Results
#'
#' @description Creates a concise textual interpretation of fitted RoBMA
#'   \code{brma} objects.
#'
#' @param object a fitted model object.
#' @param ... additional arguments passed to methods. The \code{brma} method
#' reserves \code{...}; unused arguments error, except deprecated
#' \code{output_scale}, which is accepted with a warning.
#'
#' @return A character vector with class \code{"interpret.brma"}. The
#'   normalized BayesTools interpretation records are stored in the
#'   \code{"records"} attribute; the \code{brma} method also stores
#'   \code{"scope"}, \code{"conditional"}, and optionally \code{"priors"}.
#'
#' @export
interpret <- function(object, ...) {

  UseMethod("interpret")
}

#' @rdname interpret
#' @export
interpret.default <- function(object, ...) {

  stop("'interpret' is available only for fitted brma objects.", call. = FALSE)
}

#' @rdname interpret
#' @param output_measure optional effect-size measure used for the pooled
#'   effect only. Supported conversion targets include \code{"SMD"},
#'   \code{"COR"}, \code{"ZCOR"}, and \code{"OR"} when the input scale allows
#'   conversion.
#' @param transform optional display transformation. \code{"EXP"} exponentiates
#'   log-scale \code{"OR"}, \code{"RR"}, \code{"HR"}, and \code{"IRR"}
#'   pooled effects; aliases for no transform are accepted.
#' @param conditional whether to summarize conditional estimates for RoBMA
#'   product-space objects. Defaults to \code{FALSE}.
#' @param scope character vector specifying sections to include. Use
#'   \code{"core"} for the default concise interpretation, \code{"all"} for all
#'   sections, or any of \code{"components"}, \code{"estimates"},
#'   \code{"moderators"}, \code{"scale"}, and \code{"bias"}.
#' @param probs two quantiles used for credible intervals. Defaults to
#'   \code{c(.025, .975)}.
#' @param central whether estimates are described by posterior mean or median.
#'   Defaults to \code{NULL}, which uses the posterior mean, except for
#'   \code{transform = "EXP"} pooled effects where the posterior mode is used.
#' @param priors whether to print prior distributions after the interpretation.
#'   Defaults to \code{FALSE}.
#' @param digits number of digits after the decimal point.
#'
#' @export
interpret.brma <- function(object, output_measure = NULL, transform = NULL,
                           conditional = FALSE, scope = "core",
                           probs = c(.025, .975),
                           central = NULL,
                           priors = FALSE,
                           digits = 3, ...) {

  if (is.null(object[["fit"]]) || length(object[["fit"]]) == 0L) {
    stop("'interpret' requires a fitted brma object.", call. = FALSE)
  }

  dots <- list(...)
  if ("output_scale" %in% names(dots)) {
    if (!is.null(output_measure)) {
      stop("Use only one of 'output_measure' and deprecated 'output_scale'.",
           call. = FALSE)
    }
    output_measure <- .interpret_output_scale(dots[["output_scale"]])
    dots[["output_scale"]] <- NULL
    warning(
      "'output_scale' is deprecated; use 'output_measure' instead.",
      call. = FALSE
    )
  }
  .check_unused_dots(dots, allowed = character(), caller = "interpret.brma()")

  BayesTools::check_bool(conditional, "conditional")
  if (conditional && !.is_RoBMA(object)) {
    stop("'conditional' interpretations are available only for RoBMA objects.",
         call. = FALSE)
  }

  scope <- .interpret_scope(scope)
  BayesTools::check_real(probs, "probs", lower = 0, upper = 1,
                         check_length = 2, allow_NA = FALSE)
  probs <- sort(probs)
  BayesTools::check_bool(priors, "priors")
  BayesTools::check_int(digits, "digits", lower = 0)

  effect_transform <- .effect_output_setup(
    object         = object,
    output_measure = output_measure,
    transform      = transform
  )
  central_explicit <- !is.null(central)
  if (is.null(central)) {
    central <- "mean"
  } else {
    central <- match.arg(central, c("mean", "median"))
  }
  effect_central <- if (
    !central_explicit && identical(effect_transform[["transform"]], "EXP")
  ) {
    "mode"
  } else {
    central
  }

  summary_object <- summary(
    object,
    probs                    = unique(c(probs[1], .5, probs[2])),
    include_mcmc_diagnostics = FALSE,
    conditional              = conditional
  )

  records <- .interpret_brma_records(
    object           = object,
    summary_object   = summary_object,
    effect_transform = effect_transform,
    conditional      = conditional,
    scope            = scope,
    probs            = probs,
    central          = central,
    effect_central   = effect_central
  )
  text <- .interpret_records_text(
    records  = records,
    digits   = digits,
    averaged = .is_RoBMA(object)
  )

  text <- text[nzchar(text)]
  class(text) <- c("interpret.brma", "character")
  attr(text, "scope")       <- scope
  attr(text, "conditional") <- conditional
  attr(text, "priors")      <- if (priors) .interpret_prior_text(object)
  attr(text, "records")     <- records

  return(text)
}

#' @rdname interpret
#' @param x an \code{interpret.brma} object.
#'
#' @export
print.interpret.brma <- function(x, ...) {

  cat(paste(unclass(x), collapse = "\n"))
  prior_text <- attr(x, "priors")
  if (length(prior_text) > 0L) {
    cat("\n\n")
    cat(paste(prior_text, collapse = "\n"))
  }
  cat("\n")

  return(invisible(x))
}

.interpret_prior_text <- function(object) {

  text <- capture.output(print_prior(object))
  text <- text[nzchar(text)]

  if (length(text) == 0L) {
    return(character())
  }

  return(c("Prior distributions:", text))
}

.interpret_scope <- function(scope) {

  choices <- c(
    "core", "all", "components", "estimates", "moderators", "scale", "bias"
  )
  BayesTools::check_char(scope, "scope", check_length = 0, allow_NA = FALSE)
  scope <- match.arg(scope, choices, several.ok = TRUE)

  if ("all" %in% scope) {
    scope <- c("components", "estimates", "moderators", "scale", "bias")
  }
  if ("core" %in% scope) {
    scope <- c(setdiff(scope, "core"), "components", "estimates")
  }

  return(unique(scope))
}

.interpret_output_scale <- function(output_scale) {

  if (is.null(output_scale)) {
    return(NULL)
  }

  BayesTools::check_char(output_scale, "output_scale", check_length = 1,
                         allow_NA = FALSE)
  output_scale <- toupper(gsub("[^A-Za-z0-9]", "", output_scale))
  aliases <- c(
    D          = "SMD",
    COHENSD    = "SMD",
    SMD        = "SMD",
    R          = "COR",
    COR        = "COR",
    CORRELATION = "COR",
    Z          = "ZCOR",
    ZCOR       = "ZCOR",
    FISHERSZ   = "ZCOR",
    LOGOR      = "OR",
    OR         = "OR",
    Y          = "GEN",
    NONE       = "GEN"
  )

  if (!output_scale %in% names(aliases)) {
    stop("Unknown deprecated 'output_scale'.", call. = FALSE)
  }

  return(unname(aliases[[output_scale]]))
}

.interpret_brma_records <- function(object, summary_object, effect_transform,
                                    conditional, scope, probs, central,
                                    effect_central) {

  sources <- .interpret_brma_sources(
    object           = object,
    summary_object   = summary_object,
    effect_transform = effect_transform,
    conditional      = conditional,
    scope            = scope,
    probs            = probs,
    central          = central,
    effect_central   = effect_central
  )
  plan <- .interpret_brma_plan(
    summary_object   = summary_object,
    effect_transform = effect_transform,
    scope            = scope,
    sources          = sources
  )

  .interpret_records_standard(
    sources = sources,
    plan    = plan
  )
}

.interpret_records_standard <- function(sources, plan) {

  if (!exists("interpret_records", envir = asNamespace("BayesTools"),
              inherits = FALSE)) {
    stop(
      "'interpret' requires a BayesTools version with interpret_records().",
      call. = FALSE
    )
  }

  return(BayesTools::interpret_records(
    sources = sources,
    plan    = plan,
    output  = "records",
    missing = "error"
  ))
}

.interpret_brma_sources <- function(object, summary_object, effect_transform,
                                    conditional, scope, probs, central,
                                    effect_central) {

  sources <- stats::setNames(list(), character())

  if (length(summary_object[["inclusion_components"]]) > 0L) {
    sources[["component_tests"]] <- .interpret_evidence_source(
      summary_object[["inclusion_components"]]
    )
  }
  if (length(summary_object[["inclusion_mods"]]) > 0L) {
    sources[["moderator_tests"]] <- .interpret_evidence_source(
      table = summary_object[["inclusion_mods"]],
      label = paste0("the location component '",
                     rownames(summary_object[["inclusion_mods"]]), "'")
    )
  }
  if (length(summary_object[["inclusion_scale"]]) > 0L) {
    sources[["scale_tests"]] <- .interpret_evidence_source(
      table = summary_object[["inclusion_scale"]],
      label = paste0("the scale component '",
                     rownames(summary_object[["inclusion_scale"]]), "'")
    )
  }

  if ("estimates" %in% scope) {
    effect <- tryCatch(
      pooled_effect(
        object         = object,
        output_measure = effect_transform[["output_measure"]],
        transform      = effect_transform[["transform"]],
        probs          = probs,
        conditional    = conditional
      ),
      error = function(e) NULL
    )
    if (!is.null(effect)) {
      sources[["pooled_effect"]] <- .interpret_record_source(
        .interpret_samples_record(
          samples      = effect,
          parameter    = effect_transform[["label"]],
          conditioning = if (conditional) {
            "conditional on effect inclusion"
          } else {
            NULL
          },
          probs        = probs,
          central      = effect_central
        )
      )
    }

    heterogeneity <- tryCatch(
      pooled_heterogeneity(
        object      = object,
        probs       = probs,
        conditional = conditional
      ),
      error = function(e) NULL
    )
    if (!is.null(heterogeneity)) {
      sources <- .interpret_add_heterogeneity_sources(
        sources       = sources,
        heterogeneity = heterogeneity,
        conditional   = conditional,
        probs         = probs,
        central       = central
      )
    }

    common_estimates <- .interpret_summary_estimates(summary_object, conditional)
    if (length(common_estimates) > 0L && "rho" %in% rownames(common_estimates)) {
      sources[["common_estimates"]] <- .interpret_estimate_source(
        table   = common_estimates,
        probs   = probs,
        central = central,
        units   = NULL
      )
    }
  }

  if ("moderators" %in% scope) {
    moderator_estimates <- .interpret_summary_estimates_mods(
      summary_object = summary_object,
      conditional    = conditional
    )
    if (length(moderator_estimates) > 0L) {
      sources[["moderator_estimates"]] <- .interpret_estimate_source(
        table   = moderator_estimates,
        probs   = probs,
        central = central,
        units   = if (.effect_output_requested(effect_transform)) {
          NULL
        } else {
          effect_transform[["label"]]
        }
      )
    }
  }

  if ("scale" %in% scope) {
    scale_estimates <- .interpret_summary_estimates_scale(
      summary_object = summary_object,
      conditional    = conditional
    )
    if (length(scale_estimates) > 0L) {
      sources[["scale_estimates"]] <- .interpret_estimate_source(
        table   = scale_estimates,
        probs   = probs,
        central = central,
        units   = NULL
      )
    }
  }

  if ("bias" %in% scope) {
    bias_estimates <- .interpret_summary_estimates_bias(
      summary_object = summary_object,
      conditional    = conditional
    )
    if (length(bias_estimates) > 0L) {
      sources[["bias_estimates"]] <- .interpret_estimate_source(
        table   = bias_estimates,
        probs   = probs,
        central = central,
        units   = NULL
      )
    }
  }

  return(sources)
}

.interpret_evidence_source <- function(table, label = NULL) {

  schema <- list(
    BF             = "inclusion_BF",
    BF_scale       = "linear",
    BF_orientation = "inclusion_over_exclusion"
  )
  if (!is.null(label)) {
    schema[["label"]] <- label
  }

  return(list(
    data   = table,
    schema = schema
  ))
}

.interpret_estimate_source <- function(table, probs, central, units) {

  row <- table[1, , drop = FALSE]
  schema <- list(
    central        = .interpret_central_column(row, central),
    central_name   = central,
    lower          = .interpret_probability_column(row, probs[1]),
    upper          = .interpret_probability_column(row, probs[2]),
    lower_prob     = probs[1],
    upper_prob     = probs[2]
  )
  if (!is.null(units)) {
    schema[["units"]] <- units
  }

  return(list(
    data   = table,
    schema = schema
  ))
}

.interpret_record_source <- function(record) {

  return(list(
    type = "record",
    data = record
  ))
}

.interpret_add_heterogeneity_sources <- function(sources, heterogeneity,
                                                 conditional, probs, central) {

  conditioning <- if (conditional) {
    "conditional on heterogeneity inclusion"
  } else {
    NULL
  }

  if (!is.list(heterogeneity) || is.matrix(heterogeneity)) {
    sources[["pooled_heterogeneity"]] <- .interpret_record_source(
      .interpret_samples_record(
        samples      = heterogeneity,
        parameter    = "tau",
        conditioning = conditioning,
        probs        = probs,
        central      = central
      )
    )
    return(sources)
  }

  component_labels <- .interpret_heterogeneity_component_labels(heterogeneity)
  component_ids    <- make.unique(make.names(component_labels))
  component_map    <- attr(sources, "heterogeneity_components")
  if (is.null(component_map)) {
    component_map <- character()
  }

  for (i in seq_along(heterogeneity)) {
    component_label <- component_labels[[i]]
    source_id       <- .interpret_heterogeneity_source_id(component_ids[[i]])
    component_map[[source_id]] <- component_label
    sources[[source_id]] <- .interpret_record_source(
      .interpret_samples_record(
        samples      = heterogeneity[[i]],
        parameter    = paste0("tau [", component_label, "]"),
        conditioning = conditioning,
        probs        = probs,
        central      = central
      )
    )
  }
  attr(sources, "heterogeneity_components") <- component_map

  return(sources)
}

.interpret_heterogeneity_component_labels <- function(heterogeneity) {

  component_labels <- names(heterogeneity)
  if (is.null(component_labels)) {
    component_labels <- paste0("component_", seq_along(heterogeneity))
  }
  missing_labels <- is.na(component_labels) | !nzchar(component_labels)
  component_labels[missing_labels] <- paste0(
    "component_",
    which(missing_labels)
  )

  return(component_labels)
}

.interpret_heterogeneity_source_id <- function(component_id) {

  paste0("pooled_heterogeneity.", component_id)
}

.interpret_heterogeneity_component_id_from_source <- function(source) {

  sub("^pooled_heterogeneity\\.", "", source)
}

.interpret_heterogeneity_component_from_source <- function(source, sources) {

  component_map <- attr(sources, "heterogeneity_components")
  if (!is.null(component_map) && source %in% names(component_map)) {
    return(unname(component_map[[source]]))
  }

  .interpret_heterogeneity_component_id_from_source(source)
}

.interpret_samples_record <- function(samples, parameter, conditioning,
                                      probs, central) {

  table <- summary(samples, probs = probs)
  if (central == "mode") {
    table <- cbind(
      Mode = .interpret_posterior_mode(as.numeric(samples[, 1L])),
      table
    )
  }
  row <- table[1, , drop = FALSE]

  central_col <- .interpret_central_column(row, central)
  lower_col   <- .interpret_probability_column(row, probs[1])
  upper_col   <- .interpret_probability_column(row, probs[2])

  return(list(
    kind           = "estimate",
    parameter      = parameter,
    central_name   = central,
    central_value  = as.numeric(row[1, central_col]),
    lower_value    = as.numeric(row[1, lower_col]),
    upper_value    = as.numeric(row[1, upper_col]),
    lower_prob     = probs[1],
    upper_prob     = probs[2],
    conditioning   = conditioning
  ))
}

.interpret_summary_estimates <- function(summary_object, conditional) {

  if (conditional) {
    return(summary_object[["estimates_conditional"]])
  }

  return(summary_object[["estimates"]])
}

.interpret_summary_estimates_mods <- function(summary_object, conditional) {

  if (conditional) {
    return(summary_object[["estimates_mods_conditional"]])
  }

  return(summary_object[["estimates_mods"]])
}

.interpret_summary_estimates_scale <- function(summary_object, conditional) {

  if (conditional) {
    return(summary_object[["estimates_scale_conditional"]])
  }

  return(summary_object[["estimates_scale"]])
}

.interpret_summary_estimates_bias <- function(summary_object, conditional) {

  if (conditional) {
    return(summary_object[["estimates_bias_conditional"]])
  }

  return(summary_object[["estimates_bias"]])
}

.interpret_brma_plan <- function(summary_object, effect_transform, scope,
                                 sources) {

  plan <- list(list(
    kind    = "header",
    section = "model",
    item_id = "header",
    order   = 0,
    text    = paste0(summary_object[["name"]], ".")
  ))

  if ("components" %in% scope && "component_tests" %in% names(sources)) {
    plan <- .interpret_add_component_plan(
      plan  = plan,
      table = summary_object[["inclusion_components"]],
      order = 10
    )
  }

  if ("estimates" %in% scope) {
    if ("pooled_effect" %in% names(sources)) {
      plan[[length(plan) + 1L]] <- list(
        kind    = "estimate",
        section = "estimates",
        item_id = "effect",
        order   = 100,
        source  = "pooled_effect",
        label   = "pooled effect"
      )
    }
    if ("pooled_heterogeneity" %in% names(sources)) {
      plan[[length(plan) + 1L]] <- list(
        kind    = "estimate",
        section = "estimates",
        item_id = "heterogeneity",
        order   = 110,
        source  = "pooled_heterogeneity",
        label   = "pooled heterogeneity"
      )
    }
    heterogeneity_sources <- grep(
      "^pooled_heterogeneity\\.",
      names(sources),
      value = TRUE
    )
    for (i in seq_along(heterogeneity_sources)) {
      source          <- heterogeneity_sources[[i]]
      component_id    <- .interpret_heterogeneity_component_id_from_source(source)
      component_label <- .interpret_heterogeneity_component_from_source(
        source  = source,
        sources = sources
      )
      plan[[length(plan) + 1L]] <- list(
        kind    = "estimate",
        section = "estimates",
        item_id = paste0("heterogeneity.", component_id),
        order   = 110 + i,
        source  = source,
        label   = paste0("Pooled heterogeneity (", component_label, ")")
      )
    }
    if ("common_estimates" %in% names(sources)) {
      plan[[length(plan) + 1L]] <- list(
        kind      = "estimate",
        section   = "estimates",
        item_id   = "rho",
        order     = 190,
        source    = "common_estimates",
        row       = "rho",
        label     = "cluster correlation",
        parameter = "rho"
      )
    }
  }

  if ("moderators" %in% scope) {
    if ("moderator_tests" %in% names(sources)) {
      plan[[length(plan) + 1L]] <- list(
        kind     = "for_each",
        section  = "moderators",
        item_id  = "location_inclusion",
        order    = 200,
        source   = "moderator_tests",
        rows     = "source_order",
        template = "evidence",
        BF_name  = "BF"
      )
    }
    if ("moderator_estimates" %in% names(sources)) {
      plan[[length(plan) + 1L]] <- list(
        kind     = "for_each",
        section  = "moderators",
        item_id  = "location_estimate",
        order    = 300,
        source   = "moderator_estimates",
        rows     = "source_order",
        template = "estimate"
      )
    }
  }

  if ("scale" %in% scope) {
    if ("scale_tests" %in% names(sources)) {
      plan[[length(plan) + 1L]] <- list(
        kind     = "for_each",
        section  = "scale",
        item_id  = "scale_inclusion",
        order    = 400,
        source   = "scale_tests",
        rows     = "source_order",
        template = "evidence",
        BF_name  = "BF"
      )
    }
    if ("scale_estimates" %in% names(sources)) {
      plan[[length(plan) + 1L]] <- list(
        kind     = "for_each",
        section  = "scale",
        item_id  = "scale_estimate",
        order    = 500,
        source   = "scale_estimates",
        rows     = "source_order",
        template = "estimate"
      )
    }
  }

  if ("bias" %in% scope) {
    if (!"components" %in% scope && "component_tests" %in% names(sources) &&
        length(summary_object[["inclusion_components"]]) > 0L &&
        "Publication Bias" %in%
        rownames(summary_object[["inclusion_components"]])) {
      plan[[length(plan) + 1L]] <- .interpret_component_plan_item(
        row   = "Publication Bias",
        order = 600
      )
      plan[[length(plan)]][["section"]] <- "bias"
    }
    if ("bias_estimates" %in% names(sources)) {
      plan[[length(plan) + 1L]] <- list(
        kind     = "for_each",
        section  = "bias",
        item_id  = "publication_bias_estimate",
        order    = 700,
        source   = "bias_estimates",
        rows     = "source_order",
        template = "estimate"
      )
    }
  }

  if (.effect_output_requested(effect_transform) &&
      !is.null(effect_transform[["note"]])) {
    plan[[length(plan) + 1L]] <- list(
      kind    = "note",
      section = "notes",
      item_id = "effect_transform",
      order   = 900,
      text    = .interpret_effect_transform_note(effect_transform)
    )
  }

  return(plan)
}

.interpret_effect_transform_note <- function(effect_transform) {

  note <- effect_transform[["note"]]
  if (is.null(note)) {
    return(NULL)
  }

  return(sub("^Effect estimates", "Pooled effect estimates", note))
}

.interpret_add_component_plan <- function(plan, table, order) {

  rows <- rownames(table)
  for (i in seq_along(rows)) {
    plan[[length(plan) + 1L]] <- .interpret_component_plan_item(
      row   = rows[[i]],
      order = order + i
    )
  }

  return(plan)
}

.interpret_component_plan_item <- function(row, order) {

  BF_names <- c(
    "Effect"           = "BF10",
    "Heterogeneity"    = "BFrf",
    "Publication Bias" = "BFpb"
  )
  labels <- c(
    "Effect"           = "the effect",
    "Heterogeneity"    = "heterogeneity",
    "Publication Bias" = "publication bias"
  )
  item_ids <- c(
    "Effect"           = "effect",
    "Heterogeneity"    = "heterogeneity",
    "Publication Bias" = "publication_bias"
  )

  return(list(
    kind    = "evidence",
    section = "components",
    item_id = if (!is.null(item_ids[[row]])) unname(item_ids[[row]]) else row,
    order   = order,
    source  = "component_tests",
    row     = row,
    label   = if (!is.null(labels[[row]])) unname(labels[[row]]) else tolower(row),
    BF_name = if (!is.null(BF_names[[row]])) unname(BF_names[[row]]) else "BF"
  ))
}

.interpret_records_text <- function(records, digits, averaged) {

  if (nrow(records) == 0L) {
    return(character())
  }

  output <- vapply(seq_len(nrow(records)), function(i) {

    record <- records[i, , drop = FALSE]
    kind   <- record[["kind"]]

    switch(
      kind,
      header   = record[["text"]],
      note     = record[["text"]],
      prior    = record[["text"]],
      footnote = record[["text"]],
      evidence = .interpret_evidence_record_sentence(record, digits),
      estimate = .interpret_estimate_record_sentence(record, digits, averaged),
      ""
    )
  }, character(1))

  names(output) <- vapply(seq_len(nrow(records)), function(i) {
    .interpret_record_output_name(records[i, , drop = FALSE])
  }, character(1))

  return(output)
}

.interpret_evidence_record_sentence <- function(record, digits) {

  BF <- record[["BF_canonical_value"]]
  BF_bound_operator <- .interpret_record_value(
    record[["BF_canonical_bound_operator"]]
  )
  BF_name <- .interpret_record_value(record[["BF_name"]])
  if (is.null(BF_name)) {
    BF_name <- "BF"
  }
  display_BF <- .interpret_display_BF(BF, BF_name, BF_bound_operator)

  text <- paste0(
    .interpret_evidence_prefix(record), ": ",
    .interpret_evidence(
      BF             = BF,
      item           = record[["label"]],
      bound_operator = BF_bound_operator
    ),
    " (", display_BF[["name"]], " ",
    if (is.null(display_BF[["operator"]])) "=" else display_BF[["operator"]],
    " ",
    .interpret_number(display_BF[["value"]], digits),
    ")."
  )

  return(text)
}

.interpret_estimate_record_sentence <- function(record, digits, averaged) {

  central_name <- record[["central_name"]]
  if (is.na(central_name) || !nzchar(central_name)) {
    central_name <- "estimate"
  }
  central_label <- paste(
    c(if (averaged) "model-averaged", central_name),
    collapse = " "
  )

  text <- paste0(
    .interpret_estimate_prefix(record), ": ",
    central_label, " ", record[["parameter"]], " = ",
    .interpret_number(record[["central_value"]], digits),
    ", ", .interpret_interval_label(c(
      record[["lower_prob"]],
      record[["upper_prob"]]
    )),
    " CrI [",
    .interpret_number(record[["lower_value"]], digits),
    ", ",
    .interpret_number(record[["upper_value"]], digits),
    "]"
  )

  units <- .interpret_record_value(record[["units"]])
  if (!is.null(units)) {
    text <- paste0(text, " on the ", units, " scale")
  }

  conditioning <- .interpret_record_value(record[["conditioning"]])
  if (!is.null(conditioning)) {
    text <- paste0(text, " (", conditioning, ")")
  }

  return(paste0(text, "."))
}

.interpret_evidence_prefix <- function(record) {

  section <- record[["section"]]
  row     <- record[["row"]]

  if (section == "components" || section == "bias") {
    return(paste0(row, " inclusion"))
  }
  if (section == "moderators") {
    return(paste0("Location inclusion (", row, ")"))
  }
  if (section == "scale") {
    return(paste0("Scale inclusion (", row, ")"))
  }

  return(record[["label"]])
}

.interpret_estimate_prefix <- function(record) {

  section <- record[["section"]]
  item_id <- record[["item_id"]]
  row     <- record[["row"]]

  if (section == "estimates" && item_id == "effect") {
    return("Pooled effect")
  }
  if (section == "estimates" && item_id == "heterogeneity") {
    return("Pooled heterogeneity")
  }
  if (section == "estimates" && startsWith(item_id, "heterogeneity.")) {
    return(record[["label"]])
  }
  if (section == "estimates" && item_id == "rho") {
    return("Cluster correlation")
  }
  if (section == "moderators") {
    return(paste0("Location estimate (", row, ")"))
  }
  if (section == "scale") {
    return(paste0("Scale estimate (", row, ")"))
  }
  if (section == "bias") {
    return(paste0("Publication-bias estimate (", row, ")"))
  }

  return(record[["label"]])
}

.interpret_record_output_name <- function(record) {

  kind    <- record[["kind"]]
  section <- record[["section"]]
  item_id <- record[["item_id"]]
  row     <- record[["row"]]

  if (kind == "header" && section == "model") {
    return("model")
  }
  if (kind == "note") {
    return("note")
  }
  if (section == "components" || (section == "bias" && kind == "evidence")) {
    return(paste0("components.", make.names(tolower(row))))
  }
  if (section == "estimates") {
    return(item_id)
  }
  if (section == "moderators" && kind == "evidence") {
    return(paste0("location.inclusion.", make.names(row)))
  }
  if (section == "moderators" && kind == "estimate") {
    return(paste0("location.estimate.", make.names(row)))
  }
  if (section == "scale" && kind == "evidence") {
    return(paste0("scale.inclusion.", make.names(row)))
  }
  if (section == "scale" && kind == "estimate") {
    return(paste0("scale.estimate.", make.names(row)))
  }
  if (section == "bias" && kind == "estimate") {
    return(paste0("publication.bias.estimate.", make.names(row)))
  }

  return(record[["record_id"]])
}

.interpret_record_value <- function(x) {

  if (length(x) == 0L || is.null(x) || is.na(x) || !nzchar(x)) {
    return(NULL)
  }

  return(x)
}

.interpret_display_BF <- function(BF, BF_name, BF_bound_operator = NULL) {

  if (is.finite(BF) && BF > 0 && BF < 1) {
    return(list(
      name     = .interpret_inverse_BF_name(BF_name),
      value    = 1 / BF,
      operator = .interpret_inverse_BF_bound_operator(BF_bound_operator)
    ))
  }
  if (identical(BF, 0)) {
    return(list(
      name     = .interpret_inverse_BF_name(BF_name),
      value    = Inf,
      operator = .interpret_inverse_BF_bound_operator(BF_bound_operator)
    ))
  }

  return(list(
    name     = BF_name,
    value    = BF,
    operator = BF_bound_operator
  ))
}

.interpret_inverse_BF_name <- function(BF_name) {

  if (BF_name == "BF") {
    return(BF_name)
  }
  if (BF_name == "Inclusion BF") {
    return("Exclusion BF")
  }
  if (BF_name == "Exclusion BF") {
    return("Inclusion BF")
  }
  if (nchar(BF_name) == 4L && substr(BF_name, 1L, 2L) == "BF") {
    return(paste0("BF", substr(BF_name, 4L, 4L), substr(BF_name, 3L, 3L)))
  }

  return(paste0(BF_name, "^-1"))
}

.interpret_inverse_BF_bound_operator <- function(bound_operator) {

  if (is.null(bound_operator)) {
    return(NULL)
  }
  if (bound_operator == "<") {
    return(">")
  }
  if (bound_operator == ">") {
    return("<")
  }

  return(bound_operator)
}

.interpret_evidence <- function(BF, item, bound_operator = NULL) {

  if (!is.finite(BF) && !is.infinite(BF)) {
    return(paste0("evidence for ", item, " unavailable"))
  }
  if (isTRUE(all.equal(BF, 1)) && is.null(bound_operator)) {
    return(paste0("no evidence for or against ", item))
  }

  strength <- if (abs(log(BF)) > log(10)) {
    "strong"
  } else if (abs(log(BF)) > log(3)) {
    "moderate"
  } else {
    "weak"
  }
  direction <- if (BF > 1 || (isTRUE(all.equal(BF, 1)) &&
                              identical(bound_operator, ">"))) {
    "in favor of"
  } else {
    "against"
  }

  return(paste(strength, "evidence", direction, item))
}

.interpret_central_column <- function(row, central) {

  if (central == "mode" && "Mode" %in% colnames(row)) {
    return("Mode")
  }
  if (central == "mean" && "Mean" %in% colnames(row)) {
    return("Mean")
  }
  if (central == "median" && "Median" %in% colnames(row)) {
    return("Median")
  }
  if (central == "median" && "0.5" %in% colnames(row)) {
    return("0.5")
  }

  stop("Could not find requested central estimate in summary table.",
       call. = FALSE)
}

.interpret_posterior_mode <- function(x) {

  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }

  unique_x <- unique(x)
  if (length(unique_x) == 1L) {
    return(unique_x)
  }

  density <- stats::density(x)
  return(density[["x"]][which.max(density[["y"]])])
}

.interpret_probability_column <- function(row, prob) {

  column <- as.character(prob)
  if (column %in% colnames(row)) {
    return(column)
  }

  numeric_columns <- suppressWarnings(as.numeric(colnames(row)))
  index <- which(
    !is.na(numeric_columns) &
      abs(numeric_columns - prob) < sqrt(.Machine$double.eps)
  )
  if (length(index) == 1L) {
    return(colnames(row)[index])
  }

  stop("Could not find requested credible interval column in summary table.",
       call. = FALSE)
}

.interpret_interval_label <- function(probs) {

  level <- round(100 * diff(range(probs)), 1)
  if (isTRUE(all.equal(level, round(level)))) {
    level <- round(level)
  }

  return(paste0(level, "%"))
}

.interpret_number <- function(x, digits) {

  if (is.na(x)) {
    return("NA")
  }
  if (is.infinite(x)) {
    return(if (x > 0) "Inf" else "-Inf")
  }
  if (abs(x) < 0.5 * 10^(-digits)) {
    x <- 0
  }

  return(formatC(x, format = "f", digits = digits))
}
