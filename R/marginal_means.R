# ============================================================================ #
# brma.marginal_means.R
# ============================================================================ #
#
# Estimated marginal means for brma objects with moderators.
#
# The expensive computation is delegated to BayesTools::as_marginal_inference()
# and stored in a small RoBMA-side result object. Summary and plot methods only
# format or render the stored BayesTools marginal-inference object.
#
# ============================================================================ #


#' @title Estimated Marginal Means
#'
#' @description Computes estimated marginal means for fitted \code{brma}
#' objects with moderators. The returned object stores the BayesTools
#' marginal-inference object and can be summarized or plotted with
#' \code{summary()} and \code{plot()}.
#'
#' @param object a fitted model object.
#' @param ... additional arguments passed to methods.
#'
#' @return A method-specific estimated marginal means object.
#'
#' @seealso [summary()], [plot()], [summary.brma()], [regplot()]
#' @export
marginal_means <- function(object, ...) {

  UseMethod("marginal_means")
}


#' @title Estimated Marginal Means for brma Objects
#'
#' @description Computes estimated marginal means for a fitted \code{brma}
#' object with moderators using \code{BayesTools::as_marginal_inference()}.
#'
#' @param object a fitted \code{brma} object with moderators.
#' @param null_hypothesis point null hypothesis used for inclusion Bayes
#' factors. Defaults to \code{0}.
#' @param normal_approximation whether prior and posterior density at the null
#' should be approximated with a normal distribution. Defaults to \code{FALSE}.
#' @param n_samples number of samples/grid points used by BayesTools for
#' marginal prior densities. Defaults to \code{10000}.
#' @param bf whether inclusion Bayes factors should be shown by default in
#' summaries. Defaults to \code{TRUE} for RoBMA/BMA objects and \code{FALSE}
#' for single-model \code{brma} objects.
#' @inheritParams predict.brma
#' @param ... additional arguments (currently ignored).
#'
#' @return A list of class \code{marginal_means.brma} containing the
#' BayesTools \code{marginal_inference} object and parameter metadata.
#'
#' @examples \dontrun{
#' fit <- brma(yi = yi, sei = sei, mods = ~ x, data = dat)
#' mm  <- marginal_means(fit)
#' summary(mm)
#' plot(mm, parameter = "x")
#' }
#'
#' @seealso [summary()], [plot()], [summary.brma()], [regplot()]
#' @export
marginal_means.brma <- function(object, null_hypothesis = 0,
                                normal_approximation = FALSE,
                                n_samples = 10000,
                                output_measure = NULL, transform = NULL,
                                bf = NULL, ...) {

  if (!.is_mods(object)) {
    stop("'marginal_means' requires a model with moderators.", call. = FALSE)
  }
  if (is.null(object[["fit"]]) || length(object[["fit"]]) == 0L) {
    stop("'marginal_means' requires a fitted brma object.", call. = FALSE)
  }

  BayesTools::check_real(null_hypothesis, "null_hypothesis", check_length = 1)
  BayesTools::check_bool(normal_approximation, "normal_approximation")
  BayesTools::check_int(n_samples, "n_samples", lower = 2)
  model_averaged <- .is_RoBMA(object)
  bf             <- .marginal_means_resolve_bf(model_averaged, bf)
  effect_transform <- .effect_output_setup(
    object         = object,
    output_measure = output_measure,
    transform      = transform
  )

  formula          <- .create_fit_formula_list(data = object[["data"]], parameter = "mods")
  terms            <- .marginal_means_terms(formula)
  parameters       <- BayesTools::JAGS_parameter_names(
    parameters        = terms,
    formula_parameter = "mu"
  )
  conditional_list <- .marginal_means_conditional_list(
    terms      = terms,
    parameters = parameters
  )

  inference <- suppressWarnings(BayesTools::as_marginal_inference(
    model                = object[["fit"]],
    marginal_parameters  = parameters,
    parameters           = parameters,
    conditional_list     = conditional_list,
    conditional_rule     = "OR",
    formula              = formula,
    null_hypothesis      = null_hypothesis,
    normal_approximation = normal_approximation,
    n_samples            = n_samples,
    silent               = TRUE,
    force_plots          = TRUE
  ))

  available_parameters <- Reduce(
    intersect,
    list(
      parameters,
      names(inference[["averaged"]]),
      names(inference[["conditional"]]),
      names(inference[["inference"]])
    )
  )

  if (length(available_parameters) == 0L) {
    stop("No marginal means are available for this model.", call. = FALSE)
  }

  term_map <- data.frame(
    term             = terms,
    parameter        = parameters,
    label            = terms,
    check.names      = FALSE,
    stringsAsFactors = FALSE
  )
  term_map <- term_map[term_map[["parameter"]] %in% available_parameters, , drop = FALSE]

  output <- list(
    inference              = inference,
    parameters             = available_parameters,
    term_map               = term_map,
    formula                = formula,
    null_hypothesis        = null_hypothesis,
    normal_approximation   = normal_approximation,
    n_samples              = n_samples,
    input_measure          = .measure(object),
    effect_transform       = effect_transform,
    model_averaged         = model_averaged,
    bf                     = bf,
    name                   = .summary.brma_model_names(object)
  )

  class(output) <- c("marginal_means.brma", "marginal_means")

  return(output)
}


#' @title Summarize Estimated Marginal Means
#'
#' @description Summarizes estimated marginal means stored in a
#' \code{marginal_means.brma} object using
#' \code{BayesTools::marginal_estimates_table()}.
#'
#' @param object a \code{marginal_means.brma} object.
#' @param type for RoBMA product-space objects, whether to summarize
#' model-averaged (\code{"averaged"}) or conditional (\code{"conditional"})
#' marginal means. Defaults to \code{"averaged"} and is available only for
#' RoBMA marginal means.
#' @param probs quantiles of the posterior distribution to be displayed.
#' Defaults to \code{c(.025, .50, .975)}.
#' @param logBF whether to show inclusion Bayes factors on the log scale.
#' Defaults to \code{FALSE}.
#' @param BF01 whether to show inverse inclusion Bayes factors. Defaults to
#' \code{FALSE}.
#' @param bf whether to show inclusion Bayes factors. Defaults to the setting
#' stored by \code{marginal_means()}.
#' @inheritParams predict.brma
#' @param ... additional arguments (currently ignored).
#'
#' @return A \code{BayesTools_table} of class
#' \code{summary.marginal_means.brma}.
#'
#' @export
summary.marginal_means.brma <- function(object, type = NULL,
                                        probs = c(.025, .50, .975),
                                        logBF = FALSE, BF01 = FALSE,
                                        bf = NULL,
                                        output_measure, transform, ...) {

  type <- .marginal_means_type(object = object, type = type)
  BayesTools::check_real(probs, "probs", allow_NULL = TRUE, check_length = 0)
  BayesTools::check_bool(logBF, "logBF")
  BayesTools::check_bool(BF01, "BF01")
  bf <- .marginal_means_resolve_bf(object[["model_averaged"]], bf, object[["bf"]])

  if (missing(output_measure) && missing(transform)) {
    effect_transform <- object[["effect_transform"]]
  } else {
    effect_transform <- .effect_output_setup_measure(
      input_measure  = object[["input_measure"]],
      output_measure = if (missing(output_measure)) NULL else output_measure,
      transform      = if (missing(transform)) NULL else transform
    )
  }

  samples    <- .transform_marginal_samples_effect(
    samples          = object[["inference"]][[type]],
    effect_transform = effect_transform
  )
  inference  <- object[["inference"]][["inference"]]
  parameters <- object[["parameters"]]
  parameters <- parameters[
    parameters %in% names(samples) &
      parameters %in% names(inference)
  ]

  if (length(parameters) == 0L) {
    stop("No marginal means are available for type = '", type, "'.",
         call. = FALSE)
  }

  estimates <- BayesTools::marginal_estimates_table(
    samples        = samples,
    inference      = inference,
    parameters     = parameters,
    probs          = probs,
    logBF          = logBF,
    BF01           = BF01,
    formula_prefix = FALSE,
    title          = .effect_output_title(
      title = if (isTRUE(object[["model_averaged"]])) {
        switch(
          type,
          "averaged"    = "Model-Averaged Marginal Means:",
          "conditional" = "Conditional Marginal Means:"
        )
      } else {
        "Marginal Means:"
      },
      effect_transform = effect_transform
    ),
    footnotes      = effect_transform[["note"]]
  )

  if (!bf) {
    estimates <- .marginal_means_drop_bf(estimates)
  }

  class(estimates) <- c("summary.marginal_means.brma", class(estimates))
  attr(estimates, "marginal_type") <- type

  return(estimates)
}


#' @title Print Estimated Marginal Means
#'
#' @description Prints the model-averaged estimated marginal means summary.
#'
#' @param x a \code{marginal_means.brma} object.
#' @param ... additional arguments passed to \code{summary()}.
#'
#' @return Returns \code{x} invisibly.
#'
#' @export
print.marginal_means.brma <- function(x, ...) {

  print(summary(x, ...))

  return(invisible(x))
}


#' @title Print Summary of Estimated Marginal Means
#'
#' @description Prints a summary table of estimated marginal means.
#'
#' @param x a \code{summary.marginal_means.brma} object.
#' @param ... additional arguments (currently ignored).
#'
#' @return Returns \code{x} invisibly.
#'
#' @export
print.summary.marginal_means.brma <- function(x, ...) {

  class(x) <- setdiff(class(x), "summary.marginal_means.brma")

  cat("\n")
  print(x)
  cat("\n")

  return(invisible(x))
}


#' @title Plot Estimated Marginal Means
#'
#' @description Plots estimated marginal means stored in a
#' \code{marginal_means.brma} object using \code{BayesTools::plot_marginal()}.
#'
#' @param x a \code{marginal_means.brma} object.
#' @param parameter moderator term to plot. Use the original term name, for
#' example \code{"measure"}, or the internal parameter name, for example
#' \code{"mu_measure"}.
#' @param type for RoBMA product-space objects, whether to plot model-averaged
#' (\code{"averaged"}) or conditional (\code{"conditional"}) marginal means.
#' Defaults to \code{"averaged"} and is available only for RoBMA marginal
#' means.
#' @param prior whether the marginal prior distribution should be added to the
#' plot. Defaults to \code{FALSE}.
#' @param plot_type whether to use base R graphics (\code{"base"}) or ggplot2
#' (\code{"ggplot"}). Defaults to \code{"base"}.
#' @param dots_prior list of additional graphical arguments passed to the prior
#' plotting function.
#' @inheritParams predict.brma
#' @param ... additional graphical arguments passed to
#' \code{BayesTools::plot_marginal()}.
#'
#' @return \code{NULL} invisibly if \code{plot_type = "base"} or a ggplot object
#' if \code{plot_type = "ggplot"}.
#'
#' @export
plot.marginal_means.brma <- function(x, parameter, type = NULL,
                                     prior = FALSE, plot_type = "base",
                                     dots_prior = NULL,
                                     output_measure, transform, ...) {

  type <- .marginal_means_type(object = x, type = type)
  BayesTools::check_bool(prior, "prior")
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))

  if (missing(output_measure) && missing(transform)) {
    effect_transform <- x[["effect_transform"]]
  } else {
    effect_transform <- .effect_output_setup_measure(
      input_measure  = x[["input_measure"]],
      output_measure = if (missing(output_measure)) NULL else output_measure,
      transform      = if (missing(transform)) NULL else transform
    )
  }

  if (missing(parameter) || is.null(parameter)) {
    stop("The 'parameter' argument must be specified. Available terms are: ",
         .marginal_means_available_terms(x), ".", call. = FALSE)
  }

  selected <- .marginal_means_select_parameter(x, parameter)
  samples  <- x[["inference"]][[type]]

  if (!selected[["parameter"]] %in% names(samples)) {
    stop("No marginal means are available for parameter '",
         selected[["term"]], "' and type = '", type, "'.", call. = FALSE)
  }

  n_levels   <- length(samples[[selected[["parameter"]]]])
  dots       <- .set_dots_plot(..., n_levels = n_levels)
  term_label <- selected[["label"]]
  if (is.null(dots[["xlab"]])) {
    dots[["xlab"]] <- .plot_parameter_label("mu", effect_transform)
  }
  prior_component_counts <- .marginal_means_prior_component_counts(
    samples[[selected[["parameter"]]]]
  )
  dots_prior <- .set_dots_prior_marginal_means(
    dots_prior             = dots_prior,
    n_levels               = n_levels,
    prior_component_counts = prior_component_counts
  )

  args                          <- dots
  args$samples                  <- samples
  args$parameter                <- selected[["parameter"]]
  args$plot_type                <- plot_type
  args$prior                    <- prior
  args$n_points                 <- 1000
  args$transformation           <- .effect_plot_transformation(effect_transform)
  args$transformation_arguments <- NULL
  args$transformation_settings  <- FALSE
  args$par_name                 <- dots[["xlab"]]
  args$dots_prior               <- dots_prior
  if (plot_type == "base" && !identical(dots[["legend"]], FALSE)) {
    args$legend <- FALSE
  }

  plot <- suppressMessages(do.call(BayesTools::plot_marginal, args))

  if (plot_type == "base") {
    .marginal_means_add_base_legend(
      samples_parameter = samples[[selected[["parameter"]]]],
      term_label        = term_label,
      dots              = dots
    )
    return(invisible(plot))
  } else if (plot_type == "ggplot") {
    return(.marginal_means_set_ggplot_legend_title(
      plot       = plot,
      term_label = term_label
    ))
  }
}


# Extract model terms used by the moderator formula.
.marginal_means_terms <- function(formula) {

  formula_terms <- stats::terms(formula)
  terms         <- c(
    if (attr(formula_terms, "intercept") == 1L) "intercept",
    attr(formula_terms, "term.labels")
  )

  if (length(terms) == 0L) {
    stop("No moderator terms found.", call. = FALSE)
  }

  return(terms)
}


# Build BayesTools conditional-list specification for marginal means.
.marginal_means_conditional_list <- function(terms, parameters) {

  intercept_parameter <- parameters[terms == "intercept"]

  conditional_list <- lapply(seq_along(parameters), function(i) {

    c(
      if (length(intercept_parameter) > 0L && terms[i] != "intercept") {
        intercept_parameter
      },
      parameters[i]
    )
  })
  names(conditional_list) <- parameters

  return(conditional_list)
}


# Resolve whether marginal-means summaries should include inclusion BFs.
.marginal_means_resolve_bf <- function(model_averaged, bf = NULL, default = NULL) {

  if (is.null(default)) {
    default <- isTRUE(model_averaged)
  }

  if (is.null(bf)) {
    return(default)
  }

  BayesTools::check_bool(bf, "bf")

  return(bf)
}


# Remove BF columns and BF-only warnings from marginal-means summaries.
.marginal_means_drop_bf <- function(table) {

  table_type <- attr(table, "type")
  if (is.null(table_type) || !any(table_type == "inclusion_BF")) {
    return(table)
  }

  keep        <- table_type != "inclusion_BF"
  table_attrs <- attributes(table)
  table       <- table[, keep, drop = FALSE]

  class(table)        <- table_attrs[["class"]]
  attr(table, "type") <- table_type[keep]

  copy_attrs <- setdiff(
    names(table_attrs),
    c("names", "row.names", "class", "type")
  )
  for (attr_name in copy_attrs) {
    attr(table, attr_name) <- table_attrs[[attr_name]]
  }

  warnings <- attr(table, "warnings")
  if (!is.null(warnings)) {
    warnings <- warnings[!grepl("Savage-Dickey", warnings, fixed = TRUE)]
    if (length(warnings) == 0L) {
      warnings <- NULL
    }
    attr(table, "warnings") <- warnings
  }

  return(table)
}


# Resolve marginal-means summary/plot type.
.marginal_means_type <- function(object, type) {

  if (is.null(type)) {
    return("averaged")
  }

  BayesTools::check_char(type, "type")
  type <- match.arg(type, c("averaged", "conditional"))

  if (!isTRUE(object[["model_averaged"]])) {
    stop("The 'type' argument is available only for RoBMA marginal means.",
         call. = FALSE)
  }

  return(type)
}


# Select a marginal-means parameter from raw or internal names.
.marginal_means_select_parameter <- function(x, parameter) {

  BayesTools::check_char(parameter, "parameter", check_length = 1)

  term_map <- x[["term_map"]]
  lookup   <- parameter

  if (lookup == "mu") {
    lookup <- "intercept"
  }

  index <- match(lookup, term_map[["parameter"]])

  if (is.na(index)) {
    index <- match(lookup, term_map[["term"]])
  }

  if (is.na(index)) {
    jags_lookup <- BayesTools::JAGS_parameter_names(
      parameters        = lookup,
      formula_parameter = "mu"
    )
    index <- match(jags_lookup, term_map[["parameter"]])
  }

  if (is.na(index)) {
    stop("Unknown marginal means parameter '", parameter,
         "'. Available terms are: ", .marginal_means_available_terms(x),
         ".", call. = FALSE)
  }

  return(as.list(term_map[index, , drop = FALSE]))
}


# Add the moderator term as the legend title for base marginal-means plots.
.marginal_means_add_base_legend <- function(samples_parameter, term_label,
                                            dots) {

  if (identical(dots[["legend"]], FALSE)) {
    return(invisible(NULL))
  }

  level_names <- .marginal_means_level_names(samples_parameter)
  if (length(level_names) == 0L) {
    return(invisible(NULL))
  }

  graphics::legend(
    if (is.null(dots[["legend_position"]])) "topright" else dots[["legend_position"]],
    legend = level_names,
    title  = term_label,
    col    = .marginal_means_recycle_plot_value(dots[["col"]], length(level_names), "black"),
    lty    = .marginal_means_recycle_plot_value(dots[["lty"]], length(level_names), 1),
    lwd    = .marginal_means_recycle_plot_value(dots[["lwd"]], length(level_names), 1),
    bty    = "n"
  )

  return(invisible(NULL))
}


# Add the moderator term as the legend title for ggplot marginal-means plots.
.marginal_means_set_ggplot_legend_title <- function(plot, term_label) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(plot)
  }

  return(plot +
    ggplot2::guides(
      color    = ggplot2::guide_legend(title = term_label),
      linetype = ggplot2::guide_legend(title = term_label)
    ) +
    ggplot2::theme(legend.title = ggplot2::element_text()))
}


# Extract and format marginal-means legend level names.
.marginal_means_level_names <- function(samples_parameter) {

  level_names <- names(samples_parameter)
  if (is.null(level_names)) {
    return(character(0L))
  }

  dif_matches <- gregexpr("[dif:", level_names, fixed = TRUE)
  has_dif     <- vapply(dif_matches, function(x) x[1] != -1L, logical(1))
  multi_dif   <- vapply(
    dif_matches,
    function(x) x[1] != -1L && length(x) > 1L,
    logical(1)
  )

  if (any(multi_dif)) {
    level_names[multi_dif] <- gsub("__xXx__", ":", level_names[multi_dif], fixed = TRUE)
  }
  if (any(has_dif & !multi_dif)) {
    single_dif <- has_dif & !multi_dif
    level_names[single_dif] <- substr(
      level_names[single_dif],
      regexpr("[dif:", level_names[single_dif], fixed = TRUE) + 5L,
      regexpr("]", level_names[single_dif], fixed = TRUE) - 1L
    )
  }

  no_dif <- !has_dif
  if (any(no_dif & grepl("[", level_names, fixed = TRUE))) {
    bracket_names <- no_dif & grepl("[", level_names, fixed = TRUE)
    level_names[bracket_names] <- substr(
      level_names[bracket_names],
      regexpr("[", level_names[bracket_names], fixed = TRUE) + 1L,
      regexpr("]", level_names[bracket_names], fixed = TRUE) - 1L
    )
  }

  return(level_names)
}


# Recycle graphical legend values to one value per marginal level.
.marginal_means_recycle_plot_value <- function(value, n, default) {

  if (is.null(value)) {
    value <- default
  }
  if (length(value) == 1L) {
    value <- rep(value, n)
  }

  return(value)
}


# Format available marginal-means terms for error messages.
.marginal_means_available_terms <- function(x) {

  term_map <- x[["term_map"]]
  terms    <- unique(term_map[["term"]])

  return(paste0("'", terms, "'", collapse = ", "))
}


# Count prior components plotted for each marginal level.
.marginal_means_prior_component_counts <- function(samples_parameter) {

  if (!is.list(samples_parameter)) {
    samples_parameter <- list(samples_parameter)
  }

  counts <- vapply(samples_parameter, function(sample) {

    prior_density <- attr(sample, "prior_density")
    if (is.null(prior_density)) {
      return(1L)
    }

    count <- 0L
    if (!is.null(prior_density[["density"]])) {
      count <- count + 1L
    }
    if (!is.null(prior_density[["points"]]) &&
        NROW(prior_density[["points"]]) > 0L) {
      count <- count + 1L
    }

    return(max(count, 1L))
  }, integer(1))

  return(counts)
}


# Configure prior-line defaults for marginal means plots.
.set_dots_prior_marginal_means <- function(dots_prior, n_levels,
                                           prior_component_counts) {

  if (is.null(dots_prior)) {
    dots_prior <- list()
  }

  n_prior_levels <- sum(prior_component_counts)
  if (is.null(dots_prior[["col"]]) && n_levels == 1L) {
    dots_prior[["col"]] <- "black"
  } else if (is.null(dots_prior[["col"]]) && n_levels > 1L) {
    level_col <- grDevices::palette.colors(
      n       = n_levels + 1L,
      palette = "Okabe-Ito"
    )[-1L]
    dots_prior[["col"]] <- rep(level_col, prior_component_counts)
  } else if (length(dots_prior[["col"]]) == 1L) {
    dots_prior[["col"]] <- rep(dots_prior[["col"]], n_prior_levels)
  } else if (length(dots_prior[["col"]]) == n_levels &&
             n_prior_levels != n_levels) {
    dots_prior[["col"]] <- rep(dots_prior[["col"]], prior_component_counts)
  }
  if (is.null(dots_prior[["lty"]])) {
    dots_prior[["lty"]] <- rep(2, n_prior_levels)
  } else if (length(dots_prior[["lty"]]) == 1L) {
    dots_prior[["lty"]] <- rep(dots_prior[["lty"]], n_prior_levels)
  } else if (length(dots_prior[["lty"]]) == n_levels &&
             n_prior_levels != n_levels) {
    dots_prior[["lty"]] <- rep(dots_prior[["lty"]], prior_component_counts)
  }
  if (!is.null(dots_prior[["linetype"]]) &&
      length(dots_prior[["linetype"]]) == 1L) {
    dots_prior[["linetype"]] <- rep(dots_prior[["linetype"]], n_prior_levels)
  } else if (!is.null(dots_prior[["linetype"]]) &&
             length(dots_prior[["linetype"]]) == n_levels &&
             n_prior_levels != n_levels) {
    dots_prior[["linetype"]] <- rep(
      dots_prior[["linetype"]],
      prior_component_counts
    )
  }

  return(dots_prior)
}
