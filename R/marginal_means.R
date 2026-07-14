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
#' @description S3 generic for estimated marginal means. The \code{brma}
#' method works with fitted moderator models and stores the BayesTools
#' marginal-inference object for \code{summary()} and \code{plot()}.
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
#' @param n_samples number of samples/grid points used by BayesTools for
#' marginal prior densities. Defaults to \code{10000}.
#' @param bf whether inclusion Bayes factors should be shown by default in
#' summaries. Defaults to \code{TRUE} for RoBMA/BMA objects and \code{FALSE}
#' for single-model \code{brma} objects.
#' @param density_method posterior density method. \code{"KDE"} uses the
#' standard BayesTools kernel density estimate. \code{"qCMDE"} attaches RoBMA
#' row-normalized q-grid conditional densities. \code{"IWMDE"} attaches
#' Chen-style moment-matched IWMDE densities for plotting. RoBMA stores
#' separate precomputed posterior ordinates at \code{null_hypothesis}; these
#' ordinates do not alter the plotting grid. qCMDE is preferred when its
#' additional normalization cost is acceptable; IWMDE can be faster but is
#' more sensitive to its fitted conditional weights. Matching is
#' case-insensitive. When an ordinate is used for a Bayes factor,
#' \code{BF_error} is a conditional Monte Carlo error estimate. It excludes
#' uncertainty in the prior ordinate and, for IWMDE, uncertainty from estimating
#' the conditional weight function. Use \code{density_diagnostics()} on a
#' subsequent \code{hypothesis()} result to inspect the ordinate diagnostics.
#' IWMDE is unavailable for binomial and Poisson GLMMs because its
#' high-dimensional conditional weights did not meet the bridge-sampling
#' certification tolerance. Use qCMDE for GLMM density curves and Bayes factors.
#' @param parameter,type,levels optional selectors restricting qCMDE/IWMDE
#' density precomputation to a marginal-means parameter, output type, or factor
#' levels. \code{parameter} and \code{levels} also restrict the conditional
#' ordinates used for inclusion Bayes factors. These ordinates are always
#' computed from conditional marginal means, even when \code{type} excludes
#' \code{"conditional"}.
#' @param density_control named list of density-estimation settings. Supported
#' entries are \code{n_points} (default \code{100}), \code{max_samples},
#' \code{initial_samples} (default \code{500}),
#' \code{target_relative_mcse} (default \code{0.05}), \code{display_grid}
#' (default \code{"adaptive"}), \code{normalization_points} (default
#' \code{NULL}, resolved to \code{max(50, n_points)}), and
#' \code{normalization_prob} (default \code{0.999}). \code{max_samples}
#' defaults to \code{500} for stored density curves and \code{Inf} for point
#' ordinates. Point ordinates adapt deterministically from
#' \code{initial_samples} until the relative-MCSE target and all BF-grade
#' reliability gates pass, all eligible rows are used, or a finite
#' \code{max_samples} cap is reached. The normalization entries are used with
#' \code{density_method = "qCMDE"} and \code{density_method = "IWMDE"}.
#' qCMDE/IWMDE are unavailable for non-known-\code{V} \code{brma.mv()}
#' random-formula models and for derived semantic random-effect quantities.
#' @inheritParams predict.brma
#' @param ... unused additional arguments. Supplied arguments trigger a warning.
#'
#' @return A list of class \code{marginal_means.brma} containing the
#' BayesTools \code{marginal_inference} object and parameter metadata.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE) &&
#'     requireNamespace("metafor", quietly = TRUE)) {
#'   data(dat.bcg, package = "metadat")
#'   dat <- metafor::escalc(
#'     measure = "RR",
#'     ai      = tpos,
#'     bi      = tneg,
#'     ci      = cpos,
#'     di      = cneg,
#'     data    = dat.bcg
#'   )
#'
#'   fit <- brma(yi = yi, vi = vi, mods = ~ alloc, data = dat, measure = "RR")
#'   mm  <- marginal_means(fit)
#'   summary(mm)
#'   plot(mm, parameter = "alloc")
#' }
#' }
#'
#' @seealso [summary()], [plot()], [summary.brma()], [regplot()]
#' @export
marginal_means.brma <- function(object, null_hypothesis = 0,
                                n_samples = 10000,
                                output_measure = NULL, transform = NULL,
                                bf = NULL,
                                parameter = NULL, type = NULL, levels = NULL,
                                density_method = c(
                                  "KDE", "qCMDE", "IWMDE"
                                ),
                                density_control = NULL, ...) {

  if (!.is_mods(object)) {
    stop("'marginal_means' requires a model with moderators.", call. = FALSE)
  }
  if (is.null(object[["fit"]]) || length(object[["fit"]]) == 0L) {
    stop("'marginal_means' requires a fitted brma object.", call. = FALSE)
  }

  BayesTools::check_real(null_hypothesis, "null_hypothesis", check_length = 1)
  BayesTools::check_int(n_samples, "n_samples", lower = 2)
  .warn_unused_dots(
    dots    = list(...),
    allowed = character(),
    caller  = "marginal_means()"
  )
  density_method <- .density_method_normalize(
    density_method = density_method
  )
  density_control_input <- density_control
  precompute_density <- .density_method_uses_precomputed(density_method)
  if (!precompute_density &&
      (!is.null(parameter) || !is.null(type) || !is.null(levels))) {
    stop(
      "'parameter', 'type', and 'levels' are only used when ",
      "'density_method' is 'qCMDE' or 'IWMDE'.",
      call. = FALSE
    )
  }
  if (precompute_density || !is.null(density_control)) {
    density_control <- .density_control_normalize(
      density_method  = density_method,
      density_control = density_control,
      purpose         = "density"
    )
  }
  ordinate_control <- if (precompute_density) {
    .density_control_normalize(
      density_method  = density_method,
      density_control = density_control_input,
      purpose         = "ordinate"
    )
  } else {
    NULL
  }
  precompute_type     <- .marginal_means_precompute_type(type)
  precompute_levels   <- .marginal_means_precompute_levels(levels)
  precompute_targeted <- !is.null(parameter) || !is.null(type) ||
    !is.null(levels)
  model_averaged           <- .is_RoBMA(object)
  bf                       <- .marginal_means_resolve_bf(model_averaged, bf)
  point_ordinate_supported <- !precompute_density ||
    .iwmde_point_ordinate_supported(object, density_method)
  if (precompute_density) {
    .iwmde_check_density_method_supported(object, density_method)
  }
  effect_transform <- .effect_output_setup(
    object         = object,
    output_measure = output_measure,
    transform      = transform
  )

  design           <- .fitted_formula_design(object, "mu", required = TRUE)
  formula          <- .fitted_formula_evaluable(object, design, "mods")
  terms            <- .marginal_means_terms(design)
  parameters       <- BayesTools::JAGS_parameter_names(
    parameters        = terms,
    formula_parameter = "mu"
  )
  parameter_setup  <- .marginal_means_drop_fixed_zero_intercept(
    model      = object[["fit"]],
    formula    = formula,
    terms      = terms,
    parameters = parameters
  )
  formula    <- parameter_setup[["formula"]]
  terms      <- parameter_setup[["terms"]]
  parameters <- parameter_setup[["parameters"]]
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
    normal_approximation = FALSE,
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
    n_samples              = n_samples,
    conditional_rule       = "OR",
    input_measure          = .measure(object),
    effect_transform       = effect_transform,
    model_averaged         = model_averaged,
    density_method         = density_method,
    point_ordinate_supported = point_ordinate_supported,
    bf                       = bf,
    name                     = .summary.brma_model_names(object),
    source_object            = object
  )

  class(output) <- c("marginal_means.brma", "marginal_means")

  if (precompute_density) {
    output <- .marginal_means_attach_iwmde(
      object                  = object,
      marginal_means_object   = output,
      n_points                = density_control[["n_points"]],
      max_samples             = density_control[["max_samples"]],
      normalization_points    = density_control[["normalization_points"]],
      normalization_prob      = density_control[["normalization_prob"]],
      density_method          = density_method,
      display_grid            = density_control[["display_grid"]],
      null_hypothesis         = null_hypothesis,
      parameter               = parameter,
      type                    = precompute_type,
      levels                  = precompute_levels,
      targeted                = precompute_targeted,
      include_ordinates       = point_ordinate_supported,
      ordinate_control        = ordinate_control
    )
  }

  return(output)
}


.marginal_means_drop_fixed_zero_intercept <- function(model, formula, terms,
                                                      parameters) {

  posterior_samples <- .get_posterior_samples(model)
  posterior_names   <- colnames(posterior_samples)
  prior_list        <- attr(model, "prior_list", exact = TRUE)
  intercept_index   <- which(parameters == "mu_intercept")

  if (length(intercept_index) != 1L ||
      .marginal_means_parameter_in_samples("mu_intercept", posterior_names)) {
    return(list(
      formula    = formula,
      terms      = terms,
      parameters = parameters
    ))
  }

  intercept_prior <- prior_list[["mu_intercept"]]
  if (is.null(intercept_prior) ||
      !BayesTools::is.prior.point(intercept_prior) ||
      !isTRUE(all.equal(as.numeric(mean(intercept_prior)), 0))) {
    return(list(
      formula    = formula,
      terms      = terms,
      parameters = parameters
    ))
  }

  keep <- seq_along(parameters) != intercept_index
  if (!any(keep)) {
    stop("No marginal means are available for this model.", call. = FALSE)
  }

  return(list(
    formula    = stats::update.formula(formula, . ~ . - 1),
    terms      = terms[keep],
    parameters = parameters[keep]
  ))
}


.marginal_means_parameter_in_samples <- function(parameter, posterior_names) {

  parameter %in% posterior_names ||
    any(startsWith(posterior_names, paste0(parameter, "[")))
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
                                        output_measure = NULL, transform = NULL, ...) {

  type <- .marginal_means_type(object = object, type = type)
  BayesTools::check_real(probs, "probs", allow_NULL = TRUE, check_length = 0)
  BayesTools::check_bool(logBF, "logBF")
  BayesTools::check_bool(BF01, "BF01")
  bf <- .marginal_means_resolve_bf(object[["model_averaged"]], bf, object[["bf"]])
  if (isTRUE(bf) && identical(object[["point_ordinate_supported"]], FALSE)) {
    .iwmde_check_point_ordinate_supported(
      object[["source_object"]],
      object[["density_method"]]
    )
  }

  if (missing(output_measure) && missing(transform)) {
    effect_transform <- object[["effect_transform"]]
  } else {
    effect_transform <- .effect_output_setup_measure(
      input_measure  = object[["input_measure"]],
      output_measure = if (missing(output_measure)) NULL else output_measure,
      transform      = if (missing(transform)) NULL else transform
    )
  }

  inference_object <- .marginal_means_refresh_iwmde_bf(
    inference            = object[["inference"]],
    parameters           = object[["parameters"]],
    null_hypothesis      = object[["null_hypothesis"]],
    density_method       = if (is.null(object[["density_method"]])) {
      "KDE"
    } else {
      object[["density_method"]]
    },
    object               = object
  )
  samples    <- .transform_marginal_samples_effect(
    samples          = inference_object[[type]],
    effect_transform = effect_transform
  )
  inference  <- inference_object[["inference"]]
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
#' @description Prints the estimated marginal means summary.
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
#' example \code{"measure"}, \code{"intercept"} for the intercept when
#' available, \code{"mu"} as an intercept alias, or the internal parameter name,
#' for example \code{"mu_measure"}.
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
#' @param density_method posterior density method. \code{"KDE"} uses the
#' standard BayesTools kernel density estimate and is the default, even when
#' qCMDE/IWMDE densities are stored on \code{x}. \code{"qCMDE"} and
#' \code{"IWMDE"} compute any missing densities for the plotted marginal means
#' before plotting and never silently mix with KDE. Plot-time qCMDE/IWMDE
#' computation can be slow and is used only for the current plot call; precompute
#' with \code{marginal_means(..., density_method = "qCMDE")} or
#' \code{marginal_means(..., density_method = "IWMDE")} to reuse densities across
#' plots. Prefer qCMDE for likelihood-aware final results; IWMDE can be faster
#' but is more sensitive to its fitted conditional weights.
#' @param density_control named list of qCMDE/IWMDE density-estimation settings.
#' Supported entries are \code{n_points} (default \code{100}),
#' \code{max_samples} (default \code{500}), \code{initial_samples} (default
#' \code{500}), \code{target_relative_mcse} (default \code{0.05}),
#' \code{display_grid} (default \code{"adaptive"}),
#' \code{normalization_points} (default \code{NULL}, resolved to
#' \code{max(50, n_points)}), and \code{normalization_prob} (default
#' \code{0.999}). \code{initial_samples} and \code{target_relative_mcse} are
#' point-ordinate controls and do not alter this fixed-budget plot. Supplying
#' \code{density_control} forces
#' recomputation for the plotted marginal means instead of reusing stored
#' densities. Increase row and normalization budgets when diagnostics report
#' concentrated contributions, dropped rows, or unstable normalization.
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
                                     output_measure = NULL, transform = NULL,
                                     density_method = c(
                                       "KDE", "qCMDE", "IWMDE"
                                     ),
                                     density_control = NULL, ...) {

  type <- .marginal_means_type(object = x, type = type)
  BayesTools::check_bool(prior, "prior")
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  density_control_requested <- !is.null(density_control)
  density_method <- .density_method_normalize(density_method)
  if (.density_method_uses_precomputed(density_method) ||
      !is.null(density_control)) {
    density_control <- .density_control_normalize(
      density_method  = density_method,
      density_control = density_control
    )
  }
  dots_raw <- list(...)
  .warn_unused_dots(
    dots    = dots_raw,
    allowed = .plot_dots_allowed(),
    caller  = "plot.marginal_means()"
  )
  dots_raw <- .keep_allowed_dots(dots_raw, .plot_dots_allowed())

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

  if (.density_method_uses_precomputed(density_method)) {
    x <- .marginal_means_plot_prepare_iwmde(
      x               = x,
      selected        = selected,
      type            = type,
      density_method  = density_method,
      density_control = density_control,
      force           = density_control_requested
    )
    samples <- x[["inference"]][[type]]
  }

  n_levels <- length(samples[[selected[["parameter"]]]])
  dots     <- do.call(.set_dots_plot, c(dots_raw, list(n_levels = n_levels)))
  if (is.null(dots[["xlab"]])) {
    dots[["xlab"]] <- .plot_parameter_label("mu", effect_transform)
  }
  if (is.null(dots[["legend_title"]])) {
    dots[["legend_title"]] <- selected[["label"]]
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
  args$density_method           <- if (
    .density_method_uses_precomputed(density_method)
  ) "precomputed" else "KDE"

  plot <- suppressMessages(do.call(BayesTools::plot_marginal, args))

  if (plot_type == "base") {
    return(invisible(plot))
  } else if (plot_type == "ggplot") {
    return(plot)
  }
}


# Extract model terms used by the moderator formula.
.marginal_means_terms <- function(formula) {

  if (inherits(formula, "BayesTools_formula_design")) {
    terms <- .formula_design_display_names(formula[["model_terms"]])
  } else {
    formula_terms <- stats::terms(formula)
    terms         <- c(
      if (attr(formula_terms, "intercept") == 1L) "intercept",
      attr(formula_terms, "term.labels")
    )
  }

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
    level_col <- .plot_level_palette(n_levels)
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
