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
#' standard BayesTools kernel density estimate. \code{"normal"} uses the
#' BayesTools normal approximation for point-null Bayes factors and does not
#' attach a plotting density. \code{"qCMDE"} attaches RoBMA row-normalized
#' q-grid conditional densities. \code{"IWMDE"}
#' attaches Chen-style moment-matched IWMDE densities for plotting. RoBMA
#' stores separate precomputed posterior ordinates at \code{null_hypothesis};
#' these ordinates do not alter the plotting grid. Matching is
#' case-insensitive.
#' @param density_control named list of density-estimation settings. Supported
#' entries are \code{n_points}, \code{max_samples}, \code{display_grid},
#' \code{normalization_points}, and \code{normalization_prob}. The
#' normalization entries are only used with \code{density_method = "qCMDE"}.
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
                                density_method = c(
                                  "KDE", "normal", "qCMDE", "IWMDE"
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
    density_method = density_method,
    allow_normal   = TRUE
  )
  if (.density_method_uses_precomputed(density_method, allow_normal = TRUE) ||
      !is.null(density_control)) {
    density_control <- .density_control_normalize(
      density_method  = density_method,
      density_control = density_control,
      allow_normal    = TRUE
    )
  }
  model_averaged <- .is_RoBMA(object)
  bf             <- .marginal_means_resolve_bf(model_averaged, bf)
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
    normal_approximation = identical(density_method, "normal"),
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
    iwmde_signature        = .iwmde_fit_signature(object),
    bf                     = bf,
    name                   = .summary.brma_model_names(object)
  )

  class(output) <- c("marginal_means.brma", "marginal_means")

  if (.density_method_uses_precomputed(density_method, allow_normal = TRUE)) {
    output <- .marginal_means_attach_iwmde(
      object                  = object,
      marginal_means_object   = output,
      n_points                = density_control[["n_points"]],
      max_samples             = density_control[["max_samples"]],
      normalization_points    = density_control[["normalization_points"]],
      normalization_prob      = density_control[["normalization_prob"]],
      density_method          = density_method,
      display_grid            = density_control[["display_grid"]],
      null_hypothesis         = null_hypothesis
    )
  }

  return(output)
}


.marginal_means_attach_iwmde <- function(object, marginal_means_object,
                                         n_points, max_samples,
                                         normalization_points,
                                         normalization_prob, density_method,
                                         display_grid, null_hypothesis) {

  if (is.null(normalization_points)) {
    normalization_points <- max(50L, n_points)
  }

  context              <- .iwmde_context(object)
  diagnostic_cache     <- .iwmde_diagnostic_cache()
  diagnostics          <- list()
  ordinate_diagnostics <- list()
  method               <- .density_method_iwmde_estimator(density_method)
  include_values       <- NULL

  for (type in c("averaged", "conditional")) {
    specs <- .iwmde_marginal_means_specs(
      marginal_means_object = marginal_means_object,
      parameter             = NULL,
      type                  = type,
      levels                = NULL
    )

    diagnostics[[type]] <- lapply(specs, function(spec) {
      .iwmde_parameter_diagnostic(
        context              = context,
        parameter            = spec[["label"]],
        n_points             = n_points,
        max_samples          = max_samples,
        normalization_points = normalization_points,
        normalization_prob   = normalization_prob,
        method               = method,
        include_values       = include_values,
        parameter_spec       = spec,
        display_grid_method  = display_grid,
        diagnostic_cache     = diagnostic_cache
      )
    })
    names(diagnostics[[type]]) <- names(specs)

    marginal_means_object <- .marginal_means_attach_iwmde_type(
      marginal_means_object = marginal_means_object,
      type                  = type,
      specs                 = specs,
      diagnostics           = diagnostics[[type]],
      density_method        = density_method
    )
  }

  specs <- .iwmde_marginal_means_specs(
    marginal_means_object = marginal_means_object,
    parameter             = NULL,
    type                  = "conditional",
    levels                = NULL
  )
  ordinate_diagnostics[["conditional"]] <- lapply(specs, function(spec) {
    .iwmde_parameter_ordinate_diagnostic(
      context              = context,
      parameter            = spec[["label"]],
      values               = null_hypothesis,
      max_samples          = max_samples,
      normalization_points = normalization_points,
      normalization_prob   = normalization_prob,
      method               = method,
      parameter_spec       = spec,
      diagnostic_cache     = diagnostic_cache
    )
  })
  names(ordinate_diagnostics[["conditional"]]) <- names(specs)
  marginal_means_object <- .marginal_means_attach_iwmde_ordinate_type(
    marginal_means_object = marginal_means_object,
    type                  = "conditional",
    specs                 = specs,
    diagnostics           = ordinate_diagnostics[["conditional"]],
    density_method        = density_method
  )

  marginal_means_object[["iwmde_diagnostics"]] <- diagnostics
  marginal_means_object[["iwmde_ordinate_diagnostics"]] <- ordinate_diagnostics
  marginal_means_object[["iwmde_settings"]] <- list(
    n_points             = n_points,
    max_samples          = max_samples,
    normalization_points = normalization_points,
    normalization_prob   = normalization_prob,
    density_method       = density_method,
    method               = method,
    display_grid         = display_grid,
    include_values       = include_values,
    ordinate_values      = null_hypothesis
  )
  marginal_means_object[["inference"]] <- .marginal_means_refresh_iwmde_bf(
    inference            = marginal_means_object[["inference"]],
    parameters           = marginal_means_object[["parameters"]],
    null_hypothesis      = marginal_means_object[["null_hypothesis"]],
    density_method       = marginal_means_object[["density_method"]],
    require_precomputed  = TRUE
  )

  return(marginal_means_object)
}


.marginal_means_attach_iwmde_type <- function(marginal_means_object, type,
                                              specs, diagnostics,
                                              density_method) {

  for (name in names(specs)) {
    diagnostic <- diagnostics[[name]]
    if (!identical(diagnostic[["status"]], "ok")) {
      next
    }

    parameter <- specs[[name]][["parameter"]]
    level     <- specs[[name]][["level"]]
    samples   <- marginal_means_object[["inference"]][[type]][[parameter]]
    if (is.null(samples)) {
      next
    }

    if (is.list(samples) && !is.null(samples[[level]])) {
      attr(samples[[level]], "posterior_density") <- .iwmde_posterior_density_attribute(
        diagnostic     = diagnostic,
        density_method = density_method
      )
    } else {
      attr(samples, "posterior_density") <- .iwmde_posterior_density_attribute(
        diagnostic     = diagnostic,
        density_method = density_method
      )
    }
    marginal_means_object[["inference"]][[type]][[parameter]] <- samples
  }

  return(marginal_means_object)
}


.marginal_means_attach_iwmde_ordinate_type <- function(marginal_means_object, type,
                                                       specs, diagnostics,
                                                       density_method) {

  for (name in names(specs)) {
    diagnostic <- diagnostics[[name]]
    if (!identical(diagnostic[["status"]], "ok")) {
      next
    }

    parameter <- specs[[name]][["parameter"]]
    level     <- specs[[name]][["level"]]
    samples   <- marginal_means_object[["inference"]][[type]][[parameter]]
    if (is.null(samples)) {
      next
    }

    if (is.list(samples) && !is.null(samples[[level]])) {
      attr(samples[[level]], "posterior_ordinate") <- .iwmde_posterior_ordinate_attribute(
        diagnostic     = diagnostic,
        density_method = density_method
      )
    } else {
      attr(samples, "posterior_ordinate") <- .iwmde_posterior_ordinate_attribute(
        diagnostic     = diagnostic,
        density_method = density_method
      )
    }
    marginal_means_object[["inference"]][[type]][[parameter]] <- samples
  }

  return(marginal_means_object)
}


.marginal_means_refresh_iwmde_bf <- function(inference, parameters,
                                             null_hypothesis,
                                             density_method = "KDE",
                                             require_precomputed = FALSE) {

  density_method <- .density_method_normalize(
    density_method = density_method,
    allow_normal   = TRUE
  )
  if (identical(density_method, "normal") ||
      is.null(inference[["conditional"]]) ||
      is.null(inference[["inference"]])) {
    return(inference)
  }

  for (parameter in parameters) {
    posterior <- inference[["conditional"]][[parameter]]
    if (is.null(posterior)) {
      next
    }
    if (!.marginal_means_has_bf_posterior_ordinate(posterior)) {
      if (isTRUE(require_precomputed)) {
        inference[["inference"]][[parameter]] <- .marginal_means_unavailable_bf(
          posterior = posterior,
          warning   = paste0(
            "Precomputed posterior ordinate was unavailable or failed ",
            "diagnostics; Savage-Dickey Bayes factor is not reported."
          )
        )
      }
      next
    }
    posterior <- .marginal_means_bf_posterior(posterior)

    inference[["inference"]][[parameter]] <- BayesTools::Savage_Dickey_BF(
      posterior            = posterior,
      null_hypothesis      = null_hypothesis,
      normal_approximation = FALSE,
      silent               = TRUE,
      density_method       = "precomputed"
    )
  }

  return(inference)
}


.marginal_means_unavailable_bf <- function(posterior, warning) {

  if (!is.list(posterior)) {
    out <- NA_real_
    attr(out, "warnings") <- warning
    return(out)
  }

  out <- lapply(posterior, function(sample) {
    bf <- NA_real_
    attr(bf, "warnings") <- warning
    return(bf)
  })
  names(out) <- names(posterior)

  return(out)
}


.marginal_means_has_posterior_density <- function(samples) {

  if (!is.list(samples)) {
    return(!is.null(attr(samples, "posterior_density")))
  }

  return(any(vapply(samples, function(sample) {
    !is.null(attr(sample, "posterior_density"))
  }, logical(1))))
}


.marginal_means_has_bf_posterior_ordinate <- function(samples) {

  if (!is.list(samples)) {
    return(.iwmde_posterior_ordinate_supports_bf(
      attr(samples, "posterior_ordinate")
    ))
  }

  has_ordinate <- vapply(samples, function(sample) {
    .iwmde_posterior_ordinate_supports_bf(attr(sample, "posterior_ordinate"))
  }, logical(1))

  return(length(has_ordinate) > 0L && all(has_ordinate))
}


.marginal_means_bf_posterior <- function(samples) {

  if (!is.list(samples)) {
    if (!.iwmde_posterior_ordinate_supports_bf(attr(samples, "posterior_ordinate"))) {
      attr(samples, "posterior_ordinate") <- NULL
      attr(samples, "posterior_density") <- NULL
    }
    return(samples)
  }

  for (i in seq_along(samples)) {
    if (!.iwmde_posterior_ordinate_supports_bf(attr(samples[[i]], "posterior_ordinate"))) {
      attr(samples[[i]], "posterior_ordinate") <- NULL
      attr(samples[[i]], "posterior_density") <- NULL
    }
  }

  return(samples)
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
    require_precomputed  = !is.null(object[["density_method"]]) &&
      .density_method_uses_precomputed(
        density_method = object[["density_method"]],
        allow_normal   = TRUE
      )
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
                                     output_measure = NULL, transform = NULL, ...) {

  type <- .marginal_means_type(object = x, type = type)
  BayesTools::check_bool(prior, "prior")
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
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
    .marginal_means_has_posterior_density(samples[[selected[["parameter"]]]])
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
