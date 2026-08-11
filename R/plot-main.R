
### basic plotting functions ----
#' @title Plots brma Object
#'
#' @description \code{plot.brma} visualizes posterior
#' (and prior) distribution a brma object.
#'
#' @param x a fitted \code{brma}, \code{BMA}, or \code{RoBMA} object.
#' @param parameter a parameter to be plotted. Defaults to \code{"mu"} for
#' the effect size, or to the meta-regression intercept when moderators are
#' present. Additional options are \code{"tau"}, \code{"rho"} for multilevel
#' models, \code{"PET"}, \code{"PEESE"}, and \code{"omega"} or
#' \code{"weightfunction"} for selection models. Use \code{plot_pet_peese()}
#' for PET/PEESE regression plots.
#' @param parameter_mods legacy moderator selector. Prefer \code{parameter}
#' with \code{component = "mods"}. Use \code{"intercept"} for the
#' adjusted effect in meta-regression models.
#' @param parameter_scale legacy scale-regression selector. Prefer
#' \code{parameter} with \code{component = "scale"}. Use
#' \code{"intercept"} for the heterogeneity intercept in location-scale models.
#' @param component parameter component. Defaults to \code{"auto"}, which
#' infers the component when possible. Use \code{"mods"} (alias
#' \code{"location"}), \code{"scale"}, \code{"random"}, or \code{"bias"} to
#' disambiguate terms used in multiple model components. The random component
#' selects semantic standard deviation, correlation, and allocation parameters
#' from `brma.mv()` random formulas.
#' @param plot_type whether to use a base plot \code{"base"}
#' or ggplot2 \code{"ggplot"} for plotting. Defaults to
#' \code{"base"}.
#' @param prior whether prior distribution should be added to
#' figure. Defaults to \code{FALSE}.
#' @param standardized_coefficients whether to plot moderator and
#' scale-regression coefficients on the standardized predictor scale. Defaults
#' to \code{FALSE}.
#' @param conditional whether to plot the conditional posterior distribution
#' for RoBMA product-space objects. Defaults to \code{FALSE}.
#' @param density_method posterior density method. \code{"KDE"} uses the
#' standard BayesTools kernel density estimate. \code{"qCMDE"} attaches RoBMA
#' row-normalized q-grid conditional densities. \code{"IWMDE"} attaches
#' Chen-style moment-matched IWMDE densities. qCMDE is preferred when its
#' additional normalization cost is acceptable; IWMDE can be faster but is
#' more sensitive to its fitted conditional weights. Matching is
#' case-insensitive. qCMDE/IWMDE are not available for non-known-\code{V}
#' \code{brma.mv()} random-formula models, derived semantic random-effect
#' quantities, or selection-weightfunction coordinates requiring joint
#' replacement. IWMDE is also unavailable for binomial and Poisson GLMMs; use
#' qCMDE for GLMM density plots. qCMDE/IWMDE densities are evaluated on the
#' fitted coefficient coordinate. If automatic predictor scaling changes the
#' requested display coordinate, use `standardized_coefficients = TRUE`;
#' RoBMA does not infer a coefficient transformation from posterior draws.
#' @param density_control named list of density-estimation settings. Supported
#' entries are \code{n_points} (default \code{100}), \code{samples}
#' (default \code{500} for qCMDE and \code{1000} for IWMDE density curves),
#' \code{target_relative_mcse} (default \code{0.05}), \code{display_grid}
#' (default \code{"adaptive"}), \code{normalization_points} (default
#' \code{NULL}, resolved to \code{max(50, n_points)}), and
#' \code{normalization_prob} (default \code{0.999}).
#' \code{samples} controls the fixed posterior-row budget for the density
#' curve. \code{target_relative_mcse} is a point-ordinate diagnostic target and
#' does not alter this fixed-budget density plot. The normalization entries are
#' used with
#' \code{density_method = "qCMDE"} and \code{density_method = "IWMDE"}.
#' Curve diagnostics apply local relative-MCSE, effective-sample-size, and
#' contribution-concentration gates over the empirical 5--95 percent bulk and
#' record the 5 and 95 percent tail checkpoints. The entire display retains an
#' absolute MCSE safeguard relative to the density peak. Requested point
#' ordinates used for Bayes factors retain separate strict local diagnostics.
#' Increase the row and normalization budgets and compare results when density
#' diagnostics report low effective sample size, concentrated contributions,
#' or unstable normalization.
#' @param transform optional plotting transformation. \code{"EXP"} exponentiates
#' effect-size location and individual meta-regression coefficients for fitted
#' log-scale OR, RR, HR, and IRR models. \code{"LOG"} displays a positive
#' heterogeneity intercept on the log scale. The transformation is applied to
#' KDE and precomputed qCMDE/IWMDE densities with the corresponding Jacobian.
#' For automatically unscaled scale-formula intercepts, use
#' \code{standardized_coefficients = TRUE} with qCMDE/IWMDE.
#' @param dots_prior list of additional graphical arguments
#' to be passed to the plotting function of the prior
#' distribution. Supported arguments are \code{lwd},
#' \code{lty}, \code{col}, and \code{col.fill}, to adjust
#' the line thickness, line type, line color, and fill color
#' of the prior distribution respectively.
#' @inheritParams predict.brma
#' @param ... list of additional graphical arguments
#' to be passed to the plotting function. Supported arguments
#' are \code{lwd}, \code{lty}, \code{col}, \code{col.fill},
#' \code{xlab}, \code{ylab}, \code{ylab2}, \code{main}, \code{xlim},
#' \code{ylim}, and \code{ylim2}
#' to adjust the line thickness, line type, line color, fill color,
#' x-label, density-axis label, probability-axis label, title, x-axis range,
#' density-axis range, and probability-axis range respectively. Set
#' \code{ylim2} on the initial \code{plot()} call; subsequent base
#' \code{lines()} calls reuse that probability mapping.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'   fit <- bPET(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
#'
#'   plot(fit, parameter = "mu")
#'   lines(fit, parameter = "mu", col = "blue", lwd = 2)
#'
#'   ggplot_fit <- plot(fit, parameter = "mu", plot_type = "ggplot")
#'   ggplot_fit + lines(fit, parameter = "mu", plot_type = "ggplot", col = "blue")
#'
#'   plot(fit, parameter = "tau", prior = TRUE)
#'   plot(fit, parameter = "PET")
#' }
#' }
#'
#'
#' @return \code{plot.brma} returns either \code{NULL} if \code{plot_type = "base"}
#' or a \code{ggplot2} object if \code{plot_type = "ggplot"}.
#' \code{lines.brma} returns \code{NULL} for \code{plot_type = "base"} and
#' ggplot2 layer(s) for \code{plot_type = "ggplot"}.
#'
#' @seealso [RoBMA()]
#' @export
plot.brma <- function(
    x, parameter = NULL, parameter_mods = NULL, parameter_scale = NULL,
    prior = FALSE, standardized_coefficients = FALSE,
    conditional = FALSE,
    output_measure = NULL, transform = NULL,
    plot_type = "base", dots_prior = NULL,
    density_method = c("KDE", "qCMDE", "IWMDE"),
    density_control = NULL, component = "auto", ...) {

  .plot_brma(
    x                         = x,
    parameter                 = parameter,
    parameter_mods            = parameter_mods,
    parameter_scale           = parameter_scale,
    prior                     = prior,
    standardized_coefficients = standardized_coefficients,
    conditional               = conditional,
    output_measure            = output_measure,
    transform                 = transform,
    plot_type                 = plot_type,
    dots_prior                = dots_prior,
    density_method            = density_method,
    density_control           = density_control,
    component                 = component,
    add                       = FALSE,
    ...
  )
}

#' @details \code{lines.brma()} adds the posterior density to an existing base
#' plot. With \code{plot_type = "ggplot"}, it returns ggplot2 layer(s) that can
#' be added to a \code{plot.brma(..., plot_type = "ggplot")} object with
#' \code{+}. For base plots containing point masses, the initial
#' \code{plot.brma()} call establishes the secondary probability axis and
#' \code{lines.brma()} reuses it across fitted objects. An overlaid point mass
#' outside the initial \code{ylim2} is clipped with a warning.
#'
#' @rdname plot.brma
#' @export
lines.brma <- function(
    x, parameter = NULL, parameter_mods = NULL, parameter_scale = NULL,
    prior = FALSE, standardized_coefficients = FALSE,
    conditional = FALSE,
    output_measure = NULL, transform = NULL,
    plot_type = "base", dots_prior = NULL,
    density_method = c("KDE", "qCMDE", "IWMDE"),
    density_control = NULL, component = "auto", ...) {

  BayesTools::check_bool(prior, "prior")
  if (isTRUE(prior)) {
    stop(
      "'lines.brma' adds posterior densities only; use 'plot_prior()' or ",
      "prior-specific 'lines()' methods for prior overlays.",
      call. = FALSE
    )
  }

  .plot_brma(
    x                         = x,
    parameter                 = parameter,
    parameter_mods            = parameter_mods,
    parameter_scale           = parameter_scale,
    prior                     = FALSE,
    standardized_coefficients = standardized_coefficients,
    conditional               = conditional,
    output_measure            = output_measure,
    transform                 = transform,
    plot_type                 = plot_type,
    dots_prior                = dots_prior,
    density_method            = density_method,
    density_control           = density_control,
    component                 = component,
    add                       = TRUE,
    ...
  )
}

.plot_brma <- function(
    x, parameter = NULL, parameter_mods = NULL, parameter_scale = NULL,
    prior = FALSE, standardized_coefficients = FALSE,
    conditional = FALSE,
    output_measure = NULL, transform = NULL,
    plot_type = "base", dots_prior = NULL,
    density_method = c("KDE", "qCMDE", "IWMDE"),
    density_control = NULL, component = "auto", add = FALSE, ...) {

  ### check user input
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_bool(prior, "prior")
  BayesTools::check_bool(standardized_coefficients, "standardized_coefficients")
  BayesTools::check_bool(conditional, "conditional")
  BayesTools::check_bool(add, "add")
  dots_raw <- list(...)
  .warn_unused_dots(
    dots    = dots_raw,
    allowed = .plot_dots_allowed(),
    caller  = if (add) "lines.brma()" else "plot.brma()"
  )
  dots_raw <- .keep_allowed_dots(dots_raw, .plot_dots_allowed())
  density_method <- .density_method_normalize(density_method)
  if (.density_method_uses_precomputed(density_method)) {
    .iwmde_check_density_method_supported(x, density_method)
  }
  if (.density_method_uses_precomputed(density_method) ||
      !is.null(density_control)) {
    density_control <- .density_control_normalize(
      density_method  = density_method,
      density_control = density_control
    )
  }
  if (conditional && !.is_RoBMA(x)) {
    stop("'conditional' plots are available only for RoBMA objects.", call. = FALSE)
  }

  ### select and validate the parameter to be plotted
  parameter <- .check_and_select_plot_parameter(
    parameter        = parameter,
    parameter_mods   = parameter_mods,
    parameter_scale  = parameter_scale,
    component        = component,
    object           = x
  )
  parameter_entry <- .brma_parameter_select_entry(x, parameter)
  plot_transform  <- .plot_output_setup(
    object          = x,
    parameter       = parameter,
    parameter_entry = parameter_entry,
    output_measure  = output_measure,
    transform       = transform
  )

  ### obtain posterior samples in the plotting format
  is_random       <- identical(parameter_entry[["component"]], "random")
  random_label    <- NULL
  if (is_random) {
    if (conditional) {
      stop(
        "Conditional product-space plots are not available for semantic ",
        "random-effect quantities.",
        call. = FALSE
      )
    }
    if (.density_method_uses_precomputed(density_method)) {
      stop(
        "qCMDE/IWMDE plots are not available for derived random-effect quantities. ",
        "Use density_method = 'KDE'.",
        call. = FALSE
      )
    }
    sample_parameter <- parameter
    samples <- .brma_random_parameter_mixed_posterior(
      object                    = x,
      parameter                 = parameter,
      standardized_coefficients = standardized_coefficients,
      prior                     = prior
    )
    density_sample_parameter <- parameter
    random_label <- parameter_entry[["term"]]
  } else {
    sample_parameter <- .as_mixed_posteriors_parameters(x, parameter)
    samples <- .brma_as_mixed_posteriors(
      object           = x,
      parameters       = sample_parameter,
      conditional      = if (conditional) parameter else NULL,
      transform_scaled = !standardized_coefficients
    )
    density_sample_parameter <- .plot_brma_density_sample_parameter(
      samples          = samples,
      parameter        = parameter,
      sample_parameter = sample_parameter
    )
    if (.density_method_uses_precomputed(density_method)) {
      parameter_spec <- .plot_brma_formula_parameter_spec(
        object                    = x,
        parameter                 = parameter,
        parameter_entry           = parameter_entry,
        standardized_coefficients = standardized_coefficients
      )
      samples <- .plot_brma_attach_iwmde(
        object                  = x,
        samples                 = samples,
        parameter               = parameter,
        sample_parameter        = density_sample_parameter,
        conditional             = if (conditional) parameter else NULL,
        n_points                = density_control[["n_points"]],
        sample_budget           = density_control[["samples"]],
        normalization_points    = density_control[["normalization_points"]],
        normalization_prob      = density_control[["normalization_prob"]],
        density_method          = density_method,
        display_grid            = density_control[["display_grid"]],
        parameter_spec          = parameter_spec
      )
    }
  }

  ### set up plotting arguments
  n_levels   <- .get_samples_n_levels(samples, parameter)
  dots       <- do.call(.set_dots_plot, c(dots_raw, list(n_levels = n_levels)))
  dots_prior <- .set_dots_prior(dots_prior)
  if (is.null(dots[["par_name"]])) {
    dots[["par_name"]] <- if (is_random) random_label else
      .plot_parameter_label(parameter, plot_transform)
  }

  # prepare the argument call
  args                          <- dots
  args$samples                  <- samples
  args$parameter                <- parameter
  args$plot_type                <- plot_type
  args$prior                    <- prior
  args$n_points                 <- 1000
  args$n_samples                <- 10000
  args$force_samples            <- FALSE
  args$dots_prior               <- dots_prior
  args$individual               <- TRUE
  args$show_figures             <- NULL
  args$add                      <- add
  args$density_method           <- if (
    .plot_brma_has_posterior_density(samples, density_sample_parameter)
  ) {
    "precomputed"
  } else {
    if (.density_method_uses_precomputed(density_method)) {
      stop(
        .plot_brma_iwmde_unavailable_message(samples, density_method),
        call. = FALSE
      )
    }
    "KDE"
  }
  if (.effect_output_requested(plot_transform)) {
    args$transformation           <- .effect_plot_transformation(plot_transform)
    args$transformation_arguments <- NULL
    args$transformation_settings  <- FALSE
  }

  # suppress messages about transformations
  plot <- suppressMessages(do.call(BayesTools::plot_posterior, args))

  # return the plots
  if(plot_type == "base"){
    return(invisible(plot))
  }else if(plot_type == "ggplot"){
    return(plot)
  }
}


.plot_brma_attach_iwmde <- function(object, samples, parameter, sample_parameter,
                                    conditional,
                                    n_points, sample_budget,
                                    normalization_points,
                                    normalization_prob, density_method,
                                    display_grid, parameter_spec = NULL) {

  if (is.null(normalization_points)) {
    normalization_points <- max(50L, n_points)
  }
  samples <- .plot_brma_clear_posterior_density(
    samples          = samples,
    sample_parameter = sample_parameter
  )
  context        <- .iwmde_context(object)
  estimate_cache <- .iwmde_estimate_cache()

  if (inherits(samples[[sample_parameter]], "mixed_posteriors.factor")) {
    return(.plot_brma_attach_iwmde_factor(
      object               = object,
      samples              = samples,
      parameter            = parameter,
      sample_parameter     = sample_parameter,
      conditional          = conditional,
      n_points             = n_points,
      sample_budget        = sample_budget,
      normalization_points = normalization_points,
      normalization_prob   = normalization_prob,
      density_method       = density_method,
      display_grid         = display_grid,
      context              = context,
      estimate_cache       = estimate_cache
    ))
  }

  plotted_samples <- .plot_brma_plotted_samples(
    samples          = samples,
    sample_parameter = sample_parameter,
    parameter        = parameter
  )
  if (is.null(plotted_samples)) {
    return(samples)
  }
  exact_parameter_spec <- !is.null(parameter_spec)
  if (!exact_parameter_spec) {
    parameter_spec <- list(type = "primitive")
  }
  parameter_spec[["conditional"]]      <- conditional
  parameter_spec[["conditional_rule"]] <- "AND"

  estimate <- .iwmde_estimate(
    context         = context,
    parameter       = parameter,
    density_method  = density_method,
    density_control = list(
      n_points             = n_points,
      samples              = sample_budget,
      normalization_points = normalization_points,
      normalization_prob   = normalization_prob,
      display_grid         = display_grid
    ),
    outputs        = "density",
    parameter_spec = parameter_spec,
    metadata       = .iwmde_posterior_metadata(
      samples   = samples[[sample_parameter]],
      parameter = parameter
    ),
    cache          = estimate_cache
  )
  diagnostic <- estimate[["diagnostics"]][["density"]]
  attr(samples, "iwmde_diagnostics") <- list(parameter = diagnostic)

  if (identical(diagnostic[["status"]], "ok") &&
      !is.null(estimate[["posterior_density"]])) {
    posterior_density <- estimate[["posterior_density"]]
    if (!exact_parameter_spec) {
      posterior_density <- .plot_brma_align_iwmde_density(
        posterior_density = posterior_density,
        raw_samples       = diagnostic[["samples"]],
        plotted_samples   = plotted_samples
      )
    }
    if (!is.null(posterior_density)) {
      attr(samples[[sample_parameter]], "posterior_density") <-
        posterior_density
    }
  }

  return(samples)
}


.plot_brma_formula_parameter_spec <- function(
    object, parameter, parameter_entry, standardized_coefficients) {

  if (standardized_coefficients) {
    return(NULL)
  }
  target <- .hypothesis_brma_formula_coefficient_target(
    object = object,
    selected = list(
      parameter = parameter,
      component = parameter_entry[["component"]],
      entry     = parameter_entry
    )
  )
  if (is.null(target)) {
    return(NULL)
  }
  route <- .hypothesis_brma_formula_transform_route(target)
  if (identical(route[["type"]], "identity")) {
    return(list(type = "primitive"))
  }
  if (identical(route[["type"]], "affine")) {
    return(list(type = "linear", weights = route[["weights"]]))
  }
  if (!is.null(route[["reason"]])) {
    stop(route[["reason"]], call. = FALSE)
  }
  stop(
    "qCMDE/IWMDE does not support the fitted nonlinear joint transform for '",
    parameter, "'. Use density_method = 'KDE' or ",
    "standardized_coefficients = TRUE.",
    call. = FALSE
  )
}


.plot_brma_clear_posterior_density <- function(samples, sample_parameter) {

  if (!is.null(samples[[sample_parameter]])) {
    attr(samples[[sample_parameter]], "posterior_density")   <- NULL
    attr(samples[[sample_parameter]], "posterior_densities") <- NULL
  }

  return(samples)
}


.plot_brma_attach_iwmde_factor <- function(object, samples, parameter,
                                           sample_parameter, conditional,
                                           n_points, sample_budget,
                                           normalization_points,
                                           normalization_prob,
                                           density_method, display_grid,
                                           context, estimate_cache) {

  sample <- samples[[sample_parameter]]
  if (is.null(colnames(sample))) {
    return(samples)
  }
  plot_samples <- BayesTools::transform_factor_samples(samples)
  plot_sample  <- plot_samples[[sample_parameter]]
  if (is.null(plot_sample) || is.null(colnames(plot_sample))) {
    return(samples)
  }

  raw_samples <- .brma_as_mixed_posteriors(
    object           = object,
    parameters       = sample_parameter,
    conditional      = conditional,
    transform_scaled = FALSE
  )
  display_posterior <- BayesTools::marginal_posterior(
    samples       = samples,
    parameter     = sample_parameter,
    prior_samples = TRUE,
    use_formula   = FALSE
  )
  raw_posterior <- BayesTools::marginal_posterior(
    samples       = raw_samples,
    parameter     = sample_parameter,
    prior_samples = TRUE,
    use_formula   = FALSE
  )
  if (!is.list(display_posterior) || !is.list(raw_posterior)) {
    return(samples)
  }

  posterior_densities <- list()
  diagnostics         <- list()
  expected_columns    <- colnames(plot_sample)
  density_columns     <- character()

  for (column_i in seq_len(ncol(plot_sample))) {
    column  <- colnames(plot_sample)[[column_i]]
    aliases <- .plot_brma_factor_density_aliases(
      parameter   = sample_parameter,
      sample      = plot_sample,
      sample_name = column,
      level_i     = column_i
    )
    level  <- .plot_brma_factor_column_level(
      parameter   = sample_parameter,
      sample      = plot_sample,
      column      = column,
      column_i    = column_i,
      level_names = names(display_posterior)
    )
    if (is.null(level) ||
        !level %in% names(raw_posterior) ||
        !level %in% names(display_posterior)) {
      next
    }

    weights        <- attr(raw_posterior[[level]], "linear_weights", exact = TRUE)
    linear_weights <- .iwmde_linear_weights(weights)
    if (is.null(linear_weights) || length(linear_weights) == 0L) {
      next
    }

    estimate <- .iwmde_estimate(
      context         = context,
      parameter       = column,
      density_method  = density_method,
      density_control = list(
        n_points             = n_points,
        samples              = sample_budget,
        normalization_points = normalization_points,
        normalization_prob   = normalization_prob,
        display_grid         = display_grid
      ),
      outputs        = "density",
      parameter_spec = .plot_brma_iwmde_parameter_spec(
        samples     = raw_posterior[[level]],
        conditional = conditional,
        type        = "linear",
        weights     = weights
      ),
      metadata       = .iwmde_posterior_metadata(
        samples   = raw_posterior[[level]],
        parameter = aliases,
        level     = level
      ),
      cache          = estimate_cache
    )
    diagnostic <- estimate[["diagnostics"]][["density"]]
    diagnostics[[column]] <- diagnostic

    if (!identical(diagnostic[["status"]], "ok") ||
        is.null(estimate[["posterior_density"]])) {
      next
    }

    posterior_density <- estimate[["posterior_density"]]
    posterior_density <- .plot_brma_align_iwmde_density(
      posterior_density = posterior_density,
      raw_samples       = diagnostic[["samples"]],
      plotted_samples   = as.numeric(plot_sample[, column_i])
    )
    if (!is.null(posterior_density)) {
      posterior_densities <- .plot_brma_add_posterior_density(
        posterior_densities = posterior_densities,
        aliases             = aliases,
        posterior_density   = posterior_density
      )
      density_columns <- c(density_columns, column)
    }
  }

  if (length(diagnostics) == 0L && ncol(sample) == 1L) {
    column  <- colnames(sample)[[1L]]
    expected_columns <- column
    aliases <- .plot_brma_factor_density_aliases(
      parameter   = sample_parameter,
      sample      = sample,
      sample_name = column,
      level_i     = 1L
    )
    estimate <- .iwmde_estimate(
      context         = context,
      parameter       = parameter,
      density_method  = density_method,
      density_control = list(
        n_points             = n_points,
        samples              = sample_budget,
        normalization_points = normalization_points,
        normalization_prob   = normalization_prob,
        display_grid         = display_grid
      ),
      outputs        = "density",
      parameter_spec = .plot_brma_iwmde_parameter_spec(
        samples     = sample,
        conditional = conditional,
        type        = "primitive"
      ),
      metadata       = .iwmde_posterior_metadata(
        samples   = sample,
        parameter = aliases
      ),
      cache          = estimate_cache
    )
    diagnostic <- estimate[["diagnostics"]][["density"]]
    diagnostics[[column]] <- diagnostic

    if (identical(diagnostic[["status"]], "ok") &&
        !is.null(estimate[["posterior_density"]])) {
      posterior_density <- estimate[["posterior_density"]]
      posterior_density <- .plot_brma_align_iwmde_density(
        posterior_density = posterior_density,
        raw_samples       = diagnostic[["samples"]],
        plotted_samples   = as.numeric(sample[, 1L])
      )
      if (!is.null(posterior_density)) {
        posterior_densities <- .plot_brma_add_posterior_density(
          posterior_densities = posterior_densities,
          aliases             = aliases,
          posterior_density   = posterior_density
        )
        density_columns <- c(density_columns, column)
      }
    }
  }

  if (length(diagnostics) > 0L) {
    attr(samples, "iwmde_diagnostics") <- diagnostics
  }
  if (.plot_brma_factor_density_complete(density_columns, expected_columns)) {
    attr(samples[[sample_parameter]], "posterior_density")   <- NULL
    attr(samples[[sample_parameter]], "posterior_densities") <- posterior_densities
  } else {
    attr(samples[[sample_parameter]], "posterior_density")   <- NULL
    attr(samples[[sample_parameter]], "posterior_densities") <- NULL
  }

  return(samples)
}


.plot_brma_iwmde_parameter_spec <- function(samples, conditional, type,
                                           weights = NULL) {

  condition_metadata <- .iwmde_sample_condition_metadata(samples)
  if (!is.null(conditional)) {
    condition_metadata[["conditional"]]      <- conditional
    condition_metadata[["conditional_rule"]] <- "AND"
  }
  if (is.null(condition_metadata[["conditional_rule"]])) {
    condition_metadata[["conditional_rule"]] <- "AND"
  }

  spec <- c(
    list(
      type    = type,
      weights = weights
    ),
    condition_metadata
  )
  spec <- spec[!vapply(spec, is.null, logical(1))]

  return(spec)
}


.plot_brma_add_posterior_density <- function(posterior_densities, aliases,
                                             posterior_density) {

  aliases <- unique(aliases[!is.na(aliases) & nzchar(aliases)])
  for (alias in aliases) {
    posterior_densities[[alias]] <- posterior_density
  }

  return(posterior_densities)
}


.plot_brma_factor_density_complete <- function(density_columns, expected_columns) {

  density_columns <- unique(as.character(density_columns))
  density_columns <- density_columns[!is.na(density_columns) & nzchar(density_columns)]

  expected_columns <- unique(as.character(expected_columns))
  expected_columns <- expected_columns[!is.na(expected_columns) & nzchar(expected_columns)]

  if (length(expected_columns) == 0L) {
    return(FALSE)
  }

  return(setequal(density_columns, expected_columns))
}


.plot_brma_factor_column_level <- function(parameter, sample, column, column_i,
                                           level_names) {

  aliases <- .plot_brma_factor_density_aliases(
    parameter   = parameter,
    sample      = sample,
    sample_name = column,
    level_i     = column_i
  )
  bracket_matches <- regmatches(
    column,
    gregexpr("\\[[^]]+\\]", column)
  )[[1]]
  bracket_aliases <- gsub("^\\[|\\]$", "", bracket_matches)
  sample_level_names <- attr(sample, "level_names", exact = TRUE)
  if (is.list(sample_level_names) &&
      length(bracket_aliases) == length(sample_level_names)) {
    aliases <- c(
      aliases,
      .plot_brma_factor_named_cell_alias(
        level_names = sample_level_names,
        values      = bracket_aliases
      )
    )
    aliases <- c(
      aliases,
      .plot_brma_factor_named_cell_alias(
        level_names = sample_level_names,
        values      = .plot_brma_factor_display_aliases(bracket_aliases)
      )
    )
  }

  level <- intersect(level_names, aliases)
  if (length(level) != 1L) {
    return(NULL)
  }

  return(level)
}


.plot_brma_factor_density_aliases <- function(parameter, sample, sample_name,
                                             level_i) {

  aliases <- sample_name
  bracket_matches <- regmatches(
    sample_name,
    gregexpr("\\[[^]]+\\]", sample_name)
  )[[1]]
  bracket_aliases <- gsub("^\\[|\\]$", "", bracket_matches)
  if (length(bracket_aliases) > 0L) {
    aliases <- c(
      aliases,
      bracket_aliases,
      paste0(parameter, "[", bracket_aliases, "]")
    )
  }
  if (length(bracket_aliases) > 1L) {
    cell_alias <- paste0(bracket_aliases, collapse = ", ")
    aliases <- c(aliases, cell_alias, paste0(parameter, "[", cell_alias, "]"))
  }

  level_names <- attr(sample, "level_names", exact = TRUE)
  if (is.list(level_names)) {
    level_names <- .plot_brma_factor_cell_labels(level_names)
  }
  if (length(level_names) == ncol(sample)) {
    aliases <- c(
      aliases,
      level_names[[level_i]],
      paste0(parameter, "[", level_names[[level_i]], "]")
    )
  }

  factor_cell_names <- attr(sample, "factor_cell_names", exact = TRUE)
  if (length(factor_cell_names) == ncol(sample)) {
    aliases <- c(
      aliases,
      factor_cell_names[[level_i]],
      paste0(parameter, "[", factor_cell_names[[level_i]], "]")
    )
  }

  aliases <- unique(as.character(aliases))
  aliases <- aliases[!is.na(aliases) & nzchar(aliases)]

  return(aliases)
}


.plot_brma_factor_cell_labels <- function(level_names) {

  if (!is.list(level_names)) {
    return(level_names)
  }

  cells <- expand.grid(
    level_names,
    KEEP.OUT.ATTRS   = FALSE,
    stringsAsFactors = FALSE
  )
  labels <- apply(cells, 1L, function(x) {
    paste0(names(level_names), "=", as.character(x), collapse = ", ")
  })

  return(labels)
}


.plot_brma_factor_named_cell_alias <- function(level_names, values) {

  values <- as.character(values)
  if (length(values) != length(level_names)) {
    return(character())
  }

  return(paste0(names(level_names), "=", values, collapse = ", "))
}


.plot_brma_factor_display_aliases <- function(values) {

  values <- as.character(values)
  values <- sub("^dif: ", "", values)

  return(values)
}


.plot_brma_density_sample_parameter <- function(samples, parameter,
                                                sample_parameter) {

  if (parameter %in% c("PET", "PEESE", "omega") &&
      "bias" %in% sample_parameter &&
      !is.null(samples[["bias"]])) {
    return("bias")
  }
  if (!is.null(samples[[parameter]])) {
    return(parameter)
  }

  return(sample_parameter[[1L]])
}


.plot_brma_plotted_samples <- function(samples, sample_parameter, parameter) {

  sample <- samples[[sample_parameter]]
  if (is.null(sample) || is.list(sample)) {
    return(NULL)
  }
  if (is.matrix(sample)) {
    columns <- colnames(sample)
    if (is.null(columns)) {
      return(NULL)
    }
    column <- match(parameter, columns)
    if (is.na(column)) {
      return(NULL)
    }
    return(as.numeric(sample[, column]))
  }
  if (!is.numeric(sample)) {
    return(NULL)
  }

  return(as.numeric(sample))
}


.plot_brma_has_posterior_density <- function(samples, sample_parameter) {

  if (is.null(samples[[sample_parameter]])) {
    return(FALSE)
  }

  return(
    !is.null(attr(samples[[sample_parameter]], "posterior_density")) ||
      length(attr(samples[[sample_parameter]], "posterior_densities")) > 0L
  )
}


.plot_brma_iwmde_unavailable_message <- function(samples, density_method) {

  reason <- .plot_brma_iwmde_unavailable_reason(samples)
  if (is.null(reason)) {
    return(paste0(density_method, " density was not available."))
  }

  return(paste0(
    density_method,
    " density was rejected by diagnostics: ",
    reason,
    "."
  ))
}


.plot_brma_iwmde_unavailable_reason <- function(samples) {

  diagnostics <- attr(samples, "iwmde_diagnostics", exact = TRUE)
  if (is.null(diagnostics) || length(diagnostics) == 0L) {
    return(NULL)
  }

  for (diagnostic in diagnostics) {
    if (!identical(diagnostic[["status"]], "ok")) {
      reason <- diagnostic[["reason"]]
      if (length(reason) == 1L && !is.na(reason) && nzchar(reason)) {
        return(reason)
      }
      next
    }
    reason <- .iwmde_diagnostics_density_failure_reason(
      diagnostic[["diagnostics"]]
    )
    if (!is.null(reason)) {
      return(reason)
    }
  }

  return(NULL)
}


.plot_brma_align_iwmde_density <- function(posterior_density, raw_samples,
                                            plotted_samples) {

  if (is.null(posterior_density) ||
      !.plot_brma_same_sample_scale(raw_samples, plotted_samples)) {
    return(NULL)
  }

  return(posterior_density)
}


.plot_brma_same_sample_scale <- function(raw_samples, plotted_samples) {

  raw_samples     <- as.numeric(raw_samples)
  plotted_samples <- as.numeric(plotted_samples)
  if (length(raw_samples) != length(plotted_samples) ||
      length(raw_samples) == 0L ||
      any(!is.finite(raw_samples)) ||
      any(!is.finite(plotted_samples))) {
    return(FALSE)
  }

  return(identical(raw_samples, plotted_samples))
}
