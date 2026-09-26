
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
#' for PET/PEESE regression plots. Factor terms select all coefficient cells;
#' use a semantic selector such as \code{"group[level]"} to plot one cell.
#' Structural point masses are retained.
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
#' case-insensitive. For semantic random-effect quantities, qCMDE/IWMDE support
#' direct scalar fitted sources, allocated component SDs backed by a scalar
#' aggregate (`tau_total` or `tau_common`), and allocation proportions or
#' multipliers
#' backed by a fitted simplex coordinate. Other nonlinear derived quantities
#' remain KDE-only. Shared inclusion gates are supported for aggregate and
#' component SDs and allocation proportions. SDs retain their excluded zero
#' branches; proportions condition on positive total heterogeneity.
#' In independently gated allocations, realized
#' `tau_total`/`tau2_total` and `tau2_prop(...)` are KDE-only because they combine
#' multiple gate and allocation coordinates. Their structural point masses are
#' displayed separately from the continuous density; `tau2_prop(...)` is
#' conditioned on positive realized total heterogeneity.
#' qCMDE/IWMDE are not available for non-known-\code{V}
#' \code{brma.mv()} random-formula models, \code{brma.mv()} models with
#' component-specific scale formulas (a named \code{scale} list), or
#' selection-weightfunction coordinates requiring joint replacement. IWMDE is
#' also unavailable for binomial and Poisson GLMMs; use
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
#' \code{normalization_prob} (default \code{0.999}; see [hypothesis()] for the
#' per-row truncation it permits).
#' \code{integration_control} (default \code{NULL}) retains the fitted
#' selection-integration settings. Supply a control created by
#' [set_selection_likelihood_control()] to change those settings for this
#' post-fit calculation, for example
#' \code{list(integration_control = set_selection_likelihood_control(max_points_per_scramble = 32768))}.
#' This entry is available for Gaussian selection models with a fitted integration plan.
#' The fitted object and posterior draws are unchanged. The maximum point budget
#' controls factor QMC fallback; analytic and deterministic quadrature rules
#' remain unchanged.
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
#' log-scale OR, RR, HR, and IRR models. For scale-regression coefficients,
#' \code{"EXP"} displays the multiplicative change in heterogeneity.
#' \code{"LOG"} displays a positive heterogeneity intercept on the log scale.
#' The transformation is applied to KDE and precomputed qCMDE/IWMDE densities
#' with the corresponding Jacobian.
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
#' @details A qCMDE/IWMDE request that the fitted model does not support stops
#' before any density is estimated with an error of class
#' `RoBMA_density_method_unavailable` and one class naming the cause:
#' `RoBMA_density_method_random_unknown_v` (`brma.mv()` random-formula models
#' without known `V`), `RoBMA_density_method_scale_components`
#' (component-specific scale formulas), or `RoBMA_density_method_glmm` (IWMDE
#' for binomial and Poisson GLMMs).
#'
#' A returned qCMDE/IWMDE estimate that fails the density availability
#' checks raises a `RoBMA_density_plot_error` identifying the plotted parameter.
#' Its `density_diagnostics` field contains the returned diagnostic records,
#' including computed density values and selected rows, without executable
#' estimate plans that could retain the fitted object.
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
  parameter_entry <- .brma_parameter_select_entry(x, parameter, allow_factor_cells = TRUE)
  plot_transform  <- .plot_output_setup(
    object          = x,
    parameter       = parameter,
    parameter_entry = parameter_entry,
    output_measure  = output_measure,
    transform       = transform
  )

  ### obtain posterior samples in the plotting format
  is_random       <- identical(parameter_entry[["component"]], "random")
  is_factor_cell  <- !is.null(parameter_entry[["parent_parameter"]])
  structural_cell <- FALSE
  random_label    <- NULL
  if (is_factor_cell) {
    cell <- .plot_brma_factor_cell_samples(
      object = x, entry = parameter_entry,
      standardized_coefficients = standardized_coefficients,
      conditional = conditional,
      precomputed = .density_method_uses_precomputed(density_method)
    )
    samples <- cell[["samples"]]
    sample_parameter <- density_sample_parameter <- parameter
    structural_cell <- cell[["structural"]]
    if (.density_method_uses_precomputed(density_method) && !structural_cell) {
      samples <- .plot_brma_attach_iwmde(
        object = x, samples = samples, parameter = parameter,
        sample_parameter = parameter,
        conditional = if (conditional) parameter_entry[["parent_parameter"]] else NULL,
        n_points = density_control[["n_points"]],
        sample_budget = density_control[["samples"]],
        normalization_points = density_control[["normalization_points"]],
        normalization_prob = density_control[["normalization_prob"]],
        integration_control = density_control[["integration_control"]],
        density_method = density_method,
        display_grid = density_control[["display_grid"]],
        parameter_spec = cell[["parameter_spec"]]
      )
    }
  } else if (is_random) {
    if (conditional && .density_method_uses_precomputed(density_method)) {
      stop("Conditional random-effect plots support 'density_method = \"KDE\"' only.",
           call. = FALSE)
    }
    sample_parameter <- parameter
    samples <- .brma_random_parameter_mixed_posterior(
      object                    = x,
      parameter                 = parameter,
      standardized_coefficients = standardized_coefficients,
      conditional               = conditional
    )
    density_sample_parameter <- parameter
    random_label <- attr(samples, "random_parameter_label", exact = TRUE)
    # A quantity without a prior density (e.g. an original-scale correlation
    # mixing SDs) is drawn without its prior curve by BayesTools, which warns
    # with class 'BayesTools_prior_curve_unavailable'.
    if (.density_method_uses_precomputed(density_method)) {
      target <- .brma_random_parameter_density_target(
        x,
        parameter,
        operation = "plots"
      )
      if (is.null(target[["parameter"]])) {
        stop(target[["reason"]], call. = FALSE)
      }
      samples <- .plot_brma_attach_iwmde(
        object               = x,
        samples              = samples,
        parameter            = target[["parameter"]],
        sample_parameter     = density_sample_parameter,
        conditional          = NULL,
        n_points             = density_control[["n_points"]],
        sample_budget        = density_control[["samples"]],
        normalization_points = density_control[["normalization_points"]],
        normalization_prob   = density_control[["normalization_prob"]],
        integration_control  = density_control[["integration_control"]],
        density_method       = density_method,
        display_grid         = density_control[["display_grid"]],
        parameter_spec       = target[["parameter_spec"]],
        display_transform    = target[["display_transform"]]
      )
    }
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
        integration_control     = density_control[["integration_control"]],
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
    dots[["par_name"]] <- if (is_random) random_label else if (is_factor_cell) {
      parameter_entry[["selection"]][["quantities"]][["display_label"]]
    } else .plot_parameter_label(
      parameter        = parameter,
      effect_transform = plot_transform,
      entry            = parameter_entry,
      object           = x
    )
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
    if (.density_method_uses_precomputed(density_method) && !structural_cell) {
      stop(.plot_brma_iwmde_unavailable_error(
        samples        = samples,
        density_method = density_method,
        parameter      = if (is_random) random_label else parameter
      ))
    }
    "KDE"
  }
  if (.effect_output_requested(plot_transform)) {
    args$transformation           <- .effect_plot_transformation(plot_transform)
    args$transformation_arguments <- NULL
    args$transformation_settings  <- TRUE
  }

  # suppress messages about transformations
  renderer <- BayesTools::plot_posterior
  if (is_factor_cell) {
    args[c("n_samples", "force_samples", "individual", "show_figures")] <- NULL
    renderer <- BayesTools::plot_marginal
  }
  plot <- suppressMessages(do.call(renderer, args))

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
                                    display_grid, parameter_spec = NULL,
                                    display_transform = NULL,
                                    integration_control = NULL) {

  if (is.null(normalization_points)) {
    normalization_points <- max(50L, n_points)
  }
  samples <- .plot_brma_clear_posterior_density(
    samples          = samples,
    sample_parameter = sample_parameter
  )
  context        <- .iwmde_context(object, integration_control)
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
      integration_control  = integration_control,
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
      integration_control  = integration_control,
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
    point_masses      <- diagnostic[["point_masses"]]
    if (!is.null(display_transform)) {
      posterior_density <- .plot_brma_transform_iwmde_density(
        posterior_density,
        display_transform
      )
      point_masses <- .plot_brma_transform_iwmde_point_masses(
        point_masses,
        display_transform
      )
    }
    if (!exact_parameter_spec) {
      posterior_density <- .plot_brma_align_iwmde_density(
        posterior_density = posterior_density,
        raw_samples       = diagnostic[["samples"]],
        plotted_samples   = plotted_samples
      )
    }
    if (!is.null(posterior_density)) {
      samples[[sample_parameter]] <- .iwmde_attach_posterior_density(
        samples[[sample_parameter]],
        posterior_density,
        point_masses = point_masses
      )
    }
  }

  return(samples)
}


.plot_brma_transform_iwmde_density <- function(
    posterior_density, display_transform) {

  if (is.null(posterior_density)) {
    return(NULL)
  }
  source_x <- posterior_density[["x"]]
  jacobian <- BayesTools::parameter_transform_jacobian(
    source_x,
    display_transform
  )
  if (any(!is.finite(jacobian) | jacobian <= 0)) {
    stop("Unsupported qCMDE/IWMDE display transform.", call. = FALSE)
  }
  posterior_density[["x"]] <- BayesTools::parameter_transform_forward(
    source_x,
    display_transform
  )
  posterior_density[["y"]] <- posterior_density[["y"]] / jacobian
  if (is.unsorted(posterior_density[["x"]])) {
    order <- order(posterior_density[["x"]])
    posterior_density[["x"]] <- posterior_density[["x"]][order]
    posterior_density[["y"]] <- posterior_density[["y"]][order]
  }
  if (length(posterior_density[["x"]]) > 1L &&
      any(diff(posterior_density[["x"]]) <= 0)) {
    stop("Unsupported qCMDE/IWMDE display transform.", call. = FALSE)
  }
  posterior_density
}


# The point masses of a qCMDE/IWMDE estimate (a table with 'x' and 'mass')
# at their display-scale locations.
.plot_brma_transform_iwmde_point_masses <- function(
    point_masses, display_transform) {

  if (is.data.frame(point_masses) && nrow(point_masses) > 0L &&
      "x" %in% names(point_masses)) {
    point_masses[["x"]] <- BayesTools::parameter_transform_forward(
      point_masses[["x"]],
      display_transform
    )
  }

  return(point_masses)
}


.plot_brma_formula_parameter_spec <- function(
    object, parameter, parameter_entry, standardized_coefficients) {

  if (standardized_coefficients) {
    return(NULL)
  }
  route <- .brma_formula_coefficient_route(
    object = object,
    selected = list(
      parameter = parameter,
      component = parameter_entry[["component"]],
      entry     = parameter_entry
    )
  )
  if (is.null(route)) {
    return(NULL)
  }
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
    BayesTools::posterior_metadata(
      samples[[sample_parameter]],
      "posterior_density"
    ) <- NULL
    BayesTools::posterior_metadata(
      samples[[sample_parameter]],
      "posterior_densities"
    ) <- NULL
  }

  return(samples)
}


.plot_brma_attach_iwmde_factor <- function(object, samples, parameter,
                                           sample_parameter, conditional,
                                           n_points, sample_budget,
                                           normalization_points,
                                           normalization_prob,
                                           density_method, display_grid,
                                           context, estimate_cache,
                                           integration_control = NULL) {

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

  # Each plotted column is the level cell of its label parts; its density is
  # stored under the column name, which BayesTools matches for the column.
  for (column_i in seq_len(ncol(plot_sample))) {
    column <- colnames(plot_sample)[[column_i]]
    level  <- .plot_brma_factor_column_level(
      sample            = plot_sample,
      column            = column,
      display_posterior = display_posterior
    )
    if (is.null(level) ||
        !level %in% names(raw_posterior) ||
        !level %in% names(display_posterior)) {
      next
    }

    weights        <- BayesTools::posterior_metadata(
      raw_posterior[[level]],
      "linear_weights"
    )
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
        integration_control  = integration_control,
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
        parameter = column,
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
      posterior_densities[[column]] <- posterior_density
      density_columns <- c(density_columns, column)
    }
  }

  if (length(diagnostics) == 0L && ncol(sample) == 1L) {
    column  <- colnames(sample)[[1L]]
    expected_columns <- column
    estimate <- .iwmde_estimate(
      context         = context,
      parameter       = parameter,
      density_method  = density_method,
      density_control = list(
        n_points             = n_points,
        samples              = sample_budget,
        normalization_points = normalization_points,
        normalization_prob   = normalization_prob,
        integration_control  = integration_control,
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
        parameter = column
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
        posterior_densities[[column]] <- posterior_density
        density_columns <- c(density_columns, column)
      }
    }
  }

  if (length(diagnostics) > 0L) {
    attr(samples, "iwmde_diagnostics") <- diagnostics
  }
  if (.plot_brma_factor_density_complete(density_columns, expected_columns)) {
    BayesTools::posterior_metadata(
      samples[[sample_parameter]],
      "posterior_density"
    ) <- NULL
    BayesTools::posterior_metadata(
      samples[[sample_parameter]],
      "posterior_densities"
    ) <- posterior_densities
  } else {
    BayesTools::posterior_metadata(
      samples[[sample_parameter]],
      "posterior_density"
    ) <- NULL
    BayesTools::posterior_metadata(
      samples[[sample_parameter]],
      "posterior_densities"
    ) <- NULL
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


# The marginal posterior of the level cell that a plotted factor column
# holds: the element whose label parts name the same level of every factor.
.plot_brma_factor_column_level <- function(sample, column, display_posterior) {

  quantities <- BayesTools::posterior_metadata(sample, "quantities")
  row        <- match(column, quantities[["column"]])
  if (is.na(row)) {
    return(NULL)
  }
  cell <- quantities[["label_parts"]][[row]][["levels"]]
  if (length(cell) == 0L) {
    return(NULL)
  }

  matches <- vapply(display_posterior, function(level_samples) {
    level_quantities <- BayesTools::posterior_metadata(
      level_samples,
      "quantities"
    )
    length(level_quantities[["label_parts"]]) == 1L &&
      identical(level_quantities[["label_parts"]][[1L]][["levels"]], cell)
  }, logical(1))
  if (sum(matches) != 1L) {
    return(NULL)
  }

  return(names(display_posterior)[matches])
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
    !is.null(BayesTools::posterior_metadata(
      samples[[sample_parameter]],
      "posterior_density"
    )) ||
      length(BayesTools::posterior_metadata(
        samples[[sample_parameter]],
        "posterior_densities"
      )) > 0L
  )
}


.plot_brma_iwmde_unavailable_error <- function(samples, density_method,
                                               parameter) {

  details <- .plot_brma_iwmde_unavailable_reason(samples)
  subject <- paste0(density_method, " density for '", parameter, "'")
  if (is.null(details)) {
    message <- paste0(subject, " was unavailable.")
  } else {
    status <- if (isTRUE(details[["rejected"]])) {
      " was rejected by diagnostics: "
    } else {
      " was unavailable: "
    }
    message <- paste0(subject, status, sub("[.]+$", "", details[["reason"]]), ".")
  }
  diagnostics <- lapply(
    attr(samples, "iwmde_diagnostics", exact = TRUE),
    function(diagnostic) {

      # Executable row-state plans contain closures that can retain the fit.
      diagnostic[["plan"]] <- NULL
      diagnostic
    }
  )

  return(structure(
    list(
      message             = message,
      call                = NULL,
      parameter           = parameter,
      density_method      = density_method,
      density_diagnostics = diagnostics
    ),
    class = c("RoBMA_density_plot_error", "error", "condition")
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
        return(list(reason = reason, rejected = FALSE))
      }
      next
    }
    reason <- .iwmde_diagnostics_density_failure_reason(
      diagnostic[["diagnostics"]]
    )
    if (!is.null(reason)) {
      return(list(reason = reason, rejected = TRUE))
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


# Name mixed factor columns by the fitted coordinates they hold. A column that
# is a catalog quantity identical to one coordinate (a direct level cell such
# as treatment 'mu_g[3]', which holds coordinate 'mu_g[2]', or a contrast
# coefficient 'mu_g{j}') takes that coordinate's name; other columns keep
# theirs. All columns are renamed at once, since a level label can equal the
# name of another level's coordinate.
# The fitted coordinates that mixed-posterior columns of a factor term hold:
# the sole dependency, with unit weight, of each column's draw metadata. Other
# columns (combinations of coordinates) name no coordinate.
.plot_brma_factor_cell_coordinate_columns <- function(samples) {

  quantities <- BayesTools::posterior_metadata(samples, "quantities")
  rows       <- match(colnames(samples), quantities[["column"]])
  if (anyNA(rows)) {
    stop("Selected factor-cell source coordinates are unavailable.", call. = FALSE)
  }
  mapped <- vapply(rows, function(row) {
    dependencies <- quantities[["dependencies"]][[row]]
    weights      <- quantities[["weights"]][[row]]
    if (length(dependencies) == 1L && identical(unname(as.numeric(weights)), 1)) {
      dependencies
    } else {
      NA_character_
    }
  }, character(1))
  if (anyDuplicated(mapped[!is.na(mapped)])) {
    stop("Selected factor-cell source coordinates are ambiguous.", call. = FALSE)
  }

  return(mapped)
}


# Prepare one semantic factor cell using its catalog extraction weights. The
# parent term supplies only the established prior and conditioning metadata.
.plot_brma_factor_cell_samples <- function(
    object, entry, standardized_coefficients, conditional, precomputed) {

  parent <- entry[["parent_parameter"]]
  selection <- entry[["selection"]]
  key <- selection[["quantities"]][["extraction_key"]][[1L]]
  weights <- switch(key[["type"]],
    coordinate = stats::setNames(1, key[["dependencies"]]),
    factor_level = stats::setNames(key[["weights"]], key[["dependencies"]]),
    NULL)
  if (is.null(weights)) {
    stop("Selected factor-cell extraction metadata are unavailable.", call. = FALSE)
  }
  if (length(weights)) weights <- .iwmde_linear_weights(weights)
  raw <- .brma_as_mixed_posteriors(
    object, parent, conditional = if (conditional) parent else NULL,
    transform_scaled = FALSE
  )
  raw_marginal <- BayesTools::marginal_posterior(
    raw, parent, prior_samples = TRUE, use_formula = FALSE
  )
  same_weights <- function(sample, target = weights) {

    candidate <- .iwmde_linear_weights(BayesTools::posterior_metadata(
      sample,
      "linear_weights"
    ))
    !is.null(candidate) && identical(unname(candidate[order(names(candidate))]),
      unname(target[order(names(target))])) && setequal(names(candidate), names(target))
  }
  matches <- which(vapply(raw_marginal, same_weights, logical(1L)))
  if (!length(matches)) {
    stop("Selected factor-cell prior metadata do not match its fitted extraction weights.",
      call. = FALSE)
  }
  # Multiple cells with the same exact weights represent the same quantity
  # (for example structural zero interaction cells); no posterior matching.
  level <- names(raw_marginal)[matches[[1L]]]
  displayed <- raw
  marginal <- raw_marginal
  if (!standardized_coefficients) {
    displayed <- .brma_as_mixed_posteriors(
      object, parent, conditional = if (conditional) parent else NULL,
      transform_scaled = TRUE
    )
    marginal <- BayesTools::marginal_posterior(
      displayed, parent, prior_samples = TRUE, use_formula = FALSE
    )
  }
  value <- marginal[[level]]
  if (!is.numeric(value) || is.null(BayesTools::posterior_metadata(
    value,
    "prior_density"
  ))) {
    stop("Selected factor-cell posterior and prior are unavailable.", call. = FALSE)
  }
  # The mixed columns hold the fitted coordinates that the cell's extraction
  # key combines (their draw metadata name them).
  coordinate_columns <- .plot_brma_factor_cell_coordinate_columns(
    displayed[[parent]]
  )
  if (!all(key[["dependencies"]] %in% coordinate_columns)) {
    stop("Selected factor-cell source coordinates are unavailable.", call. = FALSE)
  }
  source_samples <- as.matrix(displayed[[parent]])[
    ,
    !is.na(coordinate_columns),
    drop = FALSE
  ]
  colnames(source_samples) <- coordinate_columns[!is.na(coordinate_columns)]
  draws <- BayesTools::parameter_draws(object[["fit"]], selection,
    model_samples = source_samples)
  if (nrow(as.matrix(draws)) != length(value)) {
    stop("Selected factor-cell draws have inconsistent row metadata.", call. = FALSE)
  }
  value <- .brma_draws_with_values(value, as.numeric(as.matrix(draws)))
  attr(value, "parameter") <- entry[["parameter"]]
  attr(value, "level_name") <- level
  class(value) <- unique(c(class(value), "marginal_posterior"))
  structural <- !length(weights)
  if (precomputed && !standardized_coefficients && !structural) {
    transform <- BayesTools::JAGS_formula_coefficient_transform(
      object[["fit"]], entry[["formula_parameter"]], target_scale = "original"
    )
    targets <- match(names(weights), transform[["target_names"]])
    if (anyNA(targets) ||
        any(transform[["output_transforms"]][names(weights)] != "identity")) {
      stop("qCMDE/IWMDE for this factor cell requires 'standardized_coefficients = TRUE'.",
        call. = FALSE)
    }
    matrix <- transform[["matrix"]][targets, , drop = FALSE]
    weights <- stats::setNames(as.numeric(crossprod(weights, matrix)), colnames(matrix))
    weights <- .iwmde_linear_weights(weights)
    if (any(transform[["source_transforms"]][names(weights)] != "identity")) {
      stop("qCMDE/IWMDE for this factor cell requires 'standardized_coefficients = TRUE'.",
        call. = FALSE)
    }
  }
  samples <- stats::setNames(list(value), entry[["parameter"]])
  list(samples = samples, structural = structural,
    parameter_spec = .plot_brma_iwmde_parameter_spec(value,
      conditional = if (conditional) parent else NULL,
      type = "linear", weights = weights))
}
