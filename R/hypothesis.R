# ============================================================================ #
# hypothesis.R
# ============================================================================ #
#
# User-facing hypothesis Bayes factors for brma objects. BayesTools owns the
# expression parser and BF algebra; RoBMA resolves model parameters and attaches
# likelihood-aware posterior ordinates when requested.
#
# ============================================================================ #


#' @title Hypothesis Bayes Factors
#'
#' @description Computes Bayes factors for scalar hypotheses written as
#' expressions, for example \code{"mu = 0"}, \code{"year > 0"}, or
#' \code{"mu_alloc[alternate] > mu_alloc[random]"}.
#'
#' @param object a fitted \code{brma}, \code{BMA}, \code{RoBMA}, or a posterior
#' object accepted by \code{BayesTools::hypothesis_BF()}.
#' @param ... additional arguments passed to methods.
#'
#' @return A \code{BayesTools_hypothesis_BF} BayesTools table.
#'
#' @examples \dontrun{
#' fit <- brma(
#'   yi      = c(0.12, 0.20, 0.05, 0.33, 0.18),
#'   sei     = c(0.10, 0.12, 0.09, 0.15, 0.11),
#'   measure = "GEN"
#' )
#'
#' # Average-effect hypotheses
#' hypothesis(fit, "mu != 0 vs mu = 0")
#' hypothesis(fit, "mu > 0 vs mu = 0")
#' hypothesis(fit, "mu > 0 vs mu < 0")
#'
#' # Heterogeneity hypothesis
#' hypothesis(fit, "tau = 0 vs tau > 0", component = "scale")
#'
#' dat <- data.frame(
#'   yi    = c(0.12, 0.20, 0.05, 0.33, 0.18, 0.41),
#'   sei   = c(0.10, 0.12, 0.09, 0.15, 0.11, 0.14),
#'   group = factor(c("A", "A", "A", "B", "B", "B"))
#' )
#' fit_mod <- brma(
#'   yi      = yi,
#'   sei     = sei,
#'   mods    = ~ group,
#'   data    = dat,
#'   measure = "GEN"
#' )
#'
#' # Marginal-mean hypothesis for a factor level
#' emm <- marginal_means(fit_mod, density_method = "qCMDE")
#' hypothesis(emm, "group[B] > 0 vs group[B] = 0", density_method = "qCMDE")
#' }
#'
#' @seealso [marginal_means()], [plot.brma()], [bridge_sampler.brma()]
#' @export
hypothesis <- function(object, ...) {

  UseMethod("hypothesis")
}


#' @rdname hypothesis
#' @export
hypothesis.default <- function(object, ...) {

  BayesTools::hypothesis_BF(posterior = object, ...)
}


#' @rdname hypothesis
#'
#' @param hypothesis character vector with scalar hypothesis statements.
#' @param component parameter component. Defaults to \code{"auto"}, which
#' infers the component when possible. Use \code{"mods"} (alias
#' \code{"location"}), \code{"scale"}, or \code{"random"} to disambiguate
#' terms used in multiple model components. The random component supports
#' interval and directional hypotheses for semantic standard deviation,
#' correlation, and allocation quantities. Point-null hypotheses are available
#' only when the induced prior and transformation define a coherent point
#' density. Publication-bias parameters are not supported.
#' @param standardized_coefficients whether moderator and scale coefficients
#' are tested on the standardized predictor scale. Defaults to \code{FALSE}.
#' @param conditional whether to use the conditional posterior for product-space
#' model-averaging objects. Defaults to \code{FALSE}.
#' @param logBF whether to display the Bayes factor on the log scale.
#' @param BF01 whether to display the inverse Bayes factor.
#' @param seed optional seed used by BayesTools for sampled prior quantities.
#' @param columns output columns passed to \code{BayesTools::hypothesis_BF()}.
#' \code{"default"} returns the compact BayesTools table with
#' \code{Alternative}, \code{Null}, \code{BF}, and \code{error\%(BF)}
#' columns. \code{"all"} also returns \code{prior}, \code{posterior}, and
#' \code{method} columns. The \code{prior} and \code{posterior} columns are
#' diagnostics rather than probabilities for all hypothesis types.
#' @param density_method posterior density method. \code{"KDE"} uses the
#' standard BayesTools kernel density estimate. \code{"normal"} uses the
#' BayesTools normal approximation for point-null hypotheses on fitted
#' \code{brma} objects. \code{"qCMDE"} and \code{"IWMDE"} attach RoBMA
#' likelihood-aware posterior ordinates for direct point-null hypotheses and
#' propagate their \code{BF_error} estimates printed as \code{error\%(BF)}. For
#' \code{marginal_means.brma} objects, \code{"normal"} is not supported, and
#' \code{"qCMDE"} and \code{"IWMDE"} compute missing point-null ordinates from
#' the stored source model. Matching is case-insensitive.
#' @param density_control named list of qCMDE/IWMDE tuning settings.
#' @param n_samples number of prior samples/grid points used by BayesTools for
#' deterministic marginal prior densities.
#'
#' @export
hypothesis.brma <- function(object, hypothesis,
                            component = "auto",
                            standardized_coefficients = FALSE,
                            conditional = FALSE,
                            logBF = FALSE, BF01 = FALSE, seed = NULL,
                            density_method = c("KDE", "normal", "qCMDE", "IWMDE"),
                            density_control = NULL,
                            n_samples = 10000,
                            columns = "default", ...) {

  if (is.null(object[["fit"]]) || length(object[["fit"]]) == 0L) {
    stop("'hypothesis' requires a fitted brma object.", call. = FALSE)
  }

  BayesTools::check_char(hypothesis, "hypothesis", check_length = 0, allow_NA = FALSE)
  BayesTools::check_bool(standardized_coefficients, "standardized_coefficients")
  BayesTools::check_bool(conditional, "conditional")
  BayesTools::check_bool(logBF, "logBF")
  BayesTools::check_bool(BF01, "BF01")
  BayesTools::check_real(seed, "seed", check_length = 1, allow_NULL = TRUE, allow_NA = FALSE)
  BayesTools::check_char(columns, "columns", check_length = 0, allow_NA = FALSE)
  BayesTools::check_int(n_samples, "n_samples", lower = 2)
  .warn_unused_dots(
    dots    = list(...),
    allowed = character(),
    caller  = "hypothesis()"
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
  if (conditional && !.is_RoBMA(object)) {
    stop("'conditional' hypotheses are available only for RoBMA objects.",
         call. = FALSE)
  }
  if (.density_method_uses_precomputed(density_method, allow_normal = TRUE)) {
    .check_iwmde_available(object, "qCMDE/IWMDE hypothesis()")
  }

  selected <- .hypothesis_brma_select_parameter(
    object     = object,
    hypothesis = hypothesis,
    component  = component
  )
  parameter <- selected[["parameter"]]
  parameter_label <- .hypothesis_brma_alias_label(
    aliases   = selected[["aliases"]],
    parameter = parameter
  )
  .hypothesis_brma_check_supported_component(selected[["component"]])
  hypothesis <- .hypothesis_brma_rewrite(
    hypothesis = hypothesis,
    aliases    = selected[["aliases"]],
    parameter  = parameter
  )

  if (identical(selected[["component"]], "random")) {
    return(.hypothesis_brma_random(
      object                    = object,
      parameter                 = parameter,
      hypothesis                = hypothesis,
      standardized_coefficients = standardized_coefficients,
      conditional               = conditional,
      logBF                     = logBF,
      BF01                      = BF01,
      seed                      = seed,
      density_method            = density_method,
      n_samples                 = n_samples,
      columns                   = columns
    ))
  }

  if (.hypothesis_brma_log_intercept_draws_available(
      object                    = object,
      selected                  = selected,
      standardized_coefficients = standardized_coefficients,
      conditional               = conditional,
      density_method            = density_method)) {
    draws <- .hypothesis_brma_log_intercept_draws(
      object    = object,
      parameter = parameter,
      n_samples = n_samples,
      seed      = seed
    )
    posterior <- draws[["posterior"]]
    if (.density_method_uses_precomputed(density_method, allow_normal = TRUE)) {
      posterior <- .hypothesis_brma_attach_iwmde_log_intercept(
        object               = object,
        posterior            = posterior,
        parameter            = parameter,
        hypothesis           = hypothesis,
        n_points             = density_control[["n_points"]],
        max_samples          = density_control[["max_samples"]],
        normalization_points = density_control[["normalization_points"]],
        normalization_prob   = density_control[["normalization_prob"]],
        density_method       = density_method
      )
    }
    out <- BayesTools::hypothesis_BF(
      posterior      = posterior,
      prior          = draws[["prior"]],
      hypothesis     = hypothesis,
      parameter      = parameter,
      logBF          = logBF,
      BF01           = BF01,
      seed           = seed,
      columns        = columns,
      density_method = if (.density_method_uses_precomputed(
          density_method, allow_normal = TRUE)) "precomputed" else density_method
    )
    if (.density_method_uses_precomputed(density_method, allow_normal = TRUE)) {
      out <- .hypothesis_brma_append_iwmde_warnings(out, posterior)
    }
    return(out)
  }

  sample_parameter <- .as_mixed_posteriors_parameters(object, parameter)
  samples <- .brma_as_mixed_posteriors(
    object           = object,
    parameters       = sample_parameter,
    conditional      = if (conditional) parameter else NULL,
    conditional_rule = "AND",
    transform_scaled = !standardized_coefficients,
    n_prior_samples  = n_samples
  )
  density_sample_parameter <- .plot_brma_density_sample_parameter(
    samples          = samples,
    parameter        = parameter,
    sample_parameter = sample_parameter
  )
  posterior <- BayesTools::marginal_posterior(
    samples       = samples,
    parameter     = density_sample_parameter,
    prior_samples = TRUE,
    use_formula   = FALSE,
    n_samples     = n_samples
  )

  if (.density_method_uses_precomputed(density_method, allow_normal = TRUE)) {
    posterior <- .hypothesis_brma_attach_iwmde(
      object                   = object,
      posterior                = posterior,
      parameter                = parameter,
      parameter_label          = parameter_label,
      hypothesis               = hypothesis,
      conditional              = if (conditional) parameter else NULL,
      n_points                 = density_control[["n_points"]],
      max_samples              = density_control[["max_samples"]],
      normalization_points     = density_control[["normalization_points"]],
      normalization_prob       = density_control[["normalization_prob"]],
      density_method           = density_method,
      n_samples                = n_samples
    )
  }

  out <- BayesTools::hypothesis_BF(
    posterior      = posterior,
    hypothesis     = hypothesis,
    parameter      = parameter,
    logBF          = logBF,
    BF01           = BF01,
    seed           = seed,
    columns        = columns,
    density_method = if (.density_method_uses_precomputed(density_method, allow_normal = TRUE)) {
      "precomputed"
    } else {
      density_method
    }
  )

  if (.density_method_uses_precomputed(density_method, allow_normal = TRUE)) {
    out <- .hypothesis_brma_append_iwmde_warnings(
      table     = out,
      posterior = posterior
    )
  }

  return(out)
}

.hypothesis_brma_random <- function(
    object, parameter, hypothesis, standardized_coefficients,
    conditional, logBF, BF01, seed, density_method, n_samples, columns) {

  if (conditional) {
    stop(
      "Conditional product-space hypotheses are not available for semantic ",
      "random-effect quantities.",
      call. = FALSE
    )
  }
  if (.density_method_uses_precomputed(density_method, allow_normal = TRUE)) {
    stop(
      "qCMDE/IWMDE hypotheses are not available for semantic random-effect ",
      "quantities. Use density_method = 'KDE' or 'normal'.",
      call. = FALSE
    )
  }

  posterior <- .brma_random_parameter_select(
    object                    = object,
    parameter                 = parameter,
    standardized_coefficients = standardized_coefficients
  )
  prior <- .brma_random_parameter_select(
    object                    = object,
    parameter                 = parameter,
    standardized_coefficients = standardized_coefficients,
    prior                     = TRUE,
    n_prior_samples           = n_samples,
    seed                      = seed
  )
  posterior_values <- as.numeric(posterior[["samples"]][, 1L])
  prior_values     <- as.numeric(prior[["samples"]][, 1L])

  point_refs <- BayesTools::hypothesis_parse_point_reference(
    hypothesis     = hypothesis,
    allow_compound = TRUE
  )
  if (nrow(point_refs) > 0L) {
    if (any(!point_refs[["direct"]])) {
      stop(
        "Point-null tests for random-effect quantities require a direct scalar ",
        "parameter reference.",
        call. = FALSE
      )
    }
    reason <- .brma_random_parameter_point_test_reason(
      posterior[["spec"]]
    )
    if (nzchar(reason)) {
      stop(reason, call. = FALSE)
    }

    support <- .brma_random_parameter_support(
      posterior[["spec"]],
      posterior[["prior"]],
      posterior[["source_prior"]]
    )
    values <- point_refs[["value"]]
    at_boundary <- (is.finite(support[1L]) & values <= support[1L]) |
      (is.finite(support[2L]) & values >= support[2L])
    if (any(at_boundary)) {
      stop(
        "Point-null Bayes factors at the support boundary are not available ",
        "for random-effect quantity '", posterior[["entry"]][["term"]], "'.",
        call. = FALSE
      )
    }
  }
  if (length(unique(prior_values[is.finite(prior_values)])) < 2L) {
    stop(
      "Hypothesis tests are not defined for fixed random-effect quantity '",
      posterior[["entry"]][["term"]], "'.",
      call. = FALSE
    )
  }

  BayesTools::hypothesis_BF(
    posterior      = posterior_values,
    prior          = prior_values,
    hypothesis     = hypothesis,
    parameter      = parameter,
    logBF          = logBF,
    BF01           = BF01,
    seed           = seed,
    columns        = columns,
    density_method = density_method
  )
}

.hypothesis_brma_log_intercept_draws_available <- function(
    object, selected, standardized_coefficients, conditional,
    density_method) {

  if (standardized_coefficients || conditional ||
      !identical(selected[["component"]], "scale")) {
    return(FALSE)
  }

  entry <- .brma_parameter_select_entry(
    object    = object,
    parameter = selected[["parameter"]],
    component = selected[["component"]]
  )
  formula_parameter <- entry[["formula_parameter"]]
  if (is.null(formula_parameter) || length(formula_parameter) != 1L ||
      is.na(formula_parameter) || !nzchar(formula_parameter)) {
    return(FALSE)
  }
  formula_scale     <- attr(object[["fit"]], "formula_scale", exact = TRUE)
  scale_info        <- formula_scale[[formula_parameter]]

  isTRUE(attr(scale_info, "log_intercept", exact = TRUE)) &&
    identical(entry[["term"]], "intercept")
}

.hypothesis_brma_log_intercept_draws <- function(object, parameter,
                                                  n_samples, seed) {

  samples <- .brma_as_mixed_posteriors(
    object           = object,
    parameters       = parameter,
    conditional      = NULL,
    conditional_rule = "AND",
    transform_scaled = TRUE,
    n_prior_samples  = n_samples
  )
  posterior <- samples[[parameter]]
  prior_density <- attr(samples, "prior_densities", exact = TRUE)[[parameter]]
  if (is.null(posterior) || is.null(prior_density)) {
    stop(
      "Could not reconstruct posterior and prior draws for transformed scale ",
      "intercept '", parameter, "'.",
      call. = FALSE
    )
  }
  attr(posterior, "parameter")     <- parameter
  attr(posterior, "prior_density") <- prior_density
  class(posterior) <- unique(c(
    class(posterior),
    "marginal_posterior.simple",
    "marginal_posterior"
  ))

  list(
    posterior = posterior,
    prior     = NULL
  )
}

#' @rdname hypothesis
#' @export
bf_hypothesis <- function(object, ...) {

  hypothesis(object, ...)
}


#' @rdname hypothesis
#' @export
BF_hypothesis <- function(object, ...) {

  hypothesis(object, ...)
}
