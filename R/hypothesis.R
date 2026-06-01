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
#' \code{"location"}) or \code{"scale"} to disambiguate terms used in
#' multiple model components. Publication-bias parameters are not supported.
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

  sample_parameter <- .as_mixed_posteriors_parameters(object, parameter)
  samples <- BayesTools::as_mixed_posteriors(
    model            = object[["fit"]],
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
      standardized_coefficients = standardized_coefficients,
      n_samples                = n_samples
    )
  }

  BayesTools::hypothesis_BF(
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
