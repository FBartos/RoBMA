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
#' @details For a continuous focal parameter \eqn{\theta}, nuisance parameters
#' \eqn{\xi}, and an unnormalized joint posterior kernel
#' \eqn{q(\theta,\xi)}, qCMDE estimates the posterior ordinate as
#' \deqn{\widehat{p}_{Q}(\theta^{*}\mid y)=\frac{1}{S}\sum_{i=1}^{S}
#' \frac{q(\theta^{*},\xi_{i})}{\int_{\Theta} q(t,\xi_{i})\,dt}.}
#' Each row is therefore normalized over the focal-parameter support before the
#' conditional densities are averaged. IWMDE instead uses
#' \deqn{\widehat{p}_{I}(\theta^{*}\mid y)=\frac{1}{S}\sum_{i=1}^{S}
#' w(\theta_{i}\mid\xi_{i})\frac{q(\theta^{*},\xi_{i})}
#' {q(\theta_{i},\xi_{i})},}
#' where \eqn{w(\theta\mid\xi)} is a normalized, fitted conditional weight.
#' Spike-mixture and product-space results additionally retain the posterior
#' mass of the continuous active branch. These estimators build on conditional
#' marginal density estimation and importance-weighted marginal density
#' estimation \insertCite{gelfand1992bayesian,chen1994importance}{RoBMA}.
#'
#' qCMDE is the preferred likelihood-aware method when it is supported and its
#' additional numerical normalization cost is acceptable. IWMDE can be faster,
#' but is more sensitive to the adequacy and concentration of its fitted
#' conditional weights. The normal approximation is useful as a rough check for
#' near-normal interior ordinates; it is not a validation reference for skewed
#' tails, bounded parameters, or one-sided support boundaries.
#' For binomial and Poisson GLMMs, the nuisance state includes the sampled
#' estimate-level heterogeneity effects and baserates or log-rates. Averaging
#' their row-conditional focal densities is an exact marginalization identity.
#' IWMDE is unavailable for GLMM density estimation because its
#' high-dimensional conditional weights did not meet bridge-sampling
#' certification. GLMM density curves and point-null hypotheses require qCMDE.
#'
#' For an unrestricted alternative against \eqn{\theta=\theta_0}, the
#' Savage-Dickey identity \insertCite{dickey1971weighted,wagenmakers2010bayesian}{RoBMA}
#' gives
#' \eqn{BF_{10}=p(\theta_0\mid H_1)/p(\theta_0\mid y,H_1)} only when the null
#' model's nuisance-parameter prior is the conditional prior induced by the
#' unrestricted model at \eqn{\theta_0}. Otherwise, fit the constrained model
#' separately and compare marginal likelihoods, for example by bridge sampling.
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
#' @references
#' \insertAllCited{}
#'
#' @seealso [density_diagnostics()], [marginal_means()], [plot.brma()],
#' [bridge_sampler.brma()]
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
#' @param hypothesis character vector with scalar hypothesis statements. For
#' \code{marginal_means.brma} objects, constants are specified on the fitted
#' linear-predictor scale, even when the object stores a display transformation.
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
#' propagate their \code{BF_error} estimates printed as \code{error\%(BF)}.
#' \code{BF_error} is a conditional Monte Carlo error estimate for the computed
#' posterior ordinate, not a complete standard error. It excludes uncertainty
#' in the prior ordinate and, for IWMDE, uncertainty from estimating the
#' conditional weight function. Use \code{density_diagnostics()} to inspect the
#' attached computation and reliability diagnostics. For
#' \code{marginal_means.brma} objects, \code{"normal"} is not supported. Their
#' default \code{NULL} reuses the density method stored by
#' \code{marginal_means()}; explicitly request \code{"KDE"} to override it.
#' \code{"qCMDE"} and \code{"IWMDE"} compute missing point-null ordinates from
#' the stored source model. Matching is case-insensitive. qCMDE/IWMDE ordinates
#' use the exact fitted-to-original structural coefficient map and induced
#' prior density stored by BayesTools. qCMDE/IWMDE supports exact linear maps.
#' For nonlinear joint maps, such as an exponentiated intercept that also
#' depends on varying slopes, use \code{density_method = "KDE"} or
#' \code{standardized_coefficients = TRUE}.
#' @param density_control named list of qCMDE/IWMDE tuning settings. Supported
#' entries are \code{n_points} (default \code{100}), \code{max_samples}
#' (default \code{Inf} for point ordinates), \code{initial_samples} (default
#' \code{500}), \code{target_relative_mcse} (default \code{0.05}),
#' \code{normalization_points} (default \code{NULL}, resolved to
#' \code{max(50, n_points)}), \code{normalization_prob} (default \code{0.999}),
#' and \code{display_grid} (default \code{"adaptive"}). Point ordinates use
#' deterministic increasing row budgets, beginning at \code{initial_samples},
#' until the precision target and every BF-grade reliability gate are met,
#' every eligible row is used, or the \code{max_samples} hard cap is reached.
#' \code{display_grid} is immaterial for
#' point-only requests. If all rows are used without meeting the target, obtain
#' more posterior draws before reporting a precise Bayes factor. Increase
#' \code{normalization_points} and \code{normalization_prob} to check numerical
#' support coverage. A point ordinate warns at relative MCSE at least 5 percent,
#' ESS below 100, maximum contribution share at least 20 percent, or fewer than
#' 100 finite contributions. It is rejected at relative MCSE at least 25
#' percent, ESS below 20, maximum share at least 50 percent, or fewer than 20
#' finite contributions. qCMDE ordinate movement warns above 2.5 percent and is
#' rejected above 5 percent; IWMDE normalization error warns above 5 percent and
#' is rejected above 10 percent. Adaptive-quadrature sensitivity warns above
#' 2.5 percent and is rejected above 5 percent.
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
      allow_normal    = TRUE,
      purpose         = "ordinate"
    )
  }
  if (conditional && !.is_RoBMA(object)) {
    stop("'conditional' hypotheses are available only for RoBMA objects.",
         call. = FALSE)
  }
  if (.density_method_uses_precomputed(density_method, allow_normal = TRUE)) {
    .check_iwmde_available(object, "qCMDE/IWMDE hypothesis()")
  }

  hypothesis <- BayesTools::hypothesis_parse(hypothesis)
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
  coefficient_target <- .hypothesis_brma_formula_coefficient_target(
    object                    = object,
    selected                  = selected,
    standardized_coefficients = standardized_coefficients
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

  sample_parameter <- .as_mixed_posteriors_parameters(object, parameter)
  samples <- .brma_as_mixed_posteriors(
    object           = object,
    parameters       = sample_parameter,
    conditional      = if (conditional) parameter else NULL,
    conditional_rule = "AND",
    transform_scaled = !standardized_coefficients,
    n_prior_samples  = n_samples
  )
  if (!is.null(coefficient_target)) {
    coefficient_target <- .hypothesis_brma_formula_prior_target(
      object      = object,
      samples     = samples,
      hypothesis  = hypothesis,
      target_info = coefficient_target
    )
    prior_densities <- attr(samples, "prior_densities", exact = TRUE)
    prior_densities[[parameter]] <- coefficient_target[["prior_density"]]
    attr(samples, "prior_densities") <- prior_densities
  }
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
  if (!is.null(coefficient_target)) {
    attr(posterior, "prior_density") <- coefficient_target[["prior_density"]]
  }

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
      initial_samples          = density_control[["initial_samples"]],
      target_relative_mcse     = density_control[["target_relative_mcse"]],
      normalization_points     = density_control[["normalization_points"]],
      normalization_prob       = density_control[["normalization_prob"]],
      density_method           = density_method,
      n_samples                = n_samples,
      parameter_spec           = if (is.null(coefficient_target)) {
        NULL
      } else {
        coefficient_target[["parameter_spec"]]
      }
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
    hypothesis     = BayesTools::hypothesis_render(hypothesis),
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

.hypothesis_brma_formula_coefficient_target <- function(
    object, selected, standardized_coefficients) {

  entry <- selected[["entry"]]
  if (standardized_coefficients ||
      selected[["component"]] %in% c("random", "bias") ||
      is.null(entry)) {
    return(NULL)
  }
  formula_parameter <- entry[["formula_parameter"]]
  if (is.null(formula_parameter) || length(formula_parameter) != 1L ||
      is.na(formula_parameter) || !nzchar(formula_parameter)) {
    return(NULL)
  }

  transform <- BayesTools::JAGS_formula_coefficient_transform(
    fit          = object[["fit"]],
    parameter    = formula_parameter,
    target_scale = "original"
  )
  if (!inherits(transform, "BayesTools_formula_coefficient_transform") ||
      !identical(transform[["schema_version"]], 1L)) {
    stop(
      "Formula coefficient transformation metadata are unsupported. Refit ",
      "the model with the current BayesTools version.",
      call. = FALSE
    )
  }
  target <- selected[["parameter"]]
  target_i <- match(target, transform[["target_names"]])
  if (is.na(target_i)) {
    stop(
      "Resolved formula coefficient '", target,
      "' is absent from the fitted coefficient transformation.",
      call. = FALSE
    )
  }

  return(list(
    formula_parameter = formula_parameter,
    target            = target,
    target_i          = target_i,
    transform         = transform
  ))
}


.hypothesis_brma_formula_prior_target <- function(
    object, samples, hypothesis, target_info) {

  prior_density <- BayesTools::JAGS_formula_prior_density(
    fit          = object[["fit"]],
    parameter    = target_info[["formula_parameter"]],
    target       = target_info[["target"]],
    target_scale = "original",
    context      = attr(samples, "prior_density_context", exact = TRUE)
  )
  refs <- .hypothesis_brma_point_refs(
    hypothesis,
    target_info[["target"]],
    require_direct = FALSE
  )
  for (value in refs[["value"]]) {
    ordinate <- BayesTools::prior_density_ordinate(prior_density, value)
    if (!isTRUE(ordinate[["exact"]])) {
      stop(
        "The induced prior ordinate for transformed coefficient '",
        target_info[["target"]], "' is not exact enough for a point-null ",
        "Bayes factor.",
        call. = FALSE
      )
    }
  }

  transform <- target_info[["transform"]]
  target_i  <- target_info[["target_i"]]
  weights   <- transform[["matrix"]][target_i, ]
  weights   <- weights[weights != 0]
  source_transforms <- transform[["source_transforms"]][names(weights)]
  output_transform  <- transform[["output_transforms"]][[
    target_info[["target"]]
  ]]
  identity <- all(source_transforms == "identity") &&
    identical(output_transform, "identity")
  unit_target <- length(weights) == 1L &&
    identical(names(weights), target_info[["target"]]) &&
    identical(unname(weights), 1)
  ordinary_unit_identity <- identity && unit_target
  log_unit_identity <- unit_target &&
    identical(source_transforms, stats::setNames("log", names(weights))) &&
    identical(output_transform, "exp")

  parameter_spec <- if (ordinary_unit_identity || log_unit_identity) {
    list(type = "primitive", prior_density = prior_density)
  } else if (identity) {
    list(type = "linear", weights = weights, prior_density = prior_density)
  } else {
    list(
      type   = "unsupported_formula_transform",
      reason = paste0(
        "qCMDE/IWMDE does not support the fitted nonlinear joint transform ",
        "for '", target_info[["target"]], "'. Use density_method = 'KDE' ",
        "or standardized_coefficients = TRUE."
      )
    )
  }

  target_info[["prior_density"]] <- prior_density
  target_info[["parameter_spec"]] <- parameter_spec
  return(target_info)
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
