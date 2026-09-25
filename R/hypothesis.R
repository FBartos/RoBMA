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
#' qCMDE is the default for fitted \code{brma} objects. Use \code{"KDE"},
#' \code{"normal"}, or \code{"IWMDE"} explicitly to select another route.
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
#' hypothesis(emm, "group[B] > 0 vs group[B] = 0")
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
#' Separate statements may reference different model parameters; each statement
#' itself must resolve to one scalar model target.
#' @param component parameter component. Defaults to \code{"auto"}, which
#' infers the component when possible. Use \code{"mods"} (alias
#' \code{"location"}), \code{"scale"}, or \code{"random"} to disambiguate
#' terms used in multiple model components. The random component supports
#' interval and directional hypotheses for semantic standard deviation,
#' correlation, and allocation quantities. Point-null hypotheses require a
#' direct parameter or level reference. Certified \code{exp(affine)}
#' fitted-scale hypotheses are available with KDE only for atom-free,
#' unconditional scalar targets. Realized gated allocation totals and
#' proportions can contain structural point masses and therefore support only
#' interval and directional hypotheses; variance proportions condition on
#' positive realized total heterogeneity. Publication-bias parameters are not
#' supported.
#' @param standardized_coefficients whether moderator and scale coefficients
#' are tested on the standardized predictor scale. Defaults to \code{FALSE}.
#' @param conditional whether to use the conditional posterior for product-space
#' model-averaging objects. Defaults to \code{FALSE}. When this argument is
#' omitted and the selected parameter has both null and alternative components,
#' a warning notes that the full ensemble is used. Pass \code{FALSE} explicitly
#' to retain that test without the warning, or \code{TRUE} to test only models
#' where the parameter is active.
#' @param logBF whether to display the Bayes factor on the log scale.
#' @param BF01 whether to display the inverse Bayes factor.
#' @param seed optional seed used by BayesTools for sampled prior quantities.
#' @param columns output columns passed to \code{BayesTools::hypothesis_BF()}.
#' \code{"default"} returns the compact BayesTools table with
#' \code{Alternative}, \code{Null}, \code{BF}, and \code{error\%(BF)}
#' columns. \code{"all"} also returns \code{prior}, \code{posterior}, and
#' \code{method} columns. The \code{prior} and \code{posterior} columns are
#' diagnostics rather than probabilities for all hypothesis types.
#' @param density_method posterior density method. Defaults to \code{"qCMDE"}
#' for fitted \code{brma} objects. \code{"KDE"} uses the standard BayesTools
#' kernel density estimate. \code{"normal"} uses the
#' BayesTools normal approximation for point-null hypotheses on fitted
#' \code{brma} objects. \code{"qCMDE"} and \code{"IWMDE"} attach RoBMA
#' likelihood-aware posterior ordinates for direct point-null hypotheses and
#' propagate their \code{BF_error} estimates printed as \code{error\%(BF)}.
#' \code{BF_error} is a Monte Carlo diagnostic for the computed posterior
#' ordinate, not a complete standard error. For a subsampled mixed posterior it
#' is a worst-correlation delta upper bound combining the selected active-row
#' and complete indicator-chain components. It excludes uncertainty in the
#' prior ordinate and, for IWMDE, uncertainty from estimating the conditional
#' weight function. At a location that is also a posterior atom, qCMDE/IWMDE
#' reports the continuous ordinate; the atom probability is a separate discrete
#' mass. Use \code{density_diagnostics()} to inspect the attached computation and
#' reliability diagnostics. For
#' \code{marginal_means.brma} objects, \code{"normal"} is not supported. Their
#' default \code{NULL} reuses the density method stored by
#' \code{marginal_means()}; explicitly request \code{"KDE"} to override it.
#' \code{"qCMDE"} and \code{"IWMDE"} compute missing point-null ordinates from
#' the stored source model. Matching is case-insensitive. qCMDE/IWMDE ordinates
#' use the exact fitted-to-original structural coefficient map and induced
#' prior density stored by BayesTools. qCMDE/IWMDE supports exact linear maps.
#' For nonlinear joint maps, such as an exponentiated intercept that also
#' depends on varying slopes, KDE point and directional hypotheses are available
#' only when structural metadata certifies an atom-free, unconditional scalar
#' target and the point is inside the open transformed support. Alternatively,
#' use \code{standardized_coefficients = TRUE}. Atom-free pairwise contrasts
#' across levels of one single-model factor are supported. Certified
#' \code{exp(affine)} point equalities are evaluated on
#' the inverse log/affine scale, where the prior and posterior Jacobians cancel;
#' calls mixing point and region statements must be evaluated separately.
#' @param density_control named list of qCMDE/IWMDE tuning settings. Supported
#' entries are \code{n_points} (default \code{100}), \code{samples} (the fixed
#' posterior-row sample size, default \code{500} for qCMDE and \code{1000} for
#' IWMDE; use \code{Inf} for the eligible-row census),
#' \code{target_relative_mcse} (default \code{0.05}),
#' \code{normalization_points} (default \code{NULL}, resolved to
#' \code{max(50, n_points)}), \code{normalization_prob} (default \code{0.999};
#' lower than 1 for qCMDE),
#' and \code{display_grid} (default \code{"adaptive"}).
#' \code{integration_control} (default \code{NULL}) retains the fitted
#' selection-integration settings. Supply a control created by
#' [set_selection_likelihood_control()] to change those settings for this
#' post-fit calculation, for example
#' \code{list(integration_control = set_selection_likelihood_control(max_points_per_scramble = 32768))}.
#' This entry is available for Gaussian selection models with a fitted integration plan.
#' The fitted object and posterior draws are unchanged. The maximum point budget
#' controls factor QMC fallback; analytic and deterministic quadrature rules
#' remain unchanged. Point ordinates use one
#' state-independent simple random sample selected before ordinate contributions
#' are evaluated. Multiple direct scalar point ordinates for one target share
#' the same conditional-normalization pass while retaining separate diagnostics.
#' A finite estimate is returned independently of its sample diagnostics. If the
#' fixed sample does not meet the precision target, a
#' warning recommends increasing \code{samples} or using the census; if the
#' census does not meet the target, obtain more posterior draws.
#' \code{display_grid} is immaterial for point-only requests. qCMDE normalizes
#' each posterior row's conditional density over a range covering at least its
#' central \code{normalization_prob} mass, on \code{normalization_points} nodes
#' and the nested grid of those nodes and their midpoints. The mass of a row
#' left outside the range, its truncation \eqn{t}, is therefore at most
#' \code{1 - normalization_prob}: exactly computed for rows whose Gaussian
#' likelihood kernel meets a normal prior, bounded for a Gaussian kernel with
#' another proper scalar prior, and estimated from the tails (each tail
#' extended until its estimate is at most \code{(1 - normalization_prob) / 2})
#' otherwise. A row normalizer missing the fraction \eqn{t} overstates that
#' row's density by the factor \eqn{1 / (1 - t)}, so the ordinate is
#' overstated by at most \eqn{t / (1 - t)} for the largest row truncation,
#' which is at most \code{(1 - normalization_prob) / normalization_prob}; at
#' the default \code{0.999} that is about 0.1 percent, and it is attained when
#' all rows share one conditional law. [density_diagnostics()] report the
#' largest row truncation (\code{normalization_truncation}), how it is known
#' (\code{normalization_truncation_status}), and this bound
#' (\code{truncation_ordinate_bound}), together with the discretization change
#' between the two nested grids (\code{ordinate_relative_change}). Set
#' \code{normalization_prob} closer to 1 to reduce the truncation, and increase
#' \code{normalization_points} to reduce the discretization. IWMDE checks its
#' normalization mass on
#' \code{normalization_points} points over the central \code{normalization_prob}
#' quantile range of the target draws, widened by 10 percent on each side.
#' A point ordinate warns at relative MCSE at least 5 percent,
#' ESS below 100, maximum contribution share at least 20 percent, or fewer than
#' 100 finite contributions. These same-sample diagnostics do not suppress a
#' finite fixed-design estimate. The qCMDE ordinate error bound (the ordinate
#' change between the nested grids plus the tail-truncation bound) warns above
#' 2.5 percent and is
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
                            density_method = "qCMDE",
                            density_control = NULL,
                            n_samples = 10000,
                            columns = "default", ...) {

  conditional_omitted <- missing(conditional)

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
  parameter_metadata <- .brma_parameter_catalog_metadata(object)
  hypothesis <- .hypothesis_brma_ast(
    hypothesis = hypothesis,
    catalog    = parameter_metadata[["catalog"]]
  )
  requested_point_refs <- BayesTools::hypothesis_parse_point_reference(
    hypothesis     = hypothesis,
    allow_compound = TRUE
  )
  if (.density_method_uses_precomputed(
      density_method,
      allow_normal = TRUE
    ) && nrow(requested_point_refs) > 0L) {
    .check_iwmde_available(object, "qCMDE/IWMDE hypothesis()")
  }
  display_hypothesis <- hypothesis
  statement_selections <- .hypothesis_brma_select_statements(
    object     = object,
    hypothesis = hypothesis,
    component  = component,
    metadata   = parameter_metadata
  )
  statement_keys <- vapply(statement_selections, function(selection) {
    paste(selection[["component"]], selection[["parameter"]], sep = "\r")
  }, character(1))
  if (length(unique(statement_keys)) > 1L) {
    return(.hypothesis_brma_multiple_parameters(
      object                    = object,
      hypothesis                = display_hypothesis,
      selections                = statement_selections,
      standardized_coefficients = standardized_coefficients,
      conditional               = conditional,
      conditional_omitted       = conditional_omitted,
      logBF                     = logBF,
      BF01                      = BF01,
      seed                      = seed,
      density_method            = density_method,
      density_control           = density_control,
      n_samples                 = n_samples,
      columns                   = columns
    ))
  }
  selected <- .hypothesis_brma_select_parameter(
    object     = object,
    hypothesis = hypothesis,
    component  = component,
    metadata   = parameter_metadata
  )
  parameter <- selected[["parameter"]]
  parameter_label <- .hypothesis_brma_alias_label(
    aliases   = selected[["aliases"]],
    parameter = parameter
  )
  .hypothesis_brma_check_supported_component(selected[["component"]])
  prior_list <- attr(object[["fit"]], "prior_list", exact = TRUE)
  if (conditional_omitted && .is_RoBMA(object) &&
      parameter %in% names(prior_list) &&
      BayesTools::is.prior.mixture(prior_list[[parameter]])) {
    prior_components <- attr(
      prior_list[[parameter]],
      "components",
      exact = TRUE
    )
    if (any(prior_components == "null") &&
        any(prior_components == "alternative")) {
      warning(
        "Model-averaged coefficient test: this hypothesis uses the full ",
        "ensemble for ", parameter, ", including models where ",
        parameter, " is fixed by its null component. To test the ",
        "hypothesis only within models where ", parameter,
        " is active, use conditional = TRUE.",
        call. = FALSE
      )
    }
  }
  hypothesis <- .hypothesis_brma_rewrite(
    hypothesis = hypothesis,
    aliases    = selected[["aliases"]],
    parameter  = parameter
  )
  level_contrast <- .hypothesis_brma_level_contrast_candidate(
    hypothesis = hypothesis,
    parameter  = parameter
  )
  point_refs <- .hypothesis_brma_point_refs(
    hypothesis     = hypothesis,
    parameter      = parameter,
    require_direct = !level_contrast
  )
  # Point statements on a level fixed by the contrast need no posterior
  # ordinate: the level's declared atom decides them.
  fixed_point_levels <- .hypothesis_brma_fixed_point_levels(
    object     = object,
    selected   = selected,
    point_refs = point_refs
  )
  only_fixed_points <- nrow(point_refs) > 0L &&
    all(!is.na(point_refs[["level"]]) &
          point_refs[["level"]] %in% fixed_point_levels)
  if (.density_method_uses_precomputed(
      density_method,
      allow_normal = TRUE
    ) && (nrow(requested_point_refs) == 0L || only_fixed_points)) {
    density_method <- "KDE"
  }
  coefficient_target <- .hypothesis_brma_formula_coefficient_target(
    object   = object,
    selected = selected
  )
  coefficient_level_targets <-
    .hypothesis_brma_formula_coefficient_level_targets(
      object     = object,
      selected   = selected,
      point_refs = point_refs
    )
  if (!is.null(coefficient_target)) {
    coefficient_target[["route"]] <-
      .hypothesis_brma_formula_transform_route(coefficient_target)
    .hypothesis_brma_check_formula_point_support(
      point_refs  = point_refs,
      target_info = coefficient_target
    )
    if (standardized_coefficients) {
      coefficient_target <- NULL
    }
  }
  if (standardized_coefficients) {
    coefficient_level_targets <- list()
  }

  if (identical(selected[["component"]], "random")) {
    out <- .hypothesis_brma_random(
      object                    = object,
      parameter                 = parameter,
      hypothesis                = hypothesis,
      standardized_coefficients = standardized_coefficients,
      conditional               = conditional,
      logBF                     = logBF,
      BF01                      = BF01,
      seed                      = seed,
      density_method            = density_method,
      density_control           = density_control,
      n_samples                 = n_samples,
      columns                   = columns,
      parameter_label           = parameter_label
    )
    return(.hypothesis_brma_restore_hypothesis_labels(
      out             = out,
      hypothesis      = display_hypothesis,
      parameter_label = parameter_label
    ))
  }

  sample_parameter <- .as_mixed_posteriors_parameters(object, parameter)
  if (!is.null(coefficient_target) &&
      !is.null(coefficient_target[["route"]][["weights"]])) {
    # The target's prior density combines the priors of every coordinate it
    # weights; their mixed posteriors carry the mixture components that
    # BayesTools needs for the target's components.
    sample_parameter <- unique(c(
      sample_parameter,
      .hypothesis_brma_target_prior_parameters(
        object  = object,
        weights = coefficient_target[["route"]][["weights"]]
      )
    ))
  }
  samples <- .brma_as_mixed_posteriors(
    object           = object,
    parameters       = sample_parameter,
    conditional      = if (conditional) parameter else NULL,
    conditional_rule = "AND",
    transform_scaled = !standardized_coefficients,
    n_prior_samples  = n_samples
  )
  if (!is.null(coefficient_target) &&
      identical(coefficient_target[["route"]][["type"]], "unsupported")) {
    stop(coefficient_target[["route"]][["reason"]], call. = FALSE)
  }
  if (!is.null(coefficient_target) &&
      identical(coefficient_target[["route"]][["type"]], "exp_affine")) {
    if (!identical(density_method, "KDE")) {
      stop(
        "The requested nonlinear fitted-scale hypothesis for '",
        parameter, "' is supported only with density_method = 'KDE'. ",
        "qCMDE/IWMDE ordinates support only direct parameter or level point ",
        "hypotheses with an exact linear fitted-scale map.",
        call. = FALSE
      )
    }
    out <- .hypothesis_brma_exp_affine_kde(
      object      = object,
      samples     = samples,
      hypothesis  = hypothesis,
      parameter   = parameter,
      target_info = coefficient_target,
      conditional = conditional,
      logBF       = logBF,
      BF01        = BF01,
      seed        = seed,
      n_samples   = n_samples,
      columns     = columns
    )
    return(.hypothesis_brma_restore_hypothesis_labels(
      out             = out,
      hypothesis      = display_hypothesis,
      parameter_label = parameter_label
    ))
  }
  if (!is.null(coefficient_target)) {
    coefficient_target <- .hypothesis_brma_formula_prior_target(
      object      = object,
      samples     = samples,
      hypothesis  = hypothesis,
      target_info = coefficient_target
    )
    prior_densities <- BayesTools::posterior_metadata(samples, "prior_densities")
    prior_densities[[parameter]] <- coefficient_target[["prior_density"]]
    BayesTools::posterior_metadata(samples, "prior_densities") <- prior_densities
  }
  for (level in names(coefficient_level_targets)) {
    coefficient_level_targets[[level]] <-
      .hypothesis_brma_formula_prior_target(
        object       = object,
        samples      = samples,
        hypothesis   = hypothesis,
        target_info  = coefficient_level_targets[[level]],
        point_values = point_refs[["value"]][
          !is.na(point_refs[["level"]]) &
            point_refs[["level"]] == level
        ],
        force_linear = TRUE
      )
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
    BayesTools::posterior_metadata(
      posterior,
      "prior_density"
    ) <- coefficient_target[["prior_density"]]
  }
  for (level in names(coefficient_level_targets)) {
    if (is.list(posterior) && level %in% names(posterior)) {
      BayesTools::posterior_metadata(posterior[[level]], "prior_density") <-
        coefficient_level_targets[[level]][["prior_density"]]
    }
  }

  if (level_contrast) {
    out <- .hypothesis_brma_level_contrast_BF(
      object                    = object,
      posterior                 = posterior,
      hypothesis                = hypothesis,
      parameter                 = parameter,
      standardized_coefficients = standardized_coefficients,
      density_method            = density_method,
      density_control           = density_control,
      logBF                     = logBF,
      BF01                      = BF01,
      seed                      = seed,
      columns                   = columns
    )
    return(.hypothesis_brma_restore_hypothesis_labels(
      out             = out,
      hypothesis      = display_hypothesis,
      parameter_label = parameter_label
    ))
  }

  if (.density_method_uses_precomputed(density_method, allow_normal = TRUE)) {
    .hypothesis_brma_check_fixed_level_points(
      object       = object,
      posterior    = posterior,
      hypothesis   = hypothesis,
      parameter    = parameter,
      fixed_levels = fixed_point_levels,
      seed         = seed
    )
    posterior <- .hypothesis_brma_attach_iwmde(
      object                   = object,
      posterior                = posterior,
      parameter                = parameter,
      parameter_label          = parameter_label,
      hypothesis               = hypothesis,
      conditional              = if (conditional) parameter else NULL,
      n_points                 = density_control[["n_points"]],
      samples                  = density_control[["samples"]],
      target_relative_mcse     = density_control[["target_relative_mcse"]],
      normalization_points     = density_control[["normalization_points"]],
      normalization_prob       = density_control[["normalization_prob"]],
      integration_control      = density_control[["integration_control"]],
      density_method           = density_method,
      n_samples                = n_samples,
      parameter_spec           = if (is.null(coefficient_target)) {
        NULL
      } else {
        coefficient_target[["parameter_spec"]]
      },
      level_parameter_specs    = lapply(
        coefficient_level_targets,
        `[[`,
        "parameter_spec"
      )
    )
  }

  out <- tryCatch(
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
    ),
    error = function(error) {
      .hypothesis_brma_stop_point_mass(object, error)
    }
  )

  if (.density_method_uses_precomputed(density_method, allow_normal = TRUE)) {
    out <- .hypothesis_brma_append_iwmde_warnings(
      table     = out,
      posterior = posterior
    )
  }

  out <- .hypothesis_brma_restore_hypothesis_labels(
    out             = out,
    hypothesis      = display_hypothesis,
    parameter_label = parameter_label
  )

  return(out)
}


.hypothesis_brma_multiple_parameters <- function(
    object, hypothesis, selections, standardized_coefficients, conditional,
    conditional_omitted, logBF, BF01, seed, density_method, density_control,
    n_samples, columns) {

  components <- vapply(selections, `[[`, character(1), "component")
  parameters <- vapply(selections, `[[`, character(1), "parameter")
  keys       <- paste(components, parameters, sep = "\r")
  group_keys <- unique(keys)
  groups     <- lapply(group_keys, function(key) which(keys == key))
  statements <- BayesTools::hypothesis_render(hypothesis)

  for (component in unique(components)) {
    .hypothesis_brma_check_supported_component(component)
  }

  results <- lapply(groups, function(rows) {

    arguments <- list(
      object                    = object,
      hypothesis                = unname(statements[rows]),
      component                 = components[[rows[[1L]]]],
      standardized_coefficients = standardized_coefficients,
      logBF                     = logBF,
      BF01                      = BF01,
      seed                      = seed,
      density_method            = density_method,
      density_control           = density_control,
      n_samples                 = n_samples,
      columns                   = columns
    )
    if (!conditional_omitted) {
      arguments[["conditional"]] <- conditional
    }

    do.call(hypothesis.brma, arguments)
  })

  .hypothesis_brma_bind_parameter_results(
    results    = results,
    groups     = groups,
    hypothesis = hypothesis
  )
}


.hypothesis_brma_result_statements <- function(result) {

  hypothesis <- attr(result, "hypothesis_ast", exact = TRUE)
  if (!inherits(hypothesis, "BayesTools_hypothesis_ast") ||
      !is.list(hypothesis[["statements"]]) ||
      length(hypothesis[["statements"]]) == 0L) {
    return(NULL)
  }

  hypothesis[["statements"]]
}


.hypothesis_brma_bind_parameter_results <- function(
    results, groups, hypothesis) {

  if (length(results) == 0L) {
    stop("Internal error: no hypothesis parameter-group results were produced.",
         call. = FALSE)
  }
  n_statements <- length(hypothesis[["statements"]])
  group_sizes  <- vapply(groups, length, integer(1))
  result_sizes <- vapply(results, nrow, integer(1))
  valid_results <- vapply(
    results,
    inherits,
    logical(1),
    what = "BayesTools_hypothesis_BF"
  )
  same_columns <- vapply(
    results,
    function(result) identical(names(result), names(results[[1L]])),
    logical(1)
  )
  result_statement_sizes <- vapply(results, function(result) {

    statements <- .hypothesis_brma_result_statements(result)
    if (is.null(statements)) -1L else length(statements)
  }, integer(1))
  if (!all(valid_results) ||
      !identical(group_sizes, result_sizes) ||
      !identical(group_sizes, result_statement_sizes) ||
      sum(group_sizes) != n_statements || !all(same_columns)) {
    stop("Internal error: hypothesis parameter-group results are misaligned.",
         call. = FALSE)
  }

  out <- do.call(rbind, results)
  concatenated_rows <- unlist(groups, use.names = FALSE)
  restore_order     <- order(concatenated_rows)
  # Row names as in BayesTools::hypothesis_BF(): rows of several statements on
  # one quantity carry their statement number in the full hypothesis, not the
  # number within their parameter group.
  row_names <- .hypothesis_brma_row_names(
    labels     = sub(
      " \\([0-9]+\\)$", "",
      unlist(lapply(results, rownames), use.names = FALSE)
    ),
    statements = concatenated_rows
  )

  raw_BF <- unlist(lapply(results, function(result) {
    attr(result, "raw_BF", exact = TRUE)
  }), use.names = FALSE)
  if (length(raw_BF) != n_statements) {
    stop("Internal error: hypothesis parameter-group metadata are misaligned.",
         call. = FALSE)
  }

  warnings <- list()
  row_offset <- 0L
  for (i in seq_along(results)) {
    result <- results[[i]]
    rows   <- row_offset + seq_len(nrow(result))
    result_warnings <- attr(result, "warnings", exact = TRUE)
    if (length(result_warnings) > 0L && !is.null(names(result_warnings))) {
      warning_rows <- match(names(result_warnings), rownames(result))
      matched      <- !is.na(warning_rows)
      names(result_warnings)[matched] <- row_names[rows[warning_rows[matched]]]
    }
    warnings[[i]] <- result_warnings
    row_offset    <- row_offset + nrow(result)
  }
  rownames(out) <- row_names

  bound_operator <- unlist(lapply(results, function(result) {
    operator <- attr(result[["BF"]], "bound_operator", exact = TRUE)
    if (is.null(operator)) {
      return(rep(NA_character_, nrow(result)))
    }
    rep_len(as.character(operator), nrow(result))
  }), use.names = FALSE)

  out <- out[restore_order, , drop = FALSE]
  attr(out, "raw_BF")         <- raw_BF[restore_order]
  attr(out, "hypothesis_ast") <- hypothesis
  warnings <- .hypothesis_brma_unique_named_warnings(
    unlist(warnings, use.names = TRUE)
  )
  if (length(warnings) > 0L && !is.null(names(warnings))) {
    warning_rows <- match(names(warnings), rownames(out))
    warnings <- c(
      warnings[is.na(warning_rows)],
      warnings[order(warning_rows, na.last = NA)]
    )
  }
  attr(out, "warnings") <- warnings
  attr(out[["BF"]], "bound_operator") <- bound_operator[restore_order]

  footnotes <- unique(unlist(lapply(results, function(result) {
    attr(result, "footnotes", exact = TRUE)
  }), use.names = FALSE))
  attr(out, "footnotes") <- if (length(footnotes) > 0L) footnotes else NULL

  diagnostics <- lapply(results, function(result) {
    attr(result, "density_diagnostics", exact = TRUE)
  })
  diagnostics <- diagnostics[vapply(
    diagnostics,
    function(diagnostic) is.data.frame(diagnostic) && nrow(diagnostic) > 0L,
    logical(1)
  )]
  if (length(diagnostics) > 0L) {
    attr(out, "density_diagnostics") <- do.call(rbind, diagnostics)
  } else {
    attr(out, "density_diagnostics") <- NULL
  }
  attr(out, "rownames") <- FALSE

  return(out)
}

.hypothesis_brma_random <- function(
    object, parameter, hypothesis, standardized_coefficients,
    conditional, logBF, BF01, seed, density_method, density_control = NULL,
    n_samples, columns, parameter_label = parameter) {

  precomputed <- .density_method_uses_precomputed(
    density_method,
    allow_normal = TRUE
  )
  if (conditional && precomputed) {
    stop(
      "Conditional random-effect hypotheses support ",
      "'density_method = \"KDE\"' only.",
      call. = FALSE
    )
  }

  selected <- .brma_random_parameter_select(
    object                    = object,
    parameter                 = parameter,
    standardized_coefficients = standardized_coefficients
  )
  if (identical(selected[["spec"]][["status"]], "structural")) {
    stop(
      "Hypothesis tests are not defined for fixed random-effect quantity '",
      selected[["entry"]][["term"]], "'.",
      call. = FALSE
    )
  }
  # The mixed posterior carries the catalog support, the canonical prior
  # density of BayesTools, and the declared inclusion- and allocation-gate
  # atoms; conditioning keeps the draws in the quantity's inclusion event.
  samples <- .brma_random_parameter_mixed_posterior(
    object                    = object,
    parameter                 = parameter,
    standardized_coefficients = standardized_coefficients,
    conditional               = conditional,
    selected                  = selected
  )
  prior_density <- BayesTools::posterior_metadata(
    samples[[parameter]],
    "prior_density"
  )
  defined_footnote <- .brma_random_parameter_defined_footnote(
    label             = selected[["spec"]][["label"]],
    samples           = samples[[parameter]],
    posterior_defined = .brma_random_parameter_defined_share(
      object,
      samples[[parameter]]
    )
  )

  point_refs <- BayesTools::hypothesis_parse_point_reference(
    hypothesis     = hypothesis,
    allow_compound = TRUE
  )
  target <- NULL
  if (nrow(point_refs) > 0L) {
    if (any(!point_refs[["direct"]])) {
      stop(
        "Point-null tests for random-effect quantities require a direct scalar ",
        "parameter reference.",
        call. = FALSE
      )
    }
    point_status <- .brma_random_parameter_point_status(
      selected      = selected,
      prior_density = prior_density,
      values        = unique(point_refs[["value"]])
    )
    zero_alternative <- if (any(point_refs[["value"]] == 0)) {
      .brma_random_parameter_zero_boundary_alternative(object, selected)
    } else {
      NULL
    }
    if (!is.null(zero_alternative)) {
      stop(
        "Point-null Bayes factors are unavailable for allocation-derived ",
        "random-effect quantity '", selected[["spec"]][["label"]],
        "' at 0 because zero is a nonregular product boundary of the common ",
        "scale and allocation weight. Test '", zero_alternative,
        "' to compare omission of this component.",
        call. = FALSE
      )
    }
    if (precomputed) {
      target <- .brma_random_parameter_density_target(
        object,
        parameter,
        operation = "point hypotheses"
      )
      if (is.null(target[["parameter"]])) {
        stop(target[["reason"]], call. = FALSE)
      }
    }
    if (precomputed && !is.null(target[["display_transform"]])) {
      source_values <- BayesTools::parameter_transform_inverse(
        point_refs[["value"]],
        target[["display_transform"]]
      )
      jacobian <- BayesTools::parameter_transform_jacobian(
        source_values,
        target[["display_transform"]]
      )
      singular <- !is.finite(source_values) | !is.finite(jacobian) |
        jacobian <= 0
      if (any(singular)) {
        value <- point_refs[["value"]][which(singular)[[1L]]]
        stop(
          "Point-null Bayes factors are unavailable for random-effect ",
          "quantity '", selected[["entry"]][["term"]], "' at ", value,
          " because its public transformation is singular at that support ",
          "boundary. Use the corresponding directly modeled scale or a ",
          "region hypothesis.",
          call. = FALSE
        )
      }
    }
    .brma_random_parameter_check_point_status(selected, point_status)
    if (!precomputed) {
      support <- .brma_random_parameter_support(selected)
      values <- point_refs[["value"]]
      at_boundary <- (is.finite(support[1L]) & values <= support[1L]) |
        (is.finite(support[2L]) & values >= support[2L])
      if (any(at_boundary)) {
        stop(
          "Point-null Bayes factors at the support boundary are not available ",
          "for random-effect quantity '", selected[["entry"]][["term"]], "'.",
          call. = FALSE
        )
      }
    }
  }

  if (is.null(prior_density)) {
    # Without a canonical prior density, region hypotheses take the prior
    # probabilities from prior draws (point hypotheses stopped above).
    out <- .hypothesis_brma_random_prior_draws(
      object                    = object,
      parameter                 = parameter,
      selected                  = selected,
      posterior                 = samples[[parameter]],
      hypothesis                = hypothesis,
      standardized_coefficients = standardized_coefficients,
      conditional               = conditional,
      logBF                     = logBF,
      BF01                      = BF01,
      seed                      = seed,
      n_samples                 = n_samples,
      columns                   = columns,
      density_method            = if (precomputed) "KDE" else density_method
    )
    defined_footnote <- attr(out, "defined_footnote", exact = TRUE)
    attr(out, "defined_footnote") <- NULL
    if (!is.null(defined_footnote)) {
      attr(out, "footnotes") <- c(attr(out, "footnotes"), defined_footnote)
    }
    return(out)
  }

  marginal <- BayesTools::marginal_posterior(
    samples       = samples,
    parameter     = parameter,
    prior_samples = TRUE,
    use_formula   = FALSE,
    n_samples     = n_samples
  )
  if (!precomputed || nrow(point_refs) == 0L) {
    out <- BayesTools::hypothesis_BF(
      posterior      = marginal,
      hypothesis     = hypothesis,
      parameter      = parameter,
      logBF          = logBF,
      BF01           = BF01,
      seed           = seed,
      columns        = columns,
      density_method = if (precomputed) "KDE" else density_method
    )
    if (!is.null(defined_footnote)) {
      attr(out, "footnotes") <- c(attr(out, "footnotes"), defined_footnote)
    }
    return(out)
  }

  if (is.null(density_control[["normalization_points"]])) {
    density_control[["normalization_points"]] <- max(
      50L,
      density_control[["n_points"]]
    )
  }
  context        <- .iwmde_context(object, density_control[["integration_control"]])
  estimate_cache <- .iwmde_estimate_cache()
  marginal <- .hypothesis_brma_attach_iwmde_scalar(
    posterior                = marginal,
    raw_posterior            = samples[[parameter]],
    context                  = context,
    estimate_cache           = estimate_cache,
    parameter                = target[["parameter"]],
    parameter_label          = parameter_label,
    value                    = unique(point_refs[["value"]]),
    conditional              = NULL,
    n_points                 = density_control[["n_points"]],
    samples                  = density_control[["samples"]],
    target_relative_mcse     = density_control[["target_relative_mcse"]],
    normalization_points     = density_control[["normalization_points"]],
    normalization_prob       = density_control[["normalization_prob"]],
    integration_control      = density_control[["integration_control"]],
    density_method           = density_method,
    parameter_spec           = target[["parameter_spec"]],
    display_transform        = target[["display_transform"]]
  )

  out <- BayesTools::hypothesis_BF(
    posterior      = marginal,
    hypothesis     = hypothesis,
    parameter      = parameter,
    logBF          = logBF,
    BF01           = BF01,
    seed           = seed,
    columns        = columns,
    density_method = "precomputed"
  )
  .hypothesis_brma_append_iwmde_warnings(
    table     = out,
    posterior = marginal,
    parameter = parameter
  )
}


# Region hypotheses on a random-effect quantity without a canonical prior
# density: prior probabilities from the quantity's prior draws.
.hypothesis_brma_random_prior_draws <- function(
    object, parameter, selected, posterior, hypothesis,
    standardized_coefficients, conditional, logBF, BF01, seed, n_samples,
    columns, density_method) {

  if (conditional) {
    stop(
      "Conditional hypotheses are unavailable for random-effect quantity '",
      selected[["spec"]][["label"]], "' because its prior density is ",
      "unavailable.",
      call. = FALSE
    )
  }
  prior <- .brma_random_parameter_select(
    object                    = object,
    parameter                 = parameter,
    standardized_coefficients = standardized_coefficients,
    prior                     = TRUE,
    n_prior_samples           = n_samples,
    seed                      = seed
  )
  prior_values  <- as.numeric(prior[["samples"]][, 1L])
  prior_defined <- .brma_random_parameter_defined_draws(
    prior_values,
    prior[["samples"]],
    selected[["spec"]][["label"]]
  )

  out <- BayesTools::hypothesis_BF(
    posterior      = as.numeric(posterior),
    prior          = prior_values[prior_defined],
    hypothesis     = hypothesis,
    parameter      = parameter,
    logBF          = logBF,
    BF01           = BF01,
    seed           = seed,
    columns        = columns,
    density_method = density_method
  )
  attr(out, "defined_footnote") <- .brma_random_parameter_defined_footnote(
    label             = selected[["spec"]][["label"]],
    samples           = prior[["samples"]],
    posterior_defined = .brma_random_parameter_defined_share(object, posterior),
    prior_defined     = prior_defined
  )
  out
}

.hypothesis_brma_formula_coefficient_target <- function(
    object, selected) {

  entry <- selected[["entry"]]
  if (selected[["component"]] %in% c("random", "bias") ||
      is.null(entry) ||
      identical(entry[["role"]], "formula_coefficient_group")) {
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
      !identical(transform[["schema_version"]], 2L)) {
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


.hypothesis_brma_formula_coefficient_level_targets <- function(
    object, selected, point_refs) {

  entry  <- selected[["entry"]]
  levels <- unique(point_refs[["level"]][!is.na(point_refs[["level"]])])
  if (selected[["component"]] %in% c("random", "bias") ||
      is.null(entry) ||
      !identical(entry[["role"]], "formula_coefficient_group") ||
      length(levels) == 0L) {
    return(list())
  }

  formula_parameter <- entry[["formula_parameter"]]
  if (is.null(formula_parameter) || length(formula_parameter) != 1L ||
      is.na(formula_parameter) || !nzchar(formula_parameter)) {
    return(list())
  }
  transform <- BayesTools::JAGS_formula_coefficient_transform(
    fit          = object[["fit"]],
    parameter    = formula_parameter,
    target_scale = "original"
  )
  if (!inherits(transform, "BayesTools_formula_coefficient_transform") ||
      !identical(transform[["schema_version"]], 2L)) {
    stop(
      "Formula coefficient transformation metadata are unsupported. Refit ",
      "the model with the current BayesTools version.",
      call. = FALSE
    )
  }

  occurrences <- selected[["resolution"]][["occurrences"]]
  quantities  <- BayesTools::parameter_catalog(object[["fit"]])[["quantities"]]
  term_label  <- .hypothesis_brma_alias_label(
    selected[["aliases"]],
    selected[["parameter"]]
  )
  out <- lapply(levels, function(level) {
    # Messages name the level by its selector, never by its backend
    # coordinate: with levels 1:4, coordinate 'mu_g[2]' is level 3.
    selector <- paste0(term_label, "[", level, "]")
    level_occurrences <- occurrences[
      !is.na(occurrences[["level"]]) & occurrences[["level"]] == level,
      ,
      drop = FALSE
    ]
    level_name <- unique(level_occurrences[["canonical_name"]])
    if (length(level_name) != 1L) {
      stop(
        "Factor level '", selector, "' is ambiguous in the fitted parameter ",
        "catalog.",
        call. = FALSE
      )
    }
    quantity <- quantities[
      quantities[["quantity_id"]] %in%
        unique(level_occurrences[["quantity_id"]]),
      ,
      drop = FALSE
    ]
    # Canonical level names are level labels, never coordinates: a level is a
    # fitted coefficient only when it is structurally one coordinate (a direct
    # level cell), never by a unit design row of mean-difference or
    # orthonormal coding.
    target <- if (nrow(quantity) == 1L) {
      .brma_catalog_level_coordinate(quantity, quantities)
    }
    target_i <- if (is.null(target)) {
      NA_integer_
    } else {
      match(target, transform[["target_names"]])
    }
    if (is.na(target_i)) {
      if (.hypothesis_brma_level_quantity_fixed(quantity)) {
        return(NULL)
      }
      key <- if (nrow(quantity) == 1L) quantity[["extraction_key"]][[1L]]
      if (is.null(target) && is.list(key) &&
          identical(key[["type"]], "factor_level") &&
          length(key[["dependencies"]]) > 0L) {
        .hypothesis_brma_stop_combined_level(
          selected   = selected,
          quantities = quantities,
          level      = level
        )
      }
      stop(
        "Resolved formula coefficient '", level_name,
        "' is absent from the fitted coefficient transformation.",
        call. = FALSE
      )
    }
    # The level's weights on the original-scale coefficients: a direct level
    # cell is its one coordinate.
    target_info <- list(
      formula_parameter = formula_parameter,
      target            = target,
      target_i          = target_i,
      transform         = transform,
      level_selector    = selector,
      level_weights     = stats::setNames(1, target)
    )
    target_info[["route"]] <-
      .hypothesis_brma_formula_transform_route(target_info)
    target_info
  })
  names(out) <- levels
  out <- out[!vapply(out, is.null, logical(1))]

  return(out)
}


# A level quantity fixed by the contrast (the treatment reference level).
.hypothesis_brma_level_quantity_fixed <- function(quantity) {

  nrow(quantity) == 1L &&
    (quantity[["status"]] %in% c("fixed", "structural") ||
       is.finite(quantity[["fixed_value"]]))
}


# Levels of the selected factor term that point hypotheses reference and that
# the contrast fixes (the treatment reference level).
.hypothesis_brma_fixed_point_levels <- function(object, selected, point_refs) {

  entry  <- selected[["entry"]]
  levels <- unique(point_refs[["level"]][!is.na(point_refs[["level"]])])
  if (length(levels) == 0L || is.null(entry) ||
      selected[["component"]] %in% c("random", "bias") ||
      !identical(entry[["role"]], "formula_coefficient_group")) {
    return(character())
  }

  occurrences <- selected[["resolution"]][["occurrences"]]
  quantities  <- BayesTools::parameter_catalog(object[["fit"]])[["quantities"]]
  fixed <- vapply(levels, function(level) {
    ids <- occurrences[["quantity_id"]][
      !is.na(occurrences[["level"]]) & occurrences[["level"]] == level
    ]
    .hypothesis_brma_level_quantity_fixed(quantities[
      quantities[["quantity_id"]] %in% unique(ids),
      ,
      drop = FALSE
    ])
  }, logical(1))

  return(levels[fixed])
}


# Point statements on a level fixed by the contrast involve no posterior
# density: the level's declared atom decides them. Evaluate them through the
# KDE route of hypothesis.brma(), which handles declared atoms, so that a
# qCMDE/IWMDE call reports the same reason as KDE instead of failing for lack
# of an ordinate.
.hypothesis_brma_check_fixed_level_points <- function(
    object, posterior, hypothesis, parameter, fixed_levels, seed) {

  refs <- BayesTools::hypothesis_parse_point_reference(
    hypothesis     = hypothesis,
    allow_compound = TRUE
  )
  statements <- unique(refs[["hypothesis"]][
    refs[["direct"]] & refs[["parameter"]] == parameter &
      !is.na(refs[["level"]]) & refs[["level"]] %in% fixed_levels
  ])
  if (length(statements) == 0L) {
    return(invisible(TRUE))
  }
  tryCatch(
    BayesTools::hypothesis_BF(
      posterior      = posterior,
      hypothesis     = statements,
      parameter      = parameter,
      seed           = seed,
      density_method = "KDE"
    ),
    error = function(error) {
      .hypothesis_brma_stop_point_mass(object, error)
    }
  )

  return(invisible(TRUE))
}


# A spike-and-slab component puts a declared point mass at the null, so the
# Savage-Dickey ratio is structurally unavailable rather than numerically
# rejected. Name the supported alternative; rethrow other errors unchanged.
# BayesTools signals prior and posterior point masses at the null with classed
# conditions; they are matched by class, never by message.
.hypothesis_brma_stop_point_mass <- function(object, error) {

  point_mass <- inherits(error, c(
    "BayesTools_point_mass_at_null",
    "BayesTools_posterior_point_mass_at_null"
  ))
  if (!inherits(object, "RoBMA") || !point_mass) {
    stop(error)
  }
  stop(
    conditionMessage(error),
    " This parameter has a null component, so its evidence against the ",
    "null is the inclusion Bayes factor reported by 'summary()' and ",
    "'summary_models()'.",
    call. = FALSE
  )
}


# A factor level that is not structurally one fitted coordinate (levels of
# mean-difference and orthonormal coding, ordered levels beyond the first
# increment) has no fitted coefficient whose prior ordinate a point hypothesis
# could use.
.hypothesis_brma_stop_combined_level <- function(selected, quantities, level) {

  label    <- .hypothesis_brma_alias_label(
    selected[["aliases"]],
    selected[["parameter"]]
  )
  selector <- paste0(label, "[", level, "]")
  members  <- unlist(
    selected[["entry"]][["member_quantity_ids"]],
    use.names = FALSE
  )
  others   <- setdiff(
    quantities[["component"]][quantities[["quantity_id"]] %in% members],
    level
  )

  stop(
    "Point hypotheses on factor level '", selector, "' are not supported: ",
    "the level is a linear combination of the fitted contrast coefficients ",
    "(mean-difference, orthonormal, or ordered contrasts), not a fitted ",
    "coefficient itself, and point hypotheses on a single level require a ",
    "level fitted as its own coefficient (treatment or independent ",
    "contrasts). Test the level with a region hypothesis such as '",
    selector, " > 0'",
    if (length(others) > 0L) {
      paste0(" or a level contrast such as '", selector, " = ", label, "[",
             others[[1L]], "]'")
    },
    ".",
    call. = FALSE
  )
}


# How messages name a formula coefficient target: a factor level by its
# selector (e.g. 'g[10]'), never by the backend coordinate that the target
# holds; a scalar coefficient by its parameter name.
.hypothesis_brma_formula_target_name <- function(target_info) {

  if (!is.null(target_info[["level_selector"]])) {
    return(target_info[["level_selector"]])
  }

  return(target_info[["target"]])
}


.hypothesis_brma_formula_target_description <- function(target_info) {

  paste0(
    if (is.null(target_info[["level_selector"]])) {
      "transformed coefficient"
    } else {
      "factor level"
    },
    " '", .hypothesis_brma_formula_target_name(target_info), "'"
  )
}


.hypothesis_brma_formula_prior_target <- function(
    object, samples, hypothesis, target_info, point_values = NULL,
    force_linear = FALSE) {

  if (is.null(target_info[["route"]])) {
    target_info[["route"]] <-
      .hypothesis_brma_formula_transform_route(target_info)
  }
  # Factor levels are weighted combinations of the original-scale
  # coefficients; scalar coefficients are their own target.
  level_weights <- target_info[["level_weights"]]
  prior_density <- BayesTools::JAGS_formula_prior_density(
    fit          = object[["fit"]],
    parameter    = target_info[["formula_parameter"]],
    target       = if (is.null(level_weights)) target_info[["target"]],
    weights      = level_weights,
    target_scale = "original",
    context      = BayesTools::posterior_metadata(samples, "prior_context")
  )
  if (is.null(point_values)) {
    refs <- .hypothesis_brma_point_refs(
      hypothesis,
      target_info[["target"]],
      require_direct = FALSE
    )
    point_values <- refs[["value"]]
  }
  # The BayesTools exactness rule classifies each point value; a prior point
  # mass or a zero, infinite or undefined ordinate stops in hypothesis_BF()
  # with its classed condition.
  if (length(point_values) > 0L) {
    status <- BayesTools::prior_ordinate_status(prior_density, point_values)
    if (any(status[["condition"]] %in% "BayesTools_inexact_ordinate")) {
      stop(
        "The induced prior ordinate for ",
        .hypothesis_brma_formula_target_description(target_info),
        " is not exact enough for a point-null Bayes factor.",
        call. = FALSE
      )
    }
  }

  route   <- target_info[["route"]]
  weights <- route[["weights"]]
  parameter_spec <- if (identical(route[["type"]], "identity") &&
                        !force_linear) {
    list(type = "primitive", prior_density = prior_density)
  } else if (route[["type"]] %in% c("identity", "affine")) {
    list(type = "linear", weights = weights, prior_density = prior_density)
  } else {
    list(
      type   = "unsupported_formula_transform",
      reason = paste0(
        "qCMDE/IWMDE does not support the fitted nonlinear joint transform ",
        "for '", .hypothesis_brma_formula_target_name(target_info), "'. Use ",
        "density_method = 'KDE' or standardized_coefficients = TRUE."
      )
    )
  }

  target_info[["prior_density"]] <- prior_density
  target_info[["parameter_spec"]] <- parameter_spec
  return(target_info)
}


# The fitted prior-list entries (formula terms) that own the coordinates a
# formula coefficient target weights, from the coordinate table and the
# formula name maps.
.hypothesis_brma_target_prior_parameters <- function(object, weights) {

  coordinates <- BayesTools::parameter_coordinates(object[["fit"]])
  rows <- match(names(weights), coordinates[["coordinate_name"]])
  if (anyNA(rows)) {
    stop(
      "Fitted coordinate metadata of the hypothesis target are unavailable. ",
      "Refit the model with the current RoBMA/BayesTools build.",
      call. = FALSE
    )
  }
  coordinates <- coordinates[rows, , drop = FALSE]

  unique(unlist(lapply(
    unique(coordinates[["formula_parameter"]]),
    function(formula_parameter) {
      name_map <- .fitted_formula_name_map(object, formula_parameter)
      fixed    <- name_map[name_map[["kind"]] == "fixed", , drop = FALSE]
      terms    <- coordinates[["term"]][
        coordinates[["formula_parameter"]] == formula_parameter
      ]
      fixed[["jags_name"]][fixed[["term"]] %in% terms]
    }
  ), use.names = FALSE))
}


# The route of a formula coefficient target: its weights on the fitted
# coordinates, and the map type and map support that BayesTools declares for
# the target (identity, affine, exp_affine, or unsupported).
.hypothesis_brma_formula_transform_route <- function(target_info) {

  transform <- target_info[["transform"]]
  target    <- target_info[["target"]]
  name      <- .hypothesis_brma_formula_target_name(target_info)
  if (!inherits(transform, "BayesTools_formula_coefficient_transform") ||
      !identical(transform[["target_scale"]], "original")) {
    return(list(
      type   = "unsupported",
      reason = paste0(
        "The fitted coefficient transform for '", name,
        "' lacks the certified structural metadata required for hypothesis testing."
      )
    ))
  }

  target_i <- target_info[["target_i"]]
  weights  <- stats::setNames(
    as.numeric(transform[["matrix"]][target_i, , drop = FALSE]),
    colnames(transform[["matrix"]])
  )
  weights  <- weights[weights != 0]
  if (length(weights) == 0L) {
    return(list(
      type   = "unsupported",
      reason = paste0(
        "The fitted coefficient '", name,
        "' is structurally fixed and has no posterior hypothesis route."
      )
    ))
  }
  targets <- transform[["targets"]]
  row     <- if (is.data.frame(targets) &&
                 all(c("target", "map_type", "support") %in% names(targets))) {
    match(target, targets[["target"]])
  } else {
    NA_integer_
  }
  if (is.na(row)) {
    return(list(
      type   = "unsupported",
      reason = paste0(
        "The fitted coefficient transform for '", name,
        "' lacks the certified structural metadata required for hypothesis testing."
      )
    ))
  }
  map_type <- targets[["map_type"]][[row]]
  if (map_type %in% c("identity", "affine", "exp_affine")) {
    return(list(
      type    = map_type,
      weights = weights,
      support = targets[["support"]][[row]]
    ))
  }

  list(
    type   = "unsupported",
    reason = paste0(
      "The fitted nonlinear joint coefficient transform for '", name,
      "' is not supported by hypothesis()."
    )
  )
}


.hypothesis_brma_check_formula_point_support <- function(point_refs,
                                                         target_info) {

  support <- target_info[["route"]][["support"]]
  if (nrow(point_refs) == 0L || is.null(support)) {
    return(invisible(TRUE))
  }
  values <- point_refs[["value"]]
  outside_or_boundary <-
    !is.finite(values) |
    (is.finite(support[[1L]]) & values <= support[[1L]]) |
    (is.finite(support[[2L]]) & values >= support[[2L]])
  if (any(outside_or_boundary)) {
    stop(
      "Point-null value ", values[which(outside_or_boundary)[[1L]]],
      " is outside or on the boundary of the open support for ",
      .hypothesis_brma_formula_target_description(target_info), ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.hypothesis_brma_exp_affine_certify <- function(samples, target,
                                                conditional) {

  if (isTRUE(conditional)) {
    stop(
      "Nonlinear fitted-scale KDE hypotheses are unavailable for conditional ",
      "product-space posteriors.",
      call. = FALSE
    )
  }
  sample <- samples[[target]]
  if (is.null(sample) ||
      !inherits(sample, "mixed_posteriors.simple")) {
    stop(
      "Nonlinear fitted-scale KDE hypotheses require a certified scalar ",
      "mixed posterior.",
      call. = FALSE
    )
  }
  # BayesTools declares averaged (unconditioned) draws and atom-free draws in
  # their metadata.
  condition <- BayesTools::posterior_metadata(sample, "condition")
  if (!isTRUE(condition[["averaged"]])) {
    stop(
      "Nonlinear fitted-scale KDE hypotheses require structural evidence for ",
      "an unconditional posterior.",
      call. = FALSE
    )
  }

  if (!BayesTools::posterior_atoms_free(sample)) {
    stop(
      "Nonlinear fitted-scale KDE hypotheses require structural evidence that ",
      "the posterior is atom-free.",
      call. = FALSE
    )
  }

  prior_densities <- BayesTools::posterior_metadata(samples, "prior_densities")
  prior_density   <- prior_densities[[target]]
  prior_points    <- prior_density[["points"]]
  prior_atom_free <- inherits(prior_density, "prior_density") &&
    is.data.frame(prior_points) &&
    all(c("x", "p") %in% names(prior_points)) &&
    nrow(prior_points) == 0L
  if (!prior_atom_free) {
    stop(
      "Nonlinear fitted-scale KDE hypotheses require structural evidence that ",
      "the prior is atom-free.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.hypothesis_brma_exp_affine_kde <- function(
    object, samples, hypothesis, parameter, target_info, conditional, logBF,
    BF01, seed, n_samples, columns) {

  if (is.null(target_info) ||
      !identical(target_info[["route"]][["type"]], "exp_affine")) {
    stop("Internal error: the exp(affine) hypothesis route was not certified.",
         call. = FALSE)
  }
  target <- target_info[["target"]]
  .hypothesis_brma_exp_affine_certify(
    samples     = samples,
    target      = target,
    conditional = conditional
  )
  route_kind <- .hypothesis_brma_exp_affine_route_kind(hypothesis)
  if (!target %in% names(samples)) {
    stop("Transformed posterior draws for '", target, "' are unavailable.",
         call. = FALSE)
  }
  prior <- BayesTools::transform_prior_samples(
    fit       = object[["fit"]],
    n_samples = n_samples,
    seed      = seed
  )
  if (is.null(colnames(prior)) || !target %in% colnames(prior)) {
    stop("Transformed prior draws for '", target, "' are unavailable.",
         call. = FALSE)
  }

  posterior_values <- as.numeric(samples[[target]])
  prior_values     <- as.numeric(prior[, target])
  if (length(posterior_values) < 2L || length(prior_values) < 2L ||
      any(!is.finite(posterior_values)) || any(!is.finite(prior_values))) {
    stop("Finite transformed prior and posterior draws are required for '",
         target, "'.", call. = FALSE)
  }
  original_hypothesis <- hypothesis
  if (identical(route_kind, "point")) {
    if (any(posterior_values <= 0) || any(prior_values <= 0)) {
      stop(
        "Positive transformed prior and posterior draws are required for ",
        "exp(affine) point hypotheses.",
        call. = FALSE
      )
    }
    posterior_values <- log(posterior_values)
    prior_values     <- log(prior_values)
    hypothesis <- .hypothesis_brma_exp_affine_log_hypothesis(
      hypothesis = hypothesis
    )
  }
  posterior <- stats::setNames(data.frame(posterior_values), parameter)
  prior     <- stats::setNames(data.frame(prior_values), parameter)

  out <- BayesTools::hypothesis_BF(
    posterior      = posterior,
    prior          = prior,
    hypothesis     = hypothesis,
    parameter      = parameter,
    logBF          = logBF,
    BF01           = BF01,
    seed           = seed,
    columns        = columns,
    density_method = "KDE"
  )
  if (identical(route_kind, "point")) {
    out <- .hypothesis_brma_restore_hypothesis_labels(
      out        = out,
      hypothesis = original_hypothesis
    )
    out <- .hypothesis_brma_exp_affine_restore_density_scale(
      out        = out,
      hypothesis = original_hypothesis
    )
  }

  return(out)
}


.hypothesis_brma_exp_affine_restore_density_scale <- function(
    out, hypothesis) {

  density_columns <- intersect(c("prior", "posterior"), names(out))
  if (length(density_columns) == 0L) {
    return(out)
  }
  points <- vapply(hypothesis[["statements"]], function(statement) {

    values <- vapply(c("left", "right"), function(side_name) {
      statement[[side_name]][["value"]]
    }, numeric(1))
    if (!isTRUE(all.equal(values[[1L]], values[[2L]])) ||
        !is.finite(values[[1L]]) || values[[1L]] <= 0) {
      stop(
        "Internal error: exp(affine) point-density scales are invalid.",
        call. = FALSE
      )
    }

    values[[1L]]
  }, numeric(1))
  if (nrow(out) != length(points)) {
    stop("Internal error: transformed density rows are misaligned.",
         call. = FALSE)
  }
  for (column in density_columns) {
    out[[column]] <- out[[column]] / points
  }

  return(out)
}


.hypothesis_brma_exp_affine_route_kind <- function(hypothesis) {

  statements <- hypothesis[["statements"]]
  side_types <- unlist(lapply(statements, function(statement) {
    c(statement[["left"]][["type"]], statement[["right"]][["type"]])
  }), use.names = FALSE)
  point_side <- side_types %in% c("point", "not_point")
  if (any(point_side) && any(!point_side)) {
    stop(
      "exp(affine) hypotheses cannot mix point and region statements. ",
      "Evaluate point and directional hypotheses in separate calls.",
      call. = FALSE
    )
  }

  if (any(point_side)) "point" else "region"
}


.hypothesis_brma_exp_affine_log_hypothesis <- function(hypothesis) {

  statements <- hypothesis[["statements"]]
  transformed <- vapply(statements, function(statement) {

    sides <- lapply(c("left", "right"), function(side_name) {
      side <- statement[[side_name]]
      operator <- switch(
        side[["type"]],
        point     = "=",
        not_point = "!=",
        stop("Internal error: expected a direct point statement.",
             call. = FALSE)
      )
      expression <- side[["expression"]][["source"]]
      if (!is.character(expression) || length(expression) != 1L ||
          !nzchar(expression)) {
        stop("Internal error: direct point expression source is unavailable.",
             call. = FALSE)
      }
      paste(expression, operator, sprintf("%.17g", log(side[["value"]])))
    })
    if (isTRUE(statement[["explicit"]])) {
      paste(sides[[1L]], "vs", sides[[2L]])
    } else {
      sides[[1L]]
    }
  }, character(1))

  BayesTools::hypothesis_parse(transformed)
}


.hypothesis_brma_restore_hypothesis_labels <- function(
    out, hypothesis, parameter_label = NULL) {

  if (!is.data.frame(out)) {
    return(out)
  }

  statements <- hypothesis[["statements"]]
  result_label <- function(statement, side_name) {

    implicit_equality <- !isTRUE(statement[["explicit"]]) &&
      identical(statement[["left"]][["type"]], "point")
    if (implicit_equality) {
      side_name <- switch(side_name, left = "right", right = "left")
    }

    statement[[side_name]][["label"]]
  }
  transformed_statements <- .hypothesis_brma_result_statements(out)
  if (!is.list(transformed_statements) ||
      length(transformed_statements) != length(statements)) {
    stop("Internal error: transformed hypothesis metadata are misaligned.",
         call. = FALSE)
  }

  transformed_left  <- vapply(
    transformed_statements,
    result_label,
    side_name = "left",
    character(1)
  )
  transformed_right <- vapply(
    transformed_statements,
    result_label,
    side_name = "right",
    character(1)
  )
  display_left       <- vapply(
    statements,
    result_label,
    side_name = "left",
    character(1)
  )
  display_right      <- vapply(
    statements,
    result_label,
    side_name = "right",
    character(1)
  )
  if (nrow(out) == length(statements)) {
    statement_i <- seq_along(statements)
  } else if (length(statements) == 1L) {
    statement_i <- rep(1L, nrow(out))
  } else if (all(c("Alternative", "Null") %in% names(out))) {
    transformed_key <- paste(transformed_left, transformed_right, sep = "\r")
    output_key      <- paste(out[["Alternative"]], out[["Null"]], sep = "\r")
    statement_i     <- match(output_key, transformed_key)
  } else if ("Alternative" %in% names(out)) {
    statement_i <- match(out[["Alternative"]], transformed_left)
  } else if ("Null" %in% names(out)) {
    statement_i <- match(out[["Null"]], transformed_right)
  } else {
    statement_i <- rep(NA_integer_, nrow(out))
  }
  if (anyNA(statement_i)) {
    stop("Internal error: transformed hypothesis result rows are misaligned.",
         call. = FALSE)
  }
  if ("Alternative" %in% names(out)) {
    out[["Alternative"]] <- display_left[statement_i]
  }
  if ("Null" %in% names(out)) {
    out[["Null"]] <- display_right[statement_i]
  }

  attr(out, "hypothesis_ast") <- hypothesis

  if (!is.null(parameter_label)) {
    old_rownames  <- rownames(out)
    bracket_start <- regexpr("[", old_rownames, fixed = TRUE)
    suffix        <- rep("", length(old_rownames))
    has_suffix    <- bracket_start > 0L
    suffix[has_suffix] <- substring(
      old_rownames[has_suffix],
      bracket_start[has_suffix]
    )
    # Statement numbers of repeated rows are renumbered for the display labels.
    suffix <- sub(" \\([0-9]+\\)$", "", suffix)
    out <- .hypothesis_brma_set_row_names(
      out       = out,
      row_names = .hypothesis_brma_row_names(
        labels     = paste0(parameter_label, suffix),
        statements = statement_i
      )
    )
  }
  attr(out, "rownames") <- FALSE

  return(out)
}


# Row names of hypothesis tables, as in BayesTools::hypothesis_BF(): rows of
# several statements on one quantity carry the statement number ("mu (1)",
# "mu (2)"). "mu1" reads like another parameter, and "mu.1" collides with the
# names `[.data.frame` gives duplicated rows.
.hypothesis_brma_row_names <- function(labels, statements) {

  repeated <- labels %in% labels[duplicated(labels)]
  labels[repeated] <- paste0(labels[repeated], " (", statements[repeated], ")")

  return(make.unique(labels, sep = " "))
}


# Rename table rows together with the table warnings keyed by them.
.hypothesis_brma_set_row_names <- function(out, row_names) {

  old_rownames  <- rownames(out)
  rownames(out) <- row_names

  warnings <- attr(out, "warnings", exact = TRUE)
  if (!is.null(warnings) && !is.null(names(warnings))) {
    warning_rows <- match(names(warnings), old_rownames)
    matched      <- !is.na(warning_rows)
    names(warnings)[matched] <- row_names[warning_rows[matched]]
    attr(out, "warnings") <- warnings
  }

  return(out)
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
