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
#' qCMDE and IWMDE evaluate one scale formula: they are unavailable for
#' \code{brma.mv()} models with component-specific scale formulas (a named
#' \code{scale} list, e.g. one scale formula per random component), whose
#' point-null hypotheses use \code{density_method = "KDE"}.
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
#' terms used in multiple model components. A statement whose references
#' still name several model parameters stops with an error of class
#' \code{RoBMA_hypothesis_ambiguous}, with the parent class
#' \code{RoBMA_hypothesis_statement}: the statement or \code{component} needs
#' to change, the test itself is not unavailable. An ambiguous reference that
#' BayesTools refuses with the class \code{BayesTools_parameter_ambiguous}
#' (for example the contrast coefficient \code{"g{1}"} of a factor term in
#' the location and the scale formula) has these classes too, before its
#' BayesTools classes. A statement whose parameters belong to another
#' component than \code{component}, including a name that is known only
#' outside \code{component} next to references of \code{component} (with or
#' without factor levels, e.g. \code{"log_tau_g[a] > mu_g[a]"} with
#' \code{component = "mods"}), stops with the classes
#' \code{RoBMA_hypothesis_statement} and \code{RoBMA_component_mismatch}
#' (the class of a component mismatch in [plot.brma()]). The exception is
#' the alias of a factor term of both the location and the scale formula
#' (for example \code{g} for \code{mods = ~ g, scale = ~ g}) next to a name
#' that is unknown within \code{component}, such as a name of the other
#' component or a display label (for example \code{"g > log_tau_intercept"}
#' with \code{component = "mods"}): the statement is resolved without
#' \code{component}, where the alias is ambiguous, and stops with
#' \code{RoBMA_hypothesis_statement} followed by the classes of the
#' BayesTools refusal of its first unresolved reference, the ambiguity of the
#' alias (\code{BayesTools_parameter_ambiguous}, without
#' \code{RoBMA_hypothesis_ambiguous}, as \code{component} is set already) or
#' the unknown name (\code{BayesTools_parameter_not_found}, for example in
#' \code{"log_tau_intercept < g"} or, with a level of the alias, in
#' \code{"g[a] > log_tau_intercept"}). The display labels
#' of the summary tables (for example \code{"(mu) intercept"},
#' \code{"exp(intercept)"}, or \code{"(mu) g[a]"}) name the parameters they
#' label, with and without \code{component}. The random
#' component supports point, interval, and directional hypotheses for
#' semantic standard deviation, variance, correlation, and allocation
#' quantities. Point-null
#' hypotheses require a direct parameter reference, a factor level, or a
#' linear combination of the levels of one factor term (for example
#' \code{"g[a] = g[b]"} or \code{"2 * g[a] = 0.1"}). Point hypotheses on a
#' random-effect variance are evaluated through its standard deviation, so
#' that both give the same Bayes factor. A statement comparing such a point
#' with a region of the variance (for example
#' \code{"tau2 = 0.09 vs tau2 > 0.09"}) evaluates the point through the
#' standard deviation and the region on the variance draws. Certified
#' \code{exp(affine)} fitted-scale hypotheses are available with KDE only for
#' atom-free, unconditional scalar targets. Values at a prior point mass (for
#' example, gated components at 0) have no Savage-Dickey Bayes factor; the
#' Component Inclusion table compares the exclusion and inclusion of gated
#' components. Point nulls at an exact support boundary (for example, a
#' variance proportion at 0) use the one-sided prior ordinate when BayesTools
#' classifies it as exact, finite, and positive. Publication-bias parameters
#' are not supported (class \code{RoBMA_hypothesis_target}), including the
#' quantities of the publication-bias prior such as the weight-function
#' coordinates \code{"omega[0,0.025]"} and \code{"bias_indicator"}. Other
#' quantities of the fitted model that are no tested model parameter, such as
#' inclusion indicators (\code{"mu_x_indicator"}), the component inclusions
#' of variance allocations (\code{"inclusion(study)"}), or latent cluster
#' effects, are refused with the class \code{RoBMA_hypothesis_target} as
#' well; [hypothesis_quantities()] lists the tested quantities.
#' @param standardized_coefficients whether moderator and scale coefficients
#' are tested on the standardized predictor scale. Defaults to \code{FALSE}.
#' @param conditional whether to use the conditional posterior for product-space
#' model-averaging objects. Defaults to \code{FALSE}. When this argument is
#' omitted and the selected parameter has both null and alternative components,
#' a warning notes that the full ensemble is used. Pass \code{FALSE} explicitly
#' to retain that test without the warning, or \code{TRUE} to test only models
#' where the parameter is active. Point hypotheses on linear combinations of
#' factor levels of model-averaged objects require \code{conditional = TRUE}:
#' the levels share the null component of the averaged prior.
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
#' likelihood-aware posterior ordinates for point-null hypotheses and
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
#' prior density stored by BayesTools. qCMDE/IWMDE supports exact linear maps:
#' coefficients, factor levels of every contrast, and linear combinations of
#' levels. For nonlinear joint maps, such as an exponentiated intercept that
#' also depends on varying slopes, KDE point and directional hypotheses are
#' available only when structural metadata certifies an atom-free,
#' unconditional scalar target and the point is inside the open transformed
#' support. Alternatively, use \code{standardized_coefficients = TRUE}.
#'
#' Every statement is planned before it is evaluated: the plan records the
#' tested target (its weights on the fitted coefficients, its BayesTools
#' prior density and declared atoms), the classification of each point value
#' by \code{BayesTools::prior_ordinate_status()}, and whether each density
#' method can evaluate it; [hypothesis_quantities()] renders the same plans.
#' A point hypothesis whose prior ordinate is not exact, finite, and regular
#' stops with the BayesTools condition class (for example
#' \code{BayesTools_point_mass_at_null}, \code{BayesTools_infinite_ordinate},
#' or \code{BayesTools_inexact_ordinate}). Other refusals of a test that
#' cannot be computed have the class \code{RoBMA_hypothesis_unavailable} and
#' one of \code{RoBMA_hypothesis_fixed} (the quantity is fixed by the fitted
#' model), \code{RoBMA_hypothesis_target} (no supported target or prior
#' density), or \code{RoBMA_hypothesis_method} (the density method is
#' unavailable for the target); linear combinations that BayesTools cannot
#' certify keep the class \code{BayesTools_linear_target_unavailable}. A
#' method refusal because the fitted model does not support qCMDE/IWMDE also
#' has the classes of the refusals of these requests by [plot.brma()] and
#' [marginal_means()]:
#' \code{RoBMA_density_method_unavailable} and one class naming the cause,
#' \code{RoBMA_density_method_random_unknown_v} (\code{brma.mv()}
#' random-formula models without known \code{V}),
#' \code{RoBMA_density_method_scale_components} (component-specific scale
#' formulas), or \code{RoBMA_density_method_glmm} (IWMDE for binomial and
#' Poisson GLMMs). qCMDE/IWMDE refusals for a cause for which [plot.brma()]
#' refuses these methods have its classes too:
#' \code{RoBMA_density_method_conditional_random} (conditional random-effect
#' statements), \code{RoBMA_density_method_random_target} (random-effect
#' quantities without a supported scalar random-component coordinate), and
#' \code{RoBMA_density_method_original_scale} (original-scale coefficients
#' whose fitted map is nonlinear, such as \code{exp(affine)} targets). As in
#' [plot.brma()], the capability of the fitted model is checked before these
#' causes and before other qCMDE/IWMDE refusals of the tested quantity: a
#' request to which both apply stops with the capability refusal.
#'
#' Statements that need to be restated, or an argument that needs to change,
#' stop with the class \code{RoBMA_hypothesis_statement} without
#' \code{RoBMA_hypothesis_unavailable}: point expressions that are no direct
#' reference (for example \code{"2 * mu = 0"}), nonlinear expressions of
#' factor levels, several point values of a random-effect variance in one
#' statement, model-averaged marginal-means statements that mix point and
#' region events or span several levels, ambiguous references
#' (\code{RoBMA_hypothesis_ambiguous}, see \code{component}), references to
#' parameters of another component (\code{RoBMA_component_mismatch}), factor
#' contrast coefficients such as \code{"g{1}"} (state the hypothesis on
#' factor levels; a selector of a coefficient that is a level, such as
#' \code{"g{1}"} of a treatment factor, also has the classes of the
#' BayesTools refusal \code{BayesTools_selector_unavailable}, which names the
#' level), references to unknown names or factor levels, and
#' statements that reference no parameter (followed by the classes of the
#' BayesTools refusal of such a statement,
#' \code{BayesTools_hypothesis_no_parameters} and
#' \code{BayesTools_parameter_resolution_error}, on fitted objects and on
#' marginal means). An unknown name or
#' level also has the classes
#' \code{BayesTools_parameter_not_found} and
#' \code{BayesTools_parameter_resolution_error}, with which BayesTools
#' refuses it (on fitted objects and in linear combinations of
#' marginal-means levels, the BayesTools condition itself with its fields
#' \code{alias} and \code{available}). Every other refusal of a statement's
#' references by BayesTools (class
#' \code{BayesTools_parameter_resolution_error}) stops \code{hypothesis()}
#' too, on fitted objects and on marginal means, as the BayesTools condition
#' with its classes and fields after \code{RoBMA_hypothesis_statement}, also
#' when BayesTools raises it only when it evaluates the statement (for
#' example the whole factor term next to one of its levels,
#' \code{"g[a] > mu_g"}, as an unknown quantity).
#' Missing or unsupported fitted
#' metadata of the tested target (a parameter catalog without RoBMA
#' parameters, the coefficient transform, the fitted coordinates, or the
#' linear weights of a factor level; a fit of an older RoBMA/BayesTools build)
#' stop with the class \code{RoBMA_refit_required}: refit the model. Its parent
#' \code{BayesTools_refit_required} is the class of every error of RoBMA and
#' BayesTools that asks for a refit; BayesTools' own refit errors stop
#' \code{hypothesis()} unconverted.
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

  # Every statement is planned before any is evaluated; a statement the
  # requested method cannot evaluate stops with its plan's refusal. The
  # planning and the evaluation have one catch point for BayesTools'
  # refusals of the statements' references, which stop as statement errors.
  .hypothesis_catch_statement_resolution({
    hypothesis <- .hypothesis_brma_ast(
      hypothesis = hypothesis,
      catalog    = parameter_metadata[["catalog"]]
    )
    plans <- .hypothesis_plans(
      object       = object,
      hypothesis   = hypothesis,
      component    = component,
      standardized = standardized_coefficients,
      conditional  = conditional,
      metadata     = parameter_metadata,
      n_samples    = n_samples
    )
    for (plan in plans) {
      .hypothesis_plan_check(plan, density_method)
    }

    keys   <- vapply(plans, `[[`, character(1), "group")
    groups <- lapply(unique(keys), function(key) which(keys == key))
    results <- lapply(groups, function(rows) {
      .hypothesis_plan_execute(
        plans               = plans[rows],
        object              = object,
        conditional_omitted = conditional_omitted,
        logBF               = logBF,
        BF01                = BF01,
        seed                = seed,
        density_method      = density_method,
        density_control     = density_control,
        columns             = columns
      )
    })
    if (length(results) == 1L) {
      results[[1L]]
    } else {
      .hypothesis_brma_bind_parameter_results(
        results    = results,
        groups     = groups,
        hypothesis = hypothesis
      )
    }
  })
}


# The statements of the plans of one group as one hypothesis AST: the
# statements as written ('field' "statement") or rewritten onto the tested
# parameter ("hypothesis").
.hypothesis_plan_group_ast <- function(plans, field) {

  ast <- plans[[1L]][[field]]
  ast[["statements"]] <- unlist(lapply(plans, function(plan) {
    plan[[field]][["statements"]]
  }), recursive = FALSE)

  ast
}


# Evaluates the statements of one group of plans, which share one route and
# one tested target, and restores the labels the statements were written in.
.hypothesis_plan_execute <- function(plans, object, conditional_omitted, logBF,
                                     BF01, seed, density_method,
                                     density_control, columns) {

  plan       <- plans[[1L]]
  display    <- .hypothesis_plan_group_ast(plans, "statement")
  hypothesis <- .hypothesis_plan_group_ast(plans, "hypothesis")
  if (conditional_omitted) {
    .hypothesis_brma_warn_model_averaged(object, plan[["parameter"]])
  }
  arguments <- list(
    plans           = plans,
    hypothesis      = hypothesis,
    object          = object,
    logBF           = logBF,
    BF01            = BF01,
    seed            = seed,
    density_method  = density_method,
    density_control = density_control,
    columns         = columns
  )
  out <- switch(
    plan[["route"]],
    scalar      = do.call(.hypothesis_plan_execute_scalar, arguments),
    levels      = do.call(.hypothesis_plan_execute_levels, arguments),
    combination = do.call(.hypothesis_plan_execute_combination, arguments),
    random      = do.call(.hypothesis_plan_execute_random, arguments),
    stop("Internal error: unknown hypothesis route.", call. = FALSE)
  )

  .hypothesis_brma_restore_hypothesis_labels(
    out             = out,
    hypothesis      = display,
    parameter_label = plan[["label"]]
  )
}


# A model-averaged coefficient test on the full ensemble includes the models
# in which the coefficient is fixed by its null component.
.hypothesis_brma_warn_model_averaged <- function(object, parameter) {

  prior_list <- attr(object[["fit"]], "prior_list", exact = TRUE)
  if (!.is_RoBMA(object) || !parameter %in% names(prior_list) ||
      !BayesTools::is.prior.mixture(prior_list[[parameter]])) {
    return(invisible(FALSE))
  }
  prior_components <- attr(prior_list[[parameter]], "components", exact = TRUE)
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

  invisible(TRUE)
}


# The method a group is evaluated with: qCMDE/IWMDE ordinates only for point
# statements (region statements are method-free and use the BayesTools
# route).
.hypothesis_plan_execution_method <- function(plans, density_method) {

  point <- any(vapply(plans, `[[`, logical(1), "point"))
  if (.density_method_uses_precomputed(density_method, allow_normal = TRUE)) {
    return(if (point) "precomputed" else "KDE")
  }

  density_method
}


.hypothesis_plan_targets <- function(plans) {

  unlist(lapply(plans, `[[`, "targets"), recursive = FALSE)
}


.hypothesis_plan_density_control <- function(density_control) {

  if (is.null(density_control[["normalization_points"]])) {
    density_control[["normalization_points"]] <- max(
      50L,
      density_control[["n_points"]]
    )
  }

  density_control
}


# Attaches the qCMDE/IWMDE ordinates of scalar point values to a posterior.
.hypothesis_plan_attach_scalar <- function(object, posterior, parameter,
                                           parameter_label, values, spec,
                                           conditional, density_method,
                                           density_control,
                                           display_transform = NULL) {

  density_control <- .hypothesis_plan_density_control(density_control)
  .iwmde_check_point_ordinate_supported(object, density_method)
  posterior <- .hypothesis_brma_keep_requested_ordinates(
    posterior  = posterior,
    point_refs = data.frame(level = NA_character_, value = values)
  )
  context <- .iwmde_context(object, density_control[["integration_control"]])

  .hypothesis_brma_attach_iwmde_scalar(
    posterior            = posterior,
    raw_posterior        = posterior,
    context              = context,
    estimate_cache       = .iwmde_estimate_cache(),
    parameter            = parameter,
    parameter_label      = parameter_label,
    value                = unique(values),
    conditional          = conditional,
    n_points             = density_control[["n_points"]],
    samples              = density_control[["samples"]],
    target_relative_mcse = density_control[["target_relative_mcse"]],
    normalization_points = density_control[["normalization_points"]],
    normalization_prob   = density_control[["normalization_prob"]],
    integration_control  = density_control[["integration_control"]],
    density_method       = density_method,
    parameter_spec       = spec,
    display_transform    = display_transform
  )
}


.hypothesis_plan_execute_scalar <- function(plans, hypothesis, object, logBF,
                                            BF01, seed, density_method,
                                            density_control, columns) {

  plan      <- plans[[1L]]
  posterior <- plan[["draws"]][["posterior"]]
  method    <- .hypothesis_plan_execution_method(plans, density_method)
  if (identical(method, "precomputed")) {
    targets   <- .hypothesis_plan_targets(plans)
    posterior <- .hypothesis_plan_attach_scalar(
      object          = object,
      posterior       = posterior,
      parameter       = plan[["parameter"]],
      parameter_label = plan[["label"]],
      values          = vapply(targets, `[[`, numeric(1), "value"),
      spec            = targets[[1L]][["spec"]],
      conditional     = if (plan[["conditional"]]) plan[["parameter"]],
      density_method  = density_method,
      density_control = density_control
    )
  }

  out <- BayesTools::hypothesis_BF(
    posterior      = posterior,
    hypothesis     = hypothesis,
    parameter      = plan[["parameter"]],
    logBF          = logBF,
    BF01           = BF01,
    seed           = seed,
    columns        = columns,
    density_method = method
  )
  if (identical(method, "precomputed")) {
    out <- .hypothesis_brma_append_iwmde_warnings(
      table     = out,
      posterior = posterior
    )
  }

  out
}


# Statements on the marginal posterior of a factor term: each point value of
# a level gets the level's planned prior density and, for qCMDE/IWMDE, its
# ordinate of the level's linear target.
.hypothesis_plan_execute_levels <- function(plans, hypothesis, object, logBF,
                                            BF01, seed, density_method,
                                            density_control, columns) {

  plan      <- plans[[1L]]
  posterior <- plan[["draws"]][["posterior"]]
  targets   <- .hypothesis_plan_targets(plans)
  for (target in targets) {
    BayesTools::posterior_metadata(posterior[[target[["level"]]]], "prior_density") <-
      target[["prior_density"]]
  }
  method <- .hypothesis_plan_execution_method(plans, density_method)
  if (identical(method, "precomputed")) {
    density_control <- .hypothesis_plan_density_control(density_control)
    .iwmde_check_point_ordinate_supported(object, density_method)
    refs <- unique(data.frame(
      level = vapply(targets, function(target) target[["level"]], character(1)),
      value = vapply(targets, `[[`, numeric(1), "value"),
      stringsAsFactors = FALSE
    ))
    posterior <- .hypothesis_brma_keep_requested_ordinates(posterior, refs)
    context        <- .iwmde_context(object, density_control[["integration_control"]])
    estimate_cache <- .iwmde_estimate_cache()
    for (i in seq_len(nrow(refs))) {
      target <- targets[[which(
        vapply(targets, function(target) target[["level"]], character(1)) ==
          refs[["level"]][[i]]
      )[[1L]]]]
      posterior <- .hypothesis_brma_attach_iwmde_level(
        posterior            = posterior,
        raw_posterior        = posterior,
        context              = context,
        estimate_cache       = estimate_cache,
        parameter            = plan[["parameter"]],
        level                = refs[["level"]][[i]],
        value                = refs[["value"]][[i]],
        conditional          = if (plan[["conditional"]]) plan[["parameter"]],
        n_points             = density_control[["n_points"]],
        samples              = density_control[["samples"]],
        target_relative_mcse = density_control[["target_relative_mcse"]],
        normalization_points = density_control[["normalization_points"]],
        normalization_prob   = density_control[["normalization_prob"]],
        integration_control  = density_control[["integration_control"]],
        density_method       = density_method,
        parameter_spec       = target[["spec"]]
      )
    }
  }

  out <- BayesTools::hypothesis_BF(
    posterior      = posterior,
    hypothesis     = hypothesis,
    parameter      = plan[["parameter"]],
    logBF          = logBF,
    BF01           = BF01,
    seed           = seed,
    columns        = columns,
    density_method = method
  )
  if (identical(method, "precomputed")) {
    out <- .hypothesis_brma_append_iwmde_warnings(
      table     = out,
      posterior = posterior
    )
  }

  out
}


# Statements on one linear combination of the levels of a factor term,
# evaluated on the scalar linear target (level contrasts, derived linear
# expressions).
.hypothesis_plan_execute_combination <- function(plans, hypothesis, object,
                                                 logBF, BF01, seed,
                                                 density_method,
                                                 density_control, columns) {

  plan   <- plans[[1L]]
  target <- if (length(plans) == 1L) {
    plan[["linear_target"]]
  } else {
    combined <- BayesTools::hypothesis_linear_target(
      posterior  = plan[["draws"]][["posterior"]],
      hypothesis = hypothesis,
      parameter  = plan[["parameter"]]
    )
    BayesTools::posterior_metadata(combined[["posterior"]], "prior_density") <-
      BayesTools::posterior_metadata(plan[["linear_target"]][["posterior"]], "prior_density")
    combined
  }
  method <- .hypothesis_plan_execution_method(plans, density_method)
  if (identical(method, "precomputed")) {
    targets <- .hypothesis_plan_targets(plans)
    target[["posterior"]] <- .hypothesis_plan_attach_scalar(
      object          = object,
      posterior       = target[["posterior"]],
      parameter       = target[["parameter"]],
      parameter_label = plan[["label"]],
      values          = vapply(targets, `[[`, numeric(1), "value"),
      spec            = targets[[1L]][["spec"]],
      conditional     = .iwmde_first_nonempty_condition(
        BayesTools::posterior_metadata(target[["posterior"]], "condition"),
        c("effective_conditional", "conditional")
      ),
      density_method  = density_method,
      density_control = density_control
    )
  }

  out <- BayesTools::hypothesis_BF(
    posterior      = target[["posterior"]],
    hypothesis     = target[["hypothesis"]],
    parameter      = target[["parameter"]],
    logBF          = logBF,
    BF01           = BF01,
    seed           = seed,
    columns        = columns,
    density_method = method
  )
  if (identical(method, "precomputed")) {
    out <- .hypothesis_brma_append_iwmde_warnings(
      table     = out,
      posterior = target[["posterior"]]
    )
  }
  # One linear target: row i is statement i.
  .hypothesis_brma_set_row_names(
    out       = out,
    row_names = .hypothesis_brma_row_names(
      labels     = rep(plan[["parameter"]], nrow(out)),
      statements = seq_len(nrow(out))
    )
  )
}


# Statements on a random-effect quantity, evaluated on its BayesTools mixed
# posterior. Point statements on a variance are evaluated through its
# standard deviation (square-root values); the table then shows the variance
# statements and densities. Statements comparing a variance point with a
# region evaluate the point part through the standard deviation and the
# region part on the variance draws.
.hypothesis_plan_execute_random <- function(plans, hypothesis, object, logBF,
                                            BF01, seed, density_method,
                                            density_control, columns) {

  plan      <- plans[[1L]]
  parameter <- plan[["parameter"]]
  samples   <- plan[["draws"]][["samples"]]
  selected  <- plan[["random"]]
  defined_footnote <- .brma_random_parameter_defined_footnote(
    label             = selected[["spec"]][["label"]],
    samples           = samples[[parameter]],
    posterior_defined = .brma_random_parameter_defined_share(
      object,
      samples[[parameter]]
    )
  )
  method <- .hypothesis_plan_execution_method(plans, density_method)

  if (isTRUE(plan[["prior_draws"]])) {
    # Without a prior density with deterministic provenance, region
    # hypotheses take the prior probabilities from prior draws (the plans
    # refused point hypotheses).
    return(.hypothesis_brma_random_prior_draws(
      object                    = object,
      parameter                 = parameter,
      selected                  = selected,
      posterior                 = samples[[parameter]],
      hypothesis                = hypothesis,
      standardized_coefficients = plan[["standardized"]],
      logBF                     = logBF,
      BF01                      = BF01,
      seed                      = seed,
      n_samples                 = plan[["n_samples"]],
      columns                   = columns,
      density_method            = method
    ))
  }

  arguments <- list(
    plan            = plan,
    object          = object,
    logBF           = logBF,
    BF01            = BF01,
    seed            = seed,
    method          = method,
    density_method  = density_method,
    density_control = density_control,
    columns         = columns
  )
  values <- vapply(.hypothesis_plan_targets(plans), `[[`, numeric(1), "value")
  out <- if (is.null(plan[["evaluation"]][["parameter"]])) {
    do.call(.hypothesis_plan_random_evaluate, c(arguments, list(
      posterior  = samples[[parameter]],
      parameter  = parameter,
      hypothesis = hypothesis,
      values     = values
    )))
  } else if (plan[["region"]]) {
    do.call(.hypothesis_plan_random_sd_region, c(arguments, list(
      posterior  = samples[[parameter]],
      hypothesis = hypothesis,
      values     = values
    )))
  } else {
    do.call(.hypothesis_plan_random_sd_points, c(arguments, list(
      hypothesis = hypothesis,
      values     = values
    )))
  }
  if (!is.null(defined_footnote)) {
    attr(out, "footnotes") <- c(attr(out, "footnotes"), defined_footnote)
  }

  out
}


# Statements on the marginal posterior 'posterior' of the random-effect
# quantity 'parameter'. For qCMDE/IWMDE, the ordinates of the point values
# 'values' of that quantity are attached first.
.hypothesis_plan_random_evaluate <- function(plan, object, posterior, parameter,
                                             hypothesis, values, logBF, BF01,
                                             seed, method, density_method,
                                             density_control, columns) {

  if (identical(method, "precomputed")) {
    target    <- plan[["density_target"]]
    posterior <- .hypothesis_plan_attach_scalar(
      object            = object,
      posterior         = posterior,
      parameter         = target[["parameter"]],
      parameter_label   = plan[["label"]],
      values            = values,
      spec              = target[["parameter_spec"]],
      conditional       = NULL,
      density_method    = density_method,
      density_control   = density_control,
      display_transform = target[["display_transform"]]
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
    density_method = method
  )
  if (identical(method, "precomputed")) {
    out <- .hypothesis_brma_append_iwmde_warnings(
      table     = out,
      posterior = posterior,
      parameter = parameter
    )
  }

  out
}


# Point statements on a variance, evaluated through its standard deviation:
# the table shows the variance statements, and the prior and posterior
# densities on the variance scale.
.hypothesis_plan_random_sd_points <- function(plan, hypothesis, values, ...) {

  evaluation <- plan[["evaluation"]]
  out <- .hypothesis_plan_random_evaluate(
    plan       = plan,
    posterior  = evaluation[["samples"]][[evaluation[["parameter"]]]],
    parameter  = evaluation[["parameter"]],
    hypothesis = .hypothesis_plan_sd_hypothesis(hypothesis, evaluation[["parameter"]]),
    values     = sqrt(values),
    ...
  )

  .hypothesis_plan_sd_restore(out, hypothesis)
}


# Statements comparing a variance point with a region of the variance, e.g.
# 'tau2 = 0.09 vs tau2 > 0.09' (one point value per statement, 'values').
# Their Bayes factor is the point part, the point against the encompassing
# model ('tau2 = 0.09 vs tau2 != 0.09', evaluated through the standard
# deviation as a point statement), over the region part, the region against
# the encompassing model (evaluated on the variance draws); the inverse with
# the region on the left. BayesTools evaluates the statements on the variance
# with a unit point part (a posterior ordinate equal to the prior ordinate,
# carrying the Monte Carlo error of the point part, which BayesTools combines
# with that of the region part); the point part then multiplies.
.hypothesis_plan_random_sd_region <- function(plan, object, posterior,
                                              hypothesis, values, logBF, BF01,
                                              seed, method, density_method,
                                              density_control, columns) {

  if (length(values) != length(hypothesis[["statements"]])) {
    stop("Internal error: variance point-region statements are misaligned.",
         call. = FALSE)
  }
  parameter <- plan[["parameter"]]
  points    <- unique(values)
  point     <- .hypothesis_plan_random_sd_points(
    plan            = plan,
    hypothesis      = BayesTools::hypothesis_rewrite(
      BayesTools::hypothesis_parse(sprintf(
        "theta = %.17g vs theta != %.17g", points, points
      )),
      c(theta = parameter)
    ),
    values          = points,
    object          = object,
    logBF           = FALSE,
    BF01            = FALSE,
    seed            = seed,
    method          = method,
    density_method  = density_method,
    density_control = density_control,
    columns         = "default"
  )
  unit_ordinate <- exp(vapply(points, function(value) {
    BayesTools::prior_density_ordinate(
      plan[["prior_density"]],
      value
    )[["log_density"]]
  }, numeric(1)))
  # The unit part carries the conditioning of the draws, so that BayesTools
  # matches it to conditional (product-space) draws.
  BayesTools::posterior_metadata(posterior, "posterior_ordinate") <- do.call(
    BayesTools::posterior_ordinate_attribute,
    c(
      list(
        value          = points,
        ordinate       = unit_ordinate,
        method         = "unit point part",
        density_method = density_method,
        diagnostics    = list(BF_error_percent = as.numeric(point[["BF_error"]]))
      ),
      .iwmde_sample_condition_metadata(posterior)
    )
  )

  out <- BayesTools::hypothesis_BF(
    posterior      = posterior,
    hypothesis     = hypothesis,
    parameter      = parameter,
    logBF          = logBF,
    BF01           = BF01,
    seed           = seed,
    columns        = columns,
    density_method = "precomputed"
  )
  point_index <- match(values, points)
  point_BF    <- attr(point, "raw_BF", exact = TRUE)[point_index]
  point_left  <- vapply(hypothesis[["statements"]], function(statement) {
    identical(statement[["left"]][["type"]], "point")
  }, logical(1))
  raw_BF <- attr(out, "raw_BF", exact = TRUE)
  raw_BF <- ifelse(point_left, raw_BF * point_BF, raw_BF / point_BF)
  attr(out, "raw_BF") <- raw_BF
  if ("BF" %in% names(out)) {
    out[["BF"]] <- BayesTools::format_BF(raw_BF, logBF = logBF, BF01 = BF01)
  }

  # The warnings and density diagnostics of the point parts, on the rows of
  # their statements.
  point_warnings <- attr(point, "warnings", exact = TRUE)
  point_rows     <- match(names(point_warnings), rownames(point))
  if (length(point_rows) != length(point_warnings)) {
    point_rows <- rep(NA_integer_, length(point_warnings))
  }
  row_warnings <- unlist(lapply(seq_along(values), function(i) {
    warnings <- as.character(point_warnings[point_rows %in% point_index[[i]]])
    stats::setNames(warnings, rep(rownames(out)[[i]], length(warnings)))
  }), use.names = TRUE)
  attr(out, "warnings") <- .hypothesis_brma_unique_named_warnings(c(
    attr(out, "warnings", exact = TRUE),
    row_warnings,
    point_warnings[is.na(point_rows)]
  ))
  attr(out, "density_diagnostics") <- attr(point, "density_diagnostics", exact = TRUE)

  out
}


# Point statements on a variance written on its standard deviation 'sd':
# each point side's value is replaced by its square root.
.hypothesis_plan_sd_hypothesis <- function(hypothesis, sd) {

  statements <- vapply(hypothesis[["statements"]], function(statement) {

    sides <- vapply(c("left", "right"), function(side_name) {
      side <- statement[[side_name]]
      operator <- switch(
        side[["type"]],
        point     = "=",
        not_point = "!=",
        stop("Internal error: expected a point statement on a variance.",
             call. = FALSE)
      )
      paste("theta", operator, sprintf("%.17g", sqrt(side[["value"]])))
    }, character(1))
    if (isTRUE(statement[["explicit"]])) {
      paste(sides[[1L]], "vs", sides[[2L]])
    } else {
      sides[[1L]]
    }
  }, character(1))

  BayesTools::hypothesis_rewrite(
    BayesTools::hypothesis_parse(unname(statements)),
    c(theta = sd)
  )
}


# The table of variance point statements evaluated through the standard
# deviation: the variance statements, and the prior and posterior densities
# on the variance scale (divided by the derivative 2 * sd of the square).
.hypothesis_plan_sd_restore <- function(out, hypothesis) {

  values <- vapply(hypothesis[["statements"]], function(statement) {
    statement[["left"]][["value"]]
  }, numeric(1))
  if (nrow(out) != length(values)) {
    stop("Internal error: variance point rows are misaligned.", call. = FALSE)
  }
  for (column in intersect(c("prior", "posterior"), names(out))) {
    out[[column]] <- out[[column]] / (2 * sqrt(values))
  }
  attr(out, "hypothesis_ast") <- hypothesis

  out
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


# Region hypotheses on a random-effect quantity without a canonical prior
# density: prior probabilities from the quantity's prior draws.
.hypothesis_brma_random_prior_draws <- function(
    object, parameter, selected, posterior, hypothesis,
    standardized_coefficients, logBF, BF01, seed, n_samples,
    columns, density_method) {

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
  footnote <- .brma_random_parameter_defined_footnote(
    label             = selected[["spec"]][["label"]],
    samples           = prior[["samples"]],
    posterior_defined = .brma_random_parameter_defined_share(object, posterior),
    prior_defined     = prior_defined
  )
  if (!is.null(footnote)) {
    attr(out, "footnotes") <- c(attr(out, "footnotes"), footnote)
  }

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
    .stop_refit_required(
      "Formula coefficient transformation metadata are unsupported. Refit ",
      "the model with the current BayesTools version."
    )
  }
  target <- selected[["parameter"]]
  target_i <- match(target, transform[["target_names"]])
  if (is.na(target_i)) {
    .stop_refit_required(
      "Resolved formula coefficient '", target,
      "' is absent from the fitted coefficient transformation."
    )
  }

  return(list(
    formula_parameter = formula_parameter,
    target            = target,
    target_i          = target_i,
    transform         = transform
  ))
}


# The fitted prior-list entries (formula terms) that own the coordinates a
# formula coefficient target weights, from the coordinate table and the
# formula name maps.
.hypothesis_brma_target_prior_parameters <- function(object, weights) {

  coordinates <- BayesTools::parameter_coordinates(object[["fit"]])
  rows <- match(names(weights), coordinates[["coordinate_name"]])
  if (anyNA(rows)) {
    .stop_refit_required(
      "Fitted coordinate metadata of the hypothesis target are unavailable. ",
      "Refit the model with the current RoBMA/BayesTools build."
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
