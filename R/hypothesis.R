# ============================================================================ #
# hypothesis.R
# ============================================================================ #
#
# User-facing hypothesis Bayes factors for brma objects and their marginal
# means. The statements are planned and evaluated by the hypothesis engine of
# BayesTools (BayesTools::hypothesis_BF_model()); RoBMA describes its models to
# it in hypothesis-provider.R and provides the likelihood-aware posterior
# ordinates qCMDE and IWMDE in hypothesis-backends.R.
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
#' For half-bounded targets, IWMDE fits the distance from the declared boundary
#' and normalizes a Normal weight over its central 95 percent positive-tail
#' interval, using the actual computed endpoints. A draw outside that compact
#' interval contributes zero while retaining its original row and chain. The
#' same compact family is used when conditional fitting is unavailable.
#' Compact focal support does not certify a finite joint second moment or
#' MCMC long-run variance; reported uncertainty is conditional on the fitted
#' proposal and does not include uncertainty from learning it on the same draws.
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
#' Selected weightfunction models near zero cumulative ordered effects can take
#' more than 10 minutes with qCMDE or IWMDE.
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
#' @param ... additional arguments passed to methods. For fitted \code{brma}
#' objects, unused arguments trigger a warning.
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
hypothesis.default <- function(object, hypothesis, ...) {

  BayesTools::hypothesis_BF(posterior = object, hypothesis = hypothesis, ...)
}


#' @rdname hypothesis
#'
#' @param hypothesis character vector with scalar hypothesis statements. For
#' \code{marginal_means.brma} objects, constants are specified on the fitted
#' linear-predictor scale, even when the object stores a display transformation.
#' Separate statements may reference different model parameters, also the
#' marginal means of different parameters of \code{marginal_means.brma}
#' objects (unless \code{parameter} is set); each statement itself must
#' resolve to one scalar model target.
#' @param component parameter component. Defaults to \code{"auto"}, which
#' infers the component when possible. Use \code{"mods"} (alias
#' \code{"location"}), \code{"scale"}, or \code{"random"} to disambiguate
#' terms used in multiple model components. A statement whose references
#' still name several model parameters stops with an error of class
#' \code{BayesTools_hypothesis_ambiguous}, with the parent class
#' \code{BayesTools_hypothesis_statement}: the statement or \code{component} needs
#' to change, the test itself is not unavailable. An ambiguous reference that
#' BayesTools refuses with the class \code{BayesTools_parameter_ambiguous}
#' (for example the contrast coefficient \code{"g{1}"} of a factor term in
#' the location and the scale formula) has these classes too, before its
#' BayesTools classes. A statement whose parameters belong to another
#' component than \code{component}, including a name that is known only
#' outside \code{component} next to references of \code{component} (with or
#' without factor levels, e.g. \code{"log_tau_g[a] > mu_g[a]"} with
#' \code{component = "mods"}), stops with the classes
#' \code{BayesTools_hypothesis_statement} and
#' \code{BayesTools_hypothesis_component_mismatch} ([plot.brma()] stops with
#' \code{RoBMA_component_mismatch} for its own component mismatches). The
#' exception is
#' the alias of a factor term of both the location and the scale formula
#' (for example \code{g} for \code{mods = ~ g, scale = ~ g}) next to a name
#' that is unknown within \code{component}, such as a name of the other
#' component or a display label (for example \code{"g > log_tau_intercept"}
#' with \code{component = "mods"}): the statement is resolved without
#' \code{component}, where the alias is ambiguous, and stops with
#' \code{BayesTools_hypothesis_statement} followed by the classes of the
#' BayesTools refusal of its first unresolved reference, the ambiguity of the
#' alias (\code{BayesTools_parameter_ambiguous}, without
#' \code{BayesTools_hypothesis_ambiguous}, as \code{component} is set already) or
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
#' that both give the same Bayes factor at positive variance values. At zero,
#' the square map is singular: the variance is tested on its own scale and
#' follows its own prior-ordinate eligibility rules. A statement comparing such a point
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
#' are not supported (class \code{BayesTools_hypothesis_target}), including the
#' quantities of the publication-bias prior such as the weight-function
#' coordinates \code{"omega[0,0.025]"} and \code{"bias_indicator"}. Other
#' quantities of the fitted model that are no tested model parameter, such as
#' inclusion indicators (\code{"mu_x_indicator"}), the component inclusions
#' of variance allocations (\code{"inclusion(study)"}), or latent cluster
#' effects, are refused with the class \code{BayesTools_hypothesis_target} as
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
#' The statements are planned and evaluated by the hypothesis engine of
#' BayesTools (\code{BayesTools::hypothesis_BF_model()}); RoBMA describes its
#' models to it (\code{BayesTools::hypothesis_provider()}) and provides the
#' density methods qCMDE and IWMDE (\code{BayesTools::hypothesis_density_backend()}).
#' Every statement is planned before it is evaluated: the plan records the
#' tested target (its weights on the fitted coefficients, its BayesTools
#' prior density and declared atoms), the classification of each point value
#' by \code{BayesTools::prior_ordinate_status()}, and whether each density
#' method can evaluate it; [hypothesis_quantities()] renders the same plans.
#' A point hypothesis whose prior ordinate is not exact, finite, and regular
#' stops with the BayesTools condition class (for example
#' \code{BayesTools_point_mass_at_null}, \code{BayesTools_infinite_ordinate},
#' or \code{BayesTools_inexact_ordinate}). Other refusals of a test that
#' cannot be computed have the class \code{BayesTools_hypothesis_unavailable} and
#' one of \code{BayesTools_hypothesis_fixed} (the quantity is fixed by the fitted
#' model), \code{BayesTools_hypothesis_target} (no supported target or prior
#' density), or \code{BayesTools_hypothesis_method} (the density method is
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
#' [plot.brma()], allocated component SDs in a predictor-scaled block require
#' `standardized_coefficients = TRUE` or `density_method = 'KDE'`; the
#' block-wide restriction also covers an unchanged centered intercept. As in
#' [plot.brma()], the capability of the fitted model is checked before these
#' causes and before other qCMDE/IWMDE refusals of the tested quantity: a
#' request to which both apply stops with the capability refusal. A qCMDE/IWMDE
#' ordinate that cannot be constructed stops with
#' \code{RoBMA_density_ordinate_error}; when a fitted chain has no eligible row
#' of the target, the error also has the class
#' \code{RoBMA_density_chain_coverage} and its message names
#' \code{density_method = "KDE"} and longer chains as remedies. The refusal for
#' fewer than 20 continuous active samples of the target takes precedence over
#' this diagnosis: such a target is refused for that reason alone, without the
#' class \code{RoBMA_density_chain_coverage}, even when a fitted chain has none
#' of its rows.
#'
#' Statements that need to be restated, or an argument that needs to change,
#' stop with the class \code{BayesTools_hypothesis_statement} without
#' \code{BayesTools_hypothesis_unavailable}: point expressions that are no direct
#' reference (for example \code{"2 * mu = 0"}), nonlinear expressions of
#' factor levels, several point values of a random-effect variance in one
#' statement, one quantity named by several names in one statement (for
#' example \code{"mu > 0 & intercept < 1"} or \code{"g[a] > mu_g[a]"}; use
#' one name), model-averaged marginal-means statements that mix point and
#' region events or span several marginal means (levels or scalar marginal
#' means), ambiguous references
#' (\code{BayesTools_hypothesis_ambiguous}, see \code{component}), references to
#' parameters of another component (\code{BayesTools_hypothesis_component_mismatch}), marginal-means
#' statements that refer to another parameter than \code{parameter}
#' (\code{BayesTools_hypothesis_parameter_mismatch}), factor
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
#' with its classes and fields after \code{BayesTools_hypothesis_statement}, also
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
#' \code{1 - normalization_prob}: zero for rows whose Gaussian likelihood
#' kernel meets a normal prior and whose evaluated values are those of the
#' resulting truncated normal (their normalizer is its closed-form mass, so no
#' lattice decides it), exactly computed for such rows that keep the lattice
#' normalizer, bounded for a Gaussian kernel with
#' another proper scalar prior, and estimated from the tails (each tail
#' extended until its estimate is at most \code{(1 - normalization_prob) / 2})
#' otherwise. The tail extension has a limit, ends at values the joint density
#' cannot evaluate, and cannot help a tail that does not decay, so the
#' estimate for such a row can exceed that target and the bound below with it.
#' A row normalizer missing the fraction \eqn{t} overstates that
#' row's density by the factor \eqn{1 / (1 - t)}, so the ordinate is
#' overstated by at most \eqn{t / (1 - t)} for the largest row truncation,
#' which is at most \code{(1 - normalization_prob) / normalization_prob}; at
#' the default \code{0.999} that is about 0.1 percent, and it is attained when
#' all rows share one conditional law. [density_diagnostics()] report the
#' largest row truncation (\code{normalization_truncation}), how it is known
#' (\code{normalization_truncation_status}), and this bound
#' (\code{truncation_ordinate_bound}), together with the change between the
#' two nested grids (\code{ordinate_relative_change}), which estimates the
#' discretization error of the lattice normalizers, and the contribution-share
#' weighted lattice-resolution estimate (\code{lattice_resolution_relative_error}).
#' The resolution estimate uses the Gaussian kernel width at the nested-grid
#' maximum, or local log curvature when no kernel is available. Neither the
#' resolution component nor the combined error estimate is a rigorous bound;
#' a sharp prior core can be narrower than its likelihood kernel. Exact
#' nonlattice normalizers have no lattice-resolution component. Set
#' \code{normalization_prob} closer to 1 to reduce the truncation, and increase
#' \code{normalization_points} to reduce the discretization. IWMDE checks its
#' normalization mass on
#' \code{normalization_points} points over the central \code{normalization_prob}
#' quantile range of the target draws, widened by 10 percent on each side.
#' A point ordinate warns at relative MCSE at least 5 percent,
#' ESS below 100, maximum contribution share at least 20 percent, or fewer than
#' 100 finite contributions. These same-sample diagnostics do not suppress a
#' finite fixed-design estimate. The qCMDE ordinate error estimate (nested-grid
#' change, lattice resolution and the tail-truncation bound) warns above
#' 2.5 percent and is
#' rejected above 5 percent; IWMDE normalization error warns above 5 percent and
#' is rejected above 10 percent. For retained-location conditional normalizers,
#' adaptive-quadrature sensitivity warns above 2.5 percent and is rejected above 5 percent.
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
  if (conditional && !.is_RoBMA(object)) {
    stop("'conditional' hypotheses are available only for RoBMA objects.",
         call. = FALSE)
  }
  # The engine would name its own argument ('namespace') in this message.
  BayesTools::check_char(component, "component", check_length = 1, allow_NULL = TRUE, allow_NA = FALSE)

  BayesTools::hypothesis_BF_model(
    object              = object,
    hypothesis          = hypothesis,
    namespace           = component,
    transform_scaled    = !standardized_coefficients,
    conditional         = conditional,
    logBF               = logBF,
    BF01                = BF01,
    seed                = seed,
    density_method      = density_method,
    density_control     = density_control,
    n_samples           = n_samples,
    columns             = columns,
    conditional_omitted = conditional_omitted
  )
}



#' @rdname hypothesis
#'
#' @param parameter optional marginal-means parameter/term used to disambiguate
#' the hypothesis expression. Without it, each statement is tested on the
#' marginal means of the parameter it names, so that separate statements may
#' refer to different parameters (for example
#' \code{c("g[a] > 0", "h[b] > 0")}); the rows are those of the statements
#' tested one at a time, in the order of the statements. One statement must
#' refer to the marginal means of one parameter: a statement that names
#' several (for example \code{"g[a] > h[b]"}), or an alias that names several
#' parameters, stops with an error of class
#' \code{BayesTools_hypothesis_ambiguous}, with the parent class
#' \code{BayesTools_hypothesis_statement}; \code{parameter} selects the parameter
#' of an ambiguous alias. With \code{parameter}, every statement is tested on
#' that parameter, and a statement that refers to the marginal means of
#' another parameter stops with an error of class
#' \code{BayesTools_hypothesis_parameter_mismatch}, with the parent class
#' \code{BayesTools_hypothesis_statement}; for a statement that names several
#' parameters (\code{"g[a] > h[b]"} with \code{parameter = "g"}) its message
#' is the one given without \code{parameter}. A statement that references no
#' marginal-means parameter stops with the class
#' \code{BayesTools_hypothesis_statement}; for a statement without parameter
#' symbols (for example \code{"1 > 0"}, also with \code{parameter}) also with
#' \code{BayesTools_hypothesis_no_parameters} and
#' \code{BayesTools_parameter_resolution_error}, and for an unknown name, or
#' an unknown level of a factor term, also with
#' \code{BayesTools_parameter_not_found} and
#' \code{BayesTools_parameter_resolution_error}, as on fitted objects.
#' @details Marginal-means hypotheses are specified on the fitted
#' linear-predictor scale. Display transformations stored by
#' \code{marginal_means()} do not transform hypothesis constants. Single-model
#' hypotheses and model-averaged region hypotheses use averaged marginal draws.
#' Model-averaged point-null hypotheses use alternative-conditioned draws to
#' avoid null atoms in averaged marginals. For model-averaged objects,
#' hypotheses that mix point and region events, and point hypotheses on more
#' than one marginal mean (several levels of a factor, or levels and scalar
#' marginal means of one parameter or of several), are not supported.
#' Separate region statements may refer to different parameters, while each
#' statement refers to the marginal means of one parameter.
#' Single-model hypotheses use one common averaged posterior and may combine
#' these events. Point hypotheses require the exact prior ordinate that
#' BayesTools attaches to the marginal mean, as for fitted objects.
#'
#' @export

hypothesis.marginal_means.brma <- function(object, hypothesis,
                                           parameter = NULL,
                                           logBF = FALSE, BF01 = FALSE,
                                           seed = NULL,
                                           density_method = NULL,
                                           density_control = NULL,
                                           columns = "default",
                                           ...) {

  BayesTools::check_char(hypothesis, "hypothesis", check_length = 0, allow_NA = FALSE)
  BayesTools::check_bool(logBF, "logBF")
  BayesTools::check_bool(BF01, "BF01")
  BayesTools::check_real(seed, "seed", check_length = 1, allow_NULL = TRUE, allow_NA = FALSE)
  BayesTools::check_char(columns, "columns", check_length = 0, allow_NA = FALSE)
  .warn_unused_dots(
    dots    = list(...),
    allowed = character(),
    caller  = "hypothesis.marginal_means()"
  )

  BayesTools::hypothesis_BF_model(
    object          = object,
    hypothesis      = hypothesis,
    parameter       = parameter,
    logBF           = logBF,
    BF01            = BF01,
    seed            = seed,
    density_method  = density_method,
    density_control = density_control,
    columns         = columns
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


#' @importFrom BayesTools hypothesis_quantities
#' @export
BayesTools::hypothesis_quantities
