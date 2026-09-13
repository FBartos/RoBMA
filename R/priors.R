#' @title Prior Distribution
#'
#' @description Create prior distribution objects used by RoBMA and brma
#' fitting functions.
#'
#' @param distribution character. Prior distribution name.
#' @param parameters list. Distribution parameters.
#' @param truncation list with \code{lower} and \code{upper} truncation bounds.
#' @param prior_weights numeric prior model weight.
#'
#' @details This is RoBMA's re-export of \code{BayesTools::prior()}; supported
#' distribution names and parameter lists follow BayesTools.
#'
#' @return An object inheriting from \code{prior}.
#'
#' @examples
#' prior("normal", list(mean = 0, sd = 0.5))
#' prior(
#'   "normal",
#'   list(mean = 0, sd = 0.5),
#'   truncation = list(lower = 0, upper = Inf)
#' )
#'
#' @export
prior <- BayesTools::prior

#' @title Empty Prior
#'
#' @description Create an empty prior used to omit a model component.
#'
#' @inheritParams prior
#'
#' @return An object inheriting from \code{prior} and \code{prior.none}.
#'
#' @export
prior_none <- BayesTools::prior_none

#' @title Factor Prior
#'
#' @description Create priors on factor contrasts.
#'
#' @inheritParams prior
#' @param contrast character. Contrast coding used for factor levels. Common
#' RoBMA model options are \code{"treatment"}, \code{"meandif"}, and
#' \code{"orthonormal"}; the standalone helper also accepts BayesTools contrast
#' aliases such as \code{"independent"}.
#'
#' @details Mean-difference and orthonormal contrasts require vector or
#' multivariate priors. Treatment/dummy and independent contrasts use univariate
#' simple priors per contrast coefficient.
#'
#' @return An object inheriting from \code{prior}, \code{prior.factor}, and a
#' contrast-specific marker class.
#'
#' @export
prior_factor <- BayesTools::prior_factor

#' @title PET Prior
#'
#' @description Create PET publication-bias regression priors.
#'
#' @inheritParams prior
#'
#' @details This forwards to \code{BayesTools::prior_PET()} and uses the same
#' distribution and parameter conventions as \code{prior()}. By default, PET
#' priors are truncated to the positive half-line.
#'
#' @return An object inheriting from \code{prior} with the \code{prior.PET}
#' marker class.
#'
#' @seealso \code{\link{publication_bias_prior_specification}}
#'
#' @export
prior_PET  <- BayesTools::prior_PET

#' @title PEESE Prior
#'
#' @description Create PEESE publication-bias regression priors.
#'
#' @inheritParams prior
#'
#' @details This forwards to \code{BayesTools::prior_PEESE()} and uses the same
#' distribution and parameter conventions as \code{prior()}.
#'
#' @return An object inheriting from \code{prior} with the \code{prior.PEESE}
#' marker class.
#'
#' @seealso \code{\link{publication_bias_prior_specification}}
#'
#' @export
prior_PEESE <- BayesTools::prior_PEESE

#' @title Weightfunction Prior
#'
#' @description Create weightfunction publication-bias priors and their
#' weight-prior helper objects.
#'
#' @param side character. Either \code{"one-sided"} or \code{"two-sided"}.
#' @param steps numeric vector of p-value cut points. These define
#' \code{length(steps) + 1} p-value bins and must be ordered values in
#' \code{(0, 1)}.
#' @param weights a weight-prior object created by \code{wf_cumulative()},
#' \code{wf_fixed()}, or \code{wf_independent()}.
#' @param reference character. Reference bin, currently
#' \code{"most_significant"}.
#' @param model fixed likelihood choices from \code{selection_model()}. The
#' default integrates estimate random effects and the complete sampling error,
#' conditions on other random effects, and uses product weights. These choices
#' do not change the weight-height prior or \code{prior_weights}.
#' @inheritParams prior
#'
#' @details Fixed weights must have one value per p-value bin
#' (\code{length(steps) + 1}), and the reference bin must have weight 1.
#'
#' @return \code{prior_weightfunction()} returns an object inheriting from
#' \code{prior} and \code{prior.weightfunction}; the \code{wf_*()} helpers
#' return \code{weightfunction_weights} helper objects with subclass markers.
#'
#' @examples
#' prior_weightfunction("one-sided", steps = 0.025)
#' prior_weightfunction(
#'   side    = "one-sided",
#'   steps   = c(0.025, 0.5),
#'   weights = wf_fixed(c(1, 0.8, 0.6))
#' )
#'
#' @seealso \code{\link{publication_bias_prior_specification}}
#'
#' @export
prior_weightfunction <- BayesTools::prior_weightfunction

#' @title Selection Model
#'
#' @description Declare the contextual sources retained by a selection model,
#' how publication weights combine, and an optional publication-group column.
#'
#' @param estimate_random_effects Either \code{"condition"} or
#' \code{"integrate"} (default), for the declared estimate-level true-effect
#' source. One-to-one grouping determines this role in multivariate models;
#' it does not imply independent covariance.
#' @param other_random_effects Either \code{"condition"} (default) or
#' \code{"integrate"}, for all remaining declared random-effect terms.
#' @param known_sampling_variance Either \code{"condition"} or
#' \code{"integrate"} (default), for the complete sampling-error vector
#' \eqn{e \sim N(0,V)}. Conditioning retains its realized value during
#' selection normalization; integration uses the full sampling covariance.
#' This choice also applies to univariate and independent sampling errors
#' supplied through \code{vi} or \code{sei}. It does not change \code{V}.
#' @param weight_rule Either \code{"product"} or \code{"best"}, specifying how
#' estimate weights combine within a publication group.
#' @param group Optional unresolved data-column name, supplied as a bare name,
#' a backtick-quoted name, or one character string. \code{NULL} requests
#' automatic group resolution when data are bound.
#'
#' @details This is RoBMA's re-export of \code{BayesTools::selection_model()}.
#' The specification is stored in \code{prior_weightfunction(model = ...)}.
#' Constructor \code{selection = selection_model(...)} settings apply only
#' to generated default priors; explicit priors retain their own specification.
#' A group name is captured without evaluating data; RoBMA binds it when the
#' model data are available. Conditioned sources remain unknown and estimated.
#' Integrating estimate variation is the default; conditioning on it defines
#' a different selection model. The specialized clustered interface retains
#' its total heterogeneity, allocation and within-/between-cluster I2 summaries.
#'
#' With positive sampling standard errors, the sampling law can be written as
#' \eqn{e = Sz}, with
#' \eqn{z \sim N(0,R)} and \eqn{V = SRS}, where \eqn{S} contains sampling
#' standard errors and \eqn{R} is their correlation matrix. Conditioning
#' retains all of \eqn{e}; no residual sampling error is integrated separately.
#' Covariance factorizations affect computation, not the selection model.
#' Selection thresholds continue to use the original sampling standard errors.
#'
#' When all outcome-generating sources are conditioned upon and weights are
#' positive almost surely, selection weights cancel and the observed law is
#' the ordinary Gaussian model. The likelihood then contains no information
#' about the selection weights. A retained context with zero acceptance
#' probability cannot complete the selection process. Integrating all sources
#' instead gives population-level selection among completed studies.
#' These choices are independent of prediction's \code{conditioning_depth}.
#' Non-unit observation \code{weights} require
#' \code{known_sampling_variance = "integrate"}; unit weights are equivalent
#' to omitting \code{weights}. Other observation-weight restrictions still apply.
#'
#' \code{weight_rule = "best"} uses the weight of the smallest p-value's bin,
#' rather than the largest weight. \code{weight_rule = "product"} multiplies
#' the estimate weights within each publication group.
#'
#' See [bselmodel.mv()] for the supported conditioning cells, source and group
#' requirements, and covariance boundaries.
#'
#' @return A \code{selection_model} object containing fixed mode choices and a
#' deferred column name or \code{NULL} for automatic group resolution.
#'
#' @examples
#' prior_weightfunction(
#'   steps = .05,
#'   model = selection_model(known_sampling_variance = "condition", group = paper_id)
#' )
#'
#' @seealso \code{\link{prior_weightfunction}}, [bselmodel.mv()], [RoBMA.mv()]
#' @md
#' @export
selection_model <- BayesTools::selection_model

#' @rdname prior_weightfunction
#' @param alpha optional positive cumulative-Dirichlet concentration parameters,
#' one per p-value bin. If \code{NULL}, \code{prior_weightfunction()} uses
#' \code{rep(1, length(steps) + 1)}. Cumulative weights encode monotone
#' decreasing publication weights relative to the most-significant bin.
#' @export
wf_cumulative <- BayesTools::wf_cumulative

#' @rdname prior_weightfunction
#' @param omega fixed publication weights, one per bin; values must be
#' non-missing, nonnegative, and match \code{length(steps) + 1} when used in
#' \code{prior_weightfunction()}.
#' @export
wf_fixed <- BayesTools::wf_fixed

#' @rdname prior_weightfunction
#' @param prior continuous simple prior distribution for each non-reference
#' weight. Point, discrete, mixture, and other non-simple priors are invalid.
#' @param scale latent scale for independent weights; either \code{"omega"},
#' \code{"log_omega"}, or the \code{"log"} alias. Direct \code{"omega"} priors
#' need nonnegative support; \code{"log"} is normalized to \code{"log_omega"}.
#' @export
wf_independent <- BayesTools::wf_independent

#' @title Informed Prior
#'
#' @description Create empirical informed prior distributions.
#'
#' @param name character. Informed prior name. Supported examples include
#' \code{"van Erp"}, \code{"Oosterwijk"}, \code{"Cochrane"}, and medicine
#' subfields from \code{BayesTools::prior_informed_medicine_names}.
#' @param parameter character. Parameter subset. Required for medicine priors
#' and not needed for psychology priors.
#' @param type character. Effect-size type. Medicine priors use values such as
#' \code{"smd"}, \code{"logOR"}, \code{"logRR"}, \code{"RD"}, and
#' \code{"logHR"}.
#'
#' @return An object of class \code{prior}.
#'
#' @details Further details can be found in \insertCite{erp2017estimates;textual}{RoBMA},
#' \insertCite{gronau2017bayesian;textual}{RoBMA},
#' \insertCite{bartos2021bayesian;textual}{RoBMA}, and
#' \insertCite{bartos2023empirical;textual}{RoBMA}.
#' @references
#' \insertAllCited{}
#' @export
prior_informed <- BayesTools::prior_informed

# needs to be completely re-exported because of some odd roxygen issues
#' @title BayesTools Contrast Matrices
#'
#' @description BayesTools provides several contrast matrix functions for Bayesian factor analysis.
#' These functions create different types of contrast matrices that can be used with factor
#' variables in Bayesian models.
#'
#' @details
#' The package includes the following contrast functions:
#' \describe{
#'   \item{\code{contr.orthonormal}}{Return a matrix of orthonormal contrasts.
#'     Code is based on \code{stanova::contr.bayes} and corresponding to description
#'     by \insertCite{rouder2012default;textual}{BayesTools}. Returns a matrix with n rows and
#'     k columns, with k = n - 1 if \code{contrasts = TRUE} and k = n if \code{contrasts = FALSE}.}
#'   \item{\code{contr.meandif}}{Return a matrix of mean difference contrasts.
#'     This is an adjustment to the \code{contr.orthonormal} that ascertains that the prior
#'     distributions on difference between the grand mean and factor level are identical independent
#'     of the number of factor levels (which does not hold for the orthonormal contrast). Furthermore,
#'     the contrast is re-scaled so the specified prior distribution exactly corresponds to the prior
#'     distribution on difference between each factor level and the grand mean -- this is approximately
#'     twice the scale of \code{contr.orthonormal}. Returns a matrix with n rows and k columns,
#'     with k = n - 1 if \code{contrasts = TRUE} and k = n if \code{contrasts = FALSE}.}
#'   \item{\code{contr.independent}}{Return a matrix of independent contrasts -- a level for each term.
#'     Returns a matrix with n rows and k columns, with k = n if \code{contrasts = TRUE} and k = n
#'     if \code{contrasts = FALSE}.}
#' }
#'
#' @param n a vector of levels for a factor, or the number of levels
#' @param contrasts logical indicating whether contrasts should be computed
#'
#' @examples
#' # Orthonormal contrasts
#' contr.orthonormal(c(1, 2))
#' contr.orthonormal(c(1, 2, 3))
#'
#' # Mean difference contrasts
#' contr.meandif(c(1, 2))
#' contr.meandif(c(1, 2, 3))
#'
#' # Independent contrasts
#' contr.independent(c(1, 2))
#' contr.independent(c(1, 2, 3))
#'
#' @references
#' \insertAllCited{}
#'
#' @aliases contr.orthonormal contr.meandif contr.independent
#' @name contr.BayesTools
NULL

#' @rdname contr.BayesTools
#' @export
contr.orthonormal <- BayesTools::contr.orthonormal

#' @rdname contr.BayesTools
#' @export
contr.meandif <- BayesTools::contr.meandif

#' @rdname contr.BayesTools
#' @export
contr.independent <- BayesTools::contr.independent
