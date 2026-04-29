#' @title Prior Distribution
#'
#' @description Thin re-export of \code{\link[BayesTools]{prior}} for creating
#' prior distribution objects.
#'
#' @param distribution character. Prior distribution name.
#' @param parameters list. Distribution parameters.
#' @param truncation list with \code{lower} and \code{upper} truncation bounds.
#' @param prior_weights numeric prior model weight.
#'
#' @return An object of class \code{prior}.
#'
#' @export
prior <- BayesTools::prior

#' @title Empty Prior
#'
#' @description Thin re-export of \code{\link[BayesTools]{prior_none}} for
#' omitting a model component.
#'
#' @inheritParams prior
#'
#' @return An object of class \code{prior}.
#'
#' @export
prior_none <- BayesTools::prior_none

#' @title Factor Prior
#'
#' @description Thin re-export of \code{\link[BayesTools]{prior_factor}} for
#' creating priors on factor contrasts.
#'
#' @inheritParams prior
#' @param contrast character. Contrast coding used for factor levels.
#'
#' @return An object of class \code{prior}.
#'
#' @export
prior_factor <- BayesTools::prior_factor

#' @title PET Prior
#'
#' @description Thin re-export of \code{\link[BayesTools]{prior_PET}} for PET
#' publication-bias regression priors.
#'
#' @inheritParams prior
#'
#' @return An object of class \code{prior}.
#'
#' @export
prior_PET  <- BayesTools::prior_PET

#' @title PEESE Prior
#'
#' @description Thin re-export of \code{\link[BayesTools]{prior_PEESE}} for
#' PEESE publication-bias regression priors.
#'
#' @inheritParams prior
#'
#' @return An object of class \code{prior}.
#'
#' @export
prior_PEESE <- BayesTools::prior_PEESE

#' @title Weightfunction Prior
#'
#' @description Thin re-export of \code{\link[BayesTools]{prior_weightfunction}}
#' and its weight-prior helper constructors.
#'
#' @param side character. Either \code{"one-sided"} or \code{"two-sided"}.
#' @param steps numeric vector of p-value cut points.
#' @param weights a weight-prior object created by \code{wf_cumulative()},
#' \code{wf_fixed()}, or \code{wf_independent()}.
#' @param reference character. Reference bin, currently
#' \code{"most_significant"}.
#' @inheritParams prior
#'
#' @return \code{prior_weightfunction()} returns an object of class
#' \code{prior}; the \code{wf_*()} helpers return weight-prior objects.
#'
#' @export
prior_weightfunction <- BayesTools::prior_weightfunction

#' @rdname prior_weightfunction
#' @param alpha positive cumulative-Dirichlet concentration parameters.
#' @export
wf_cumulative <- BayesTools::wf_cumulative

#' @rdname prior_weightfunction
#' @param omega fixed publication weights, one per bin.
#' @export
wf_fixed <- BayesTools::wf_fixed

#' @rdname prior_weightfunction
#' @param prior prior distribution for each non-reference weight.
#' @param scale latent scale for independent weights; either \code{"omega"},
#' \code{"log_omega"}, or the \code{"log"} alias.
#' @export
wf_independent <- BayesTools::wf_independent

#' @title Informed Prior
#'
#' @description Thin re-export of \code{\link[BayesTools]{prior_informed}} for
#' empirical informed prior distributions.
#'
#' @param name character. Informed prior name.
#' @param parameter character. Optional parameter subset.
#' @param type character. Effect-size type.
#'
#' @return An object of class \code{prior}.
#'
#' @details Further details can be found in \insertCite{erp2017estimates;textual}{RoBMA},
#' \insertCite{gronau2017bayesian;textual}{RoBMA}, and
#' \insertCite{bartos2021bayesian;textual}{RoBMA}.
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
#'     distributions on difference between the gran mean and factor level are identical independent
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
