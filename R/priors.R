#' @name prior
#' @inherit BayesTools::prior
#' @export
prior <- BayesTools::prior

#' @name prior_none
#' @inherit BayesTools::prior_none
#' @export
prior_none <- BayesTools::prior_none

#' @name prior_factor
#' @inherit BayesTools::prior_factor
#' @export
prior_factor <- BayesTools::prior_factor

#' @name prior_PET
#' @inherit BayesTools::prior_PET
#' @export
prior_PET  <- BayesTools::prior_PET

#' @name prior_PEESE
#' @inherit BayesTools::prior_PEESE
#' @export
prior_PEESE <- BayesTools::prior_PEESE

#' @name prior_weightfunction
#' @inherit BayesTools::prior_weightfunction
#' @export
prior_weightfunction <- BayesTools::prior_weightfunction

#' @name prior_informed
#' @inherit BayesTools::prior_informed
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
