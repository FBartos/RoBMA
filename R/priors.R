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

#' @title Orthornomal contrast matrix
#'
#' @description Return a matrix of orthornomal contrasts.
#' Code is based on \code{stanova::contr.bayes} and corresponding to description
#' by \insertCite{rouder2012default;textual}{BayesTools}
#'
#' @param n a vector of levels for a factor, or the number of levels
#' @param contrasts logical indicating whether contrasts should be computed
#'
#' @examples
#' contr.orthonormal(c(1, 2))
#' contr.orthonormal(c(1, 2, 3))
#'
#' @references
#' \insertAllCited{}
#'
#' @return A matrix with n rows and k columns, with k = n - 1 if \code{contrasts = TRUE} and k = n
#' if \code{contrasts = FALSE}.
#'
#' @export
contr.orthonormal <- BayesTools::contr.orthonormal

#' @title Mean difference contrast matrix
#'
#' @description Return a matrix of mean difference contrasts.
#' This is an adjustment to the \code{contr.orthonormal} that ascertains that the prior
#' distributions on difference between the gran mean and factor level are identical independent
#' of the number of factor levels (which does not hold for the orthonormal contrast). Furthermore,
#' the contrast is re-scaled so the specified prior distribution exactly corresponds to the prior
#' distribution on difference between each factor level and the grand mean -- this is approximately
#' twice the scale of \code{contr.orthonormal}.
#'
#' @param n a vector of levels for a factor, or the number of levels
#' @param contrasts logical indicating whether contrasts should be computed
#'
#' @examples
#' contr.meandif(c(1, 2))
#' contr.meandif(c(1, 2, 3))
#'
#' @references
#' \insertAllCited{}
#'
#' @return A matrix with n rows and k columns, with k = n - 1 if \code{contrasts = TRUE} and k = n
#' if \code{contrasts = FALSE}.
#'
#' @export
contr.meandif <- BayesTools::contr.meandif

#' @title Independent contrast matrix
#'
#' @description Return a matrix of independent contrasts -- a level for each term.
#'
#' @param n a vector of levels for a factor, or the number of levels
#' @param contrasts logical indicating whether contrasts should be computed
#'
#' @examples
#' contr.independent(c(1, 2))
#' contr.independent(c(1, 2, 3))
#'
#' @references
#' \insertAllCited{}
#'
#' @return A matrix with n rows and k columns, with k = n if \code{contrasts = TRUE} and k = n
#' if \code{contrasts = FALSE}.
#'
#' @export
contr.independent <- function(n, contrasts = TRUE){

  if(length(n) <= 1L){
    if(is.numeric(n) && length(n) == 1L && n >= 1L){
      return(TRUE)
    }else{
      stop("Not enough degrees of freedom to define contrasts.")
    }
  }else{
    n <- length(n)
  }

  cont <- diag(x = 1, nrow = n, ncol = n)

  return(cont)
}


#' @title Set default prior distributions
#'
#' @description Set default prior distributions for RoBMA models.
#'
#' @param parameter a character string specifying the parameter for
#' which the prior distribution should be set. Available options are
#' "effect", "heterogeneity", "bias", "hierarchical", "covariates",
#' "factors".
#' @param null a logical indicating whether the prior distribution
#' should be set for the null hypothesis. Defaults to \code{FALSE}.
#' @param rescale a numeric value specifying the re-scaling factor
#' for the default prior distributions. Defaults to 1. Allows
#' convenient re-scaling of prior distributions simultaneously.
#'
#' @details The default prior distributions corresponds to the
#' specification of RoBMA-PSMA and RoBMA-regression outlined in
#' \insertCite{bartos2021no;textual}{RoBMA} and
#' \insertCite{bartos2023robust;textual}{RoBMA}.
#'
#' Specifically, the prior distributions are:
#'
#' **For the alternative hypothesis:**
#' \itemize{
#'   \item \strong{Effect:} Normal distribution with mean 0 and standard deviation 1.
#'   \item \strong{Heterogeneity:} Inverse gamma distribution with shape 1 and scale 0.15.
#'   \item \strong{Bias:} A list of 8 prior distributions defining the publication bias adjustments:
#'   \itemize{
#'     \item Two-sided: Weight function with steps 0.05.
#'     \item Two-sided: Weight function with steps 0.05 and 0.1.
#'     \item One-sided: Weight function with steps 0.05.
#'     \item One-sided: Weight function with steps 0.025 and 0.05.
#'     \item One-sided: Weight function with steps 0.05 and 0.5.
#'     \item One-sided: Weight function with steps 0.025, 0.05, and 0.5.
#'     \item PET-type model with regression coefficient: Cauchy distribution with location 0 and scale 1.
#'     \item PEESE-type model with regression coefficient: Cauchy distribution with location 0 and scale 5.
#'   }
#'   All weight functions use a unit cumulative Dirichlet prior distribution on relative prior probabilities.
#'   \item \strong{Standardized continuous covariates:} Normal distribution with mean 0 and standard deviation 0.25.
#'   \item \strong{Factors (via by-level differences from the grand mean):} Normal distribution with mean 0 and standard deviation 0.25.
#' }
#'
#' **For the null hypothesis:**
#' \itemize{
#'   \item \strong{Effect:} Point distribution at 0.
#'   \item \strong{Heterogeneity:} Point distribution at 0.
#'   \item \strong{Bias:} No prior distribution.
#'   \item \strong{Standardized continuous covariates:} Point distribution at 0.
#'   \item \strong{Factors (via by-level differences from the grand mean):} Point distribution at 0.
#' }
#'
#' The rescaling factor adjusts the width of the effect, heterogeneity, covariates, factor, and PEESE-style model prior distributions.
#' PET-style and weight function prior distributions are scale-invariant.
#'
#' @examples
#'
#' set_default_priors("effect")
#' set_default_priors("heterogeneity")
#' set_default_priors("bias")
#'
#' @return A prior distribution object or a list of prior distribution
#' objects.
#'
#' @export
set_default_priors <- function(parameter, null = FALSE, rescale = 1){

  BayesTools::check_char(parameter, "parameter", allow_values = c("effect", "heterogeneity", "bias", "hierarchical", "covariates", "factors"))
  BayesTools::check_bool(null, "null")
  BayesTools::check_real(rescale, "rescale", lower = 0)

  if(null){
    return(switch(
      parameter,
      effect        = prior(distribution = "point", parameters = list(location = 0)),
      heterogeneity = prior(distribution = "point", parameters = list(location = 0)),
      bias          = prior_none(),
      hierarchical  = NULL,
      covariates    = prior(distribution = "point", parameters = list(location = 0)),
      factors       = prior_factor("spike", parameters = list(location = 0), contrast = "meandif")
    ))
  }else{
    return(switch(
      parameter,
      effect        = prior(distribution = "normal",   parameters = list(mean  = 0, sd = 1 * rescale)),
      heterogeneity = prior(distribution = "invgamma", parameters = list(shape = 1, scale = 0.15 * rescale)),
      bias          = list(
        prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1),       steps = c(0.05)),             prior_weights = 1/12),
        prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.05, 0.1)),        prior_weights = 1/12),
        prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1),       steps = c(0.05)),             prior_weights = 1/12),
        prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.025, 0.05)),      prior_weights = 1/12),
        prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.05, 0.5)),        prior_weights = 1/12),
        prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1, 1), steps = c(0.025, 0.05, 0.5)), prior_weights = 1/12),
        prior_PET(distribution = "Cauchy",   parameters = list(0, 1),            truncation = list(0, Inf), prior_weights = 1/4),
        prior_PEESE(distribution = "Cauchy", parameters = list(0, 5 / rescale),  truncation = list(0, Inf), prior_weights = 1/4)
        ),
      hierarchical  = prior("beta", parameters = list(alpha = 1, beta = 1)),
      covariates    = prior("normal", parameters = list(mean = 0, sd = 0.25 * rescale)),
      factors       = prior_factor("mnormal", parameters = list(mean = 0, sd = 0.25 * rescale), contrast = "meandif")
    ))
  }
}

#' @title Set default prior distributions for binomial meta-analytic models
#'
#' @description Set default prior distributions for BiBMA models.
#'
#' @param parameter a character string specifying the parameter for
#' which the prior distribution should be set. Available options are
#' "effect", "heterogeneity", "baseline", "covariates",
#' "factors".
#' @param null a logical indicating whether the prior distribution
#' should be set for the null hypothesis. Defaults to \code{FALSE}.
#' @param rescale a numeric value specifying the re-scaling factor
#' for the default prior distributions. Defaults to 1. Allows
#' convenient re-scaling of prior distributions simultaneously.
#'
#' @details The default prior are based on the binary outcome meta-analyses
#' in the Cochrane Database of Systematic Reviews outlined in
#'  \insertCite{bartos2023empirical;textual}{RoBMA}.
#'
#' Specifically, the prior distributions are:
#'
#' **For the alternative hypothesis:**
#' \itemize{
#'   \item \strong{Effect:} T distribution with mean 0, scale 0.58, and 4 degrees of freedom.
#'   \item \strong{Heterogeneity:} Inverse gamma distribution with shape 1.77 and scale 0.55.
#'   \item \strong{Baseline:} No prior distribution.
#'   \item \strong{Standardized continuous covariates:} Normal distribution with mean 0 and standard deviation 0.29.
#'   \item \strong{Factors (via by-level differences from the grand mean):} Normal distribution with mean 0 and standard deviation 0.29.
#' }
#'
#' **For the null hypothesis:**
#' \itemize{
#'   \item \strong{Effect:} Point distribution at 0.
#'   \item \strong{Heterogeneity:} Point distribution at 0.
#'   \item \strong{Baseline:} Independent uniform distributions.
#'   \item \strong{Standardized continuous covariates:} Point distribution at 0.
#'   \item \strong{Factors (via by-level differences from the grand mean):} Point distribution at 0.
#' }
#'
#' The rescaling factor adjusts the width of the effect, heterogeneity, covariates, factor, and PEESE-style model prior distributions.
#' PET-style and weight function prior distributions are scale-invariant.
#'
#' @examples
#'
#' set_default_binomial_priors("effect")
#' set_default_binomial_priors("heterogeneity")
#' set_default_binomial_priors("baseline")
#'
#' @return A prior distribution object or a list of prior distribution
#' objects.
#'
#' @export
set_default_binomial_priors <- function(parameter, null = FALSE, rescale = 1){

  BayesTools::check_char(parameter, "parameter", allow_values = c("effect", "heterogeneity", "baseline", "covariates", "factors", "hierarchical"))
  BayesTools::check_bool(null, "null")
  BayesTools::check_real(rescale, "rescale", lower = 0)


  if(null){
    return(switch(
      parameter,
      effect        = prior(distribution = "point", parameters = list(location = 0)),
      heterogeneity = prior(distribution = "point", parameters = list(location = 0)),
      baseline      = prior_factor("beta", parameters = list(alpha = 1, beta = 1), contrast = "independent"),
      hierarchical  = NULL,

      covariates    = prior(distribution = "point", parameters = list(location = 0)),
      factors       = prior_factor("spike", parameters = list(location = 0), contrast = "meandif")
    ))
  }else{
    return(switch(
      parameter,
      effect        = prior(distribution = "student",   parameters = list(location = 0, scale = 0.58 * rescale, df = 4)),
      heterogeneity = prior(distribution = "invgamma",  parameters = list(shape = 1.77, scale = 0.55 * rescale)),
      baseline      = NULL,
      hierarchical  = prior("beta", parameters = list(alpha = 1, beta = 1)),

      covariates    = prior("normal", parameters = list(mean = 0, sd = 0.58 * (1/2) * rescale)),
      factors       = prior_factor("mnormal", parameters = list(mean = 0, sd = 0.58 * (1/2) * rescale), contrast = "meandif")
    ))
  }
}
