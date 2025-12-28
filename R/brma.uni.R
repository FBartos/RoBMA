### RoBMA 4.0.0


#' @title Bayesian Random-Effects Meta-Analysis
#'
#' @description
#'
#'
#' @details
#' ## Prior distributions
#' Prior distributions must be specified for all model parameters.
#' This typically includes the pooled effect \eqn{\mu} and between-study heterogeneity
#' \eqn{\tau}. In the case of meta-regression, the pooled effect \eqn{\mu} corresponds to
#' the intercept, and additional prior distributions for the regression coefficients
#' are required. In the case of a location-scale model, the between-study heterogeneity
#' corresponds to the intercept of the scale regression, and additional prior
#' distributions for the scale regression coefficients are required.
#'
#' There are several ways to specify the prior distributions: \enumerate{
#'    \item{via a standardized effect size `measure` with known unit information standard deviation,}
#'    \item{by estimating unit information standard deviation using sample sizes `ni`,}
#'    \item{by manually setting `prior_unit_information_sd`,}
#'    \item{by specifying informed empirical prior distributions via `prior_informed_field`
#'    and `prior_informed_subfield`,}
#'    \item{or via fully custom specification using the `prior_effect`, `prior_heterogeneity`,
#'    `prior_mods`, `prior_scale`, and `prior_heterogeneity_allocation` arguments.}
#' }
#' In all cases, the prior behavior can be further modified by the `rescale_priors`,
#' `standardize_continuous_predictors`, and `set_contrast_factor_predictors` arguments.
#'
#' ### (1) Specifying prior distributions via standardized effect size measures with known
#' unit information standard deviation
#' This is the easiest way to specify prior distributions. The width of prior
#' distributions is based on a fraction of the known unit information standard deviation (`UISD`)
#' \insertCite{rover2021weakly}{RoBMA}. The default prior distributions for the parameters
#' are set as follows:\tabular{ll}{
#'    effect size:               \tab Normal(0, \eqn{\frac{1}{2}} `UISD`)  \cr
#'    heterogeneity:             \tab Normal+(0, \eqn{\frac{1}{4}} `UISD`) \cr
#'    effect moderation:         \tab Normal(0, \eqn{\frac{1}{4}} `UISD`)  \cr
#'    heterogeneity moderation:  \tab Normal(0, \eqn{\frac{1}{4}})
#' }
#' The heterogeneity moderation parameters are multiplicative, as such they are independent of `UISD`.
#'
#' The default fraction of the `UISD` can be changed via `RoBMA.options()` using one of
#' the following arguments: `"default_UISD.effect"`, `"default_UISD.heterogeneity"`,
#' `"default_UISD.mods"`, `"default_UISD.scale"`.
#'
#' The known `UISD` for standardized effect size measures (`measure`) are set as follows:
#' \tabular{ll}{
#'   `"SMD"`:   \tab \eqn{\sqrt{2}} \cr
#'   `"ZCOR"`:  \tab \eqn{1}        \cr
#'   `"RR"`:    \tab \eqn{\sqrt{4}} \cr
#'   `"OR"`:    \tab \eqn{\sqrt{4}} \cr
#'   `"HR"`:    \tab \eqn{\sqrt{4}}
#' }
#' See Chapter 2.4 in \insertCite{spiegelhalter2004bayesian;textual}{RoBMA} and
#' Chapter 1 in \insertCite{grieve2022hybrid;textual}{RoBMA}.
#'
#' ### (2) Estimating unit information standard deviation using sample sizes
#' When effect sizes are on a non-standardized scale (`measure = "GEN"`) or use a
#' standardized effect size without known `UISD`, the `UISD` can be estimated from
#' sample sizes (`ni`) and standard errors (`sei`) following Equation 6 in
#' \insertCite{rover2021weakly;textual}{RoBMA}. The estimated `UISD` is then used to
#' scale the default prior distributions as described in section (1).
#'
#' Note that the known `UISD` for standardized effect size measures (section (1)) is
#' used if available, even when `ni` is provided.
#'
#' ### (3) Manually setting unit information standard deviation
#' Alternatively, the `UISD` can be specified directly via the `prior_unit_information_sd`
#' argument. This is useful when the appropriate scale for prior distributions is known
#' a priori or when multiple analyses are to be performed on subsets of the same data
#' (re-estimating `UISD` on different subsets of the data can lead to slightly different
#' prior distributions for different subsets; see [estimate_unit_information_sd()]).
#' The specified `prior_unit_information_sd` is then used to scale the default prior
#' distributions as described in section (1).
#'
#' Note that the manually specified `prior_unit_information_sd` takes precedence over
#' the estimated `UISD` from `ni` (section (2)) and the known `UISD` from `measure`
#' (section (1)).
#'
#' ### (4) Specifying informed empirical prior distributions
#' Informed prior distributions can be specified via the `prior_informed_field` and
#' `prior_informed_subfield` arguments. Currently, only `prior_informed_field = "medicine"`
#' with subfields defined in `BayesTools::prior_informed_medicine_names` is supported,
#' which uses empirically derived prior distributions from medical meta-analyses as described in
#' \insertCite{bartos2021bayesian;textual}{RoBMA} and
#' \insertCite{bartos2023empirical;textual}{RoBMA}.
#'
#' When `prior_informed_field = "medicine"`, the default `prior_informed_subfield` is
#' `"Cochrane"` (i.e., using the whole CDSR database as a reference). The informed
#' prior distributions are available for the following effect size measures:
#' \tabular{ll}{
#'   `"SMD"`:  \tab standardized mean difference \cr
#'   `"OR"`:   \tab log odds ratio               \cr
#'   `"RR"`:   \tab log risk ratio               \cr
#'   `"RD"`:   \tab risk difference              \cr
#'   `"HR"`:   \tab log hazard ratio
#' }
#'
#' Note that informed prior distributions are only available for the effect size (\eqn{\mu})
#' and heterogeneity (\eqn{\tau}) parameters. For effect moderation (`prior_mods`), the
#' informed effect prior is scaled by a factor of `RoBMA.get_option("default_informed_priors.mods")`.
#' For heterogeneity moderation (`prior_scale`), a normal prior with standard deviation
#' specified by `RoBMA.get_option("default_informed_priors.scale")` is used.
#'
#' ### (5) Fully custom prior distributions
#' Prior distributions can be fully customized by directly specifying the `prior_effect`,
#' `prior_heterogeneity`, `prior_mods`, `prior_scale`, and `prior_heterogeneity_allocation`
#' arguments. These should be prior distribution objects created via `BayesTools::prior()`
#' or related functions (e.g., `prior_factor()`).
#'
#' ### Rescaling prior distributions
#' The `rescale_priors` argument allows rescaling of all prior distributions by a
#' multiplicative factor. For example, `rescale_priors = 2` doubles the standard
#' deviations/scales of all prior distributions, making them more diffuse. This applies
#' regardless of how the priors were specified.
#'
#' ### Handling of continuous and factor predictors
#' TODO: explain `standardize_continuous_predictors` and `set_contrast_factor_predictors`
#'
#' @references
#' \insertAllCited{}
#'
#' @export
brma.uni <- function(
    # input specification
    yi, vi, sei, ni, weights,
    mods, scale, study_ids,
    data, slab, subset,
    measure = "GEN",

    # prior specification
    prior_effect, prior_heterogeneity, prior_mods, prior_scale,
    prior_heterogeneity_allocation,
    standardize_continuous_predictors = TRUE,
    set_contrast_factor_predictors = "treatment",
    prior_unit_information_sd, rescale_priors = 1,
    prior_informed_field, prior_informed_subfield,

    # MCMC fitting settings
    sample = 5000, burnin = 2000, adapt = 500,
    chains = 3, thin = 1, parallel = FALSE,
    autofit = FALSE, autofit_control = set_autofit_control(),
    convergence_checks = set_convergence_checks(),

    # additional settings
    seed = NULL, silent = TRUE, ...
    ){

  ### create the output object
  time.start   <- proc.time()
  dots         <- list(...)
  object       <- .createObject(
    dots = dots, class = "brma.uni",
    # MCMC and fitting settings
    chains = chains, adapt = adapt, burnin = burnin, sample = sample, thin = thin,
    autofit = autofit, parallel = parallel, silent = silent, seed = seed,
    autofit_control = autofit_control, convergence_checks = convergence_checks,
    # additional options
    standardize_continuous_predictors = standardize_continuous_predictors
  )

  ### check and store the data
  object$data <- .check_and_list_data(.call = match.call(), .envir = parent.frame())
  if (isTRUE(dots[["only_data"]]))
    return(object)

  ### check and store priors
  object$priors <- .check_and_list_priors.brma(
    prior_effect = prior_effect, prior_heterogeneity = prior_heterogeneity,
    prior_mods = prior_mods, prior_scale = prior_scale,
    prior_heterogeneity_allocation = prior_heterogeneity_allocation,
    set_contrast_factor_predictors    = set_contrast_factor_predictors,
    rescale_priors                    = rescale_priors,
    prior_unit_information_sd         = prior_unit_information_sd,
    prior_informed_field              = prior_informed_field,
    prior_informed_subfield           = prior_informed_subfield,
    data = object[["data"]], measure = measure)
  if (isTRUE(dots[["only_priors"]]))
    return(object)

  ### fit the model
  object$fit <- .fit(object)

  ### store simple summary & coefficients
  object$summary       <- summary(object)
  object$coefficients  <- coefficients(object)

  return(object)
}

#' @export
summary.brma.uni <- function(object, ...) {

  estimates <- BayesTools::JAGS_estimates_table(
    fit               = object[["fit"]],
    transform_factors = TRUE
  )

  return(estimates)
}
#' @export
coefficients.brma.uni <- function(object, ...) {

  estimates            <- object[["summary"]][,"Mean"]
  colnames(estimates)  <- rownames(object[["summary"]])

  return(estimates)
}
#' @export
print.brma.uni <- function(x, ...) {
  return(x[["summary"]])
}
