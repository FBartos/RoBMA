#' @title Prior specification
#' @name prior_specification
#'
#' @description
#' Prior distributions must be specified for all model parameters.
#' This typically includes the pooled effect \eqn{\mu} and between-study heterogeneity
#' \eqn{\tau}. In the case of meta-regression, the pooled effect \eqn{\mu} corresponds to
#' the intercept, and additional prior distributions for the regression coefficients
#' are required. In the case of a location-scale model, the between-study heterogeneity
#' corresponds to the intercept of the scale regression, and additional prior
#' distributions for the scale regression coefficients are required.
#'
#' @param prior_effect prior distribution for the effect size (\eqn{\mu}) parameter
#' (the intercept). Defaults to `NULL`.
#' @param prior_heterogeneity prior distribution for the heterogeneity (\eqn{\tau})
#' parameter. Defaults to `NULL`.
#' @param prior_mods prior distribution for the moderators (\eqn{\beta}) parameters.
#' Defaults to `NULL`.
#' @param prior_scale prior distribution for the scale (\eqn{\delta}) parameters.
#' Defaults to `NULL`.
#' @param prior_heterogeneity_allocation prior distribution for the fraction of variance
#' explained by the moderators across studies (\eqn{\nu}). Defaults to `NULL`.
#' @param prior_bias prior distribution for the publication bias adjustment parameters
#' (selection models) or the publication bias indicator (PET-PEESE). Defaults to `NULL`.
#' @param prior_unit_information_sd numeric. The unit information standard deviation (\eqn{\sigma_{unit}}).
#' Defaults to `NULL`.
#' @param rescale_priors numeric. A scaling factor for the prior distributions. Defaults to 1.
#' @param standardize_continuous_predictors logical. Whether to standardize continuous predictors.
#' Defaults to `TRUE`.
#' @param set_contrast_factor_predictors character. How to set contrast for factor predictors.
#' Defaults to `"treatment"`.
#' @param prior_informed_field character. The field of the informed prior distributions.
#' Defaults to `NULL`.
#' @param prior_informed_subfield character. The subfield of the informed prior distributions.
#' Defaults to `NULL`.
#'
#' @details
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
#'    heterogeneity moderation:  \tab Normal(0, \eqn{\frac{1}{2}})
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
#'   `"HR"`:    \tab \eqn{\sqrt{4}} \cr
#'   `"IRR"`:   \tab \eqn{\sqrt{4}}
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
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link[BayesTools]{prior}}, \code{\link{RoBMA.options}}, \code{\link{brma}}
#' @aliases prior_specification
NULL


#' @title Prior specification for model-averaging
#' @name RoBMA_prior_specification
#'
#' @description
#' The \code{\link{RoBMA}} function extends the prior specification described in
#' \code{\link{prior_specification}} by allowing mixture priors that define competing
#' hypotheses for Bayesian model-averaging. This enables multimodel inference via the
#' product space method \insertCite{carlin1995bayesian,lodewyckx2011tutorial}{RoBMA},
#' where a single MCMC chain explores a joint model space containing all competing models.
#'
#' @param prior_effect,prior_heterogeneity,prior_mods,prior_scale,prior_heterogeneity_allocation,prior_bias
#' Prior distribution(s) for the alternative hypothesis component(s). Can be:
#' \itemize{
#'   \item A single prior distribution object (creates a mixture with one alternative)
#'   \item A list of prior distributions (creates a mixture with multiple alternatives)
#'   \item `NULL` or `FALSE` (omits the alternative hypothesis component)
#' }
#' See \code{\link{prior_specification}} for details on specifying individual priors.
#'
#' @param prior_effect_null,prior_heterogeneity_null,prior_mods_null,prior_scale_null,prior_heterogeneity_allocation_null,prior_bias_null
#' Prior distribution(s) for the null hypothesis component(s). Can be:
#' \itemize{
#'   \item A single prior distribution object (creates a mixture with one null)
#'   \item A list of prior distributions (creates a mixture with multiple nulls)
#'   \item `NULL` or `FALSE` (omits the null hypothesis component)
#' }
#' Defaults to a point mass (spike) at zero for effect and heterogeneity parameters.
#'
#' @param model_type Character string specifying predefined publication bias model
#' ensembles. One of:
#' \itemize{
#'   \item `"PSMA"` (default): Full RoBMA-PSMA ensemble with 6 weight functions + PET + PEESE
#'   \item `"6w"`: Six weight function models
#'   \item `"2w"`: Two weight function models 
#'   \item `"PP"`: PET-PEESE models only
#' }
#' Custom `prior_bias` and `prior_bias_null` override this setting.
#'
#' @details
#' ## The product space method for model-averaging
#'
#' RoBMA implements Bayesian model-averaging using the product space method
#' \insertCite{carlin1995bayesian,lodewyckx2011tutorial}{RoBMA}. Rather than fitting
#' each model separately and combining results after the fitting is complete, 
#' this approach embeds all competing models within a single joint model space. 
#' A discrete model indicator variable selects which model component is "active" 
#' at each MCMC iteration, and the posterior probability of each model is estimated 
#' from the proportion of iterations spent in that model's subspace.
#'
#' This is achieved by specifying **mixture priors** that combine:
#' \itemize{
#'   \item **Null hypothesis component(s)**: Typically point masses (spikes) representing
#'         the absence of an effect, heterogeneity, or publication bias
#'   \item **Alternative hypothesis component(s)**: Continuous distributions representing
#'         the presence of the parameter
#' }
#'
#' The mixture structure enables direct computation of:
#' \itemize{
#'   \item **Inclusion Bayes factors**: Evidence for the presence vs. absence of each
#'         parameter (e.g., effect, heterogeneity, bias)
#'   \item **Model-averaged estimates**: Posterior estimates that account for model
#'         uncertainty by weighting across all models
#'   \item **Conditional estimates**: Posterior estimates given that a particular
#'         hypothesis is true
#' }
#'
#' ## Default mixture prior structure
#'
#' By default, \code{\link{RoBMA}} creates the following mixture priors:
#'
#' \tabular{lll}{
#'   **Parameter**      \tab **Null component**         \tab **Alternative component**   \cr
#'   Effect (\eqn{\mu})        \tab Spike(0)                   \tab Normal(0, UISD-scaled)      \cr
#'   Heterogeneity (\eqn{\tau}) \tab Spike(0)                   \tab Normal+(0, UISD-scaled)     \cr
#'   Publication bias   \tab No bias (prior_none)       \tab Weight functions + PET + PEESE \cr
#'   Allocation (\eqn{\rho})   \tab (none by default)          \tab Beta(1, 1)
#' }
#'
#' See \code{\link{prior_specification}} for details on how the UISD scaling is determined.
#'
#' ## Specifying custom mixture priors
#'
#' ### Single prior vs. list of priors
#'
#' Each `prior_*` and `prior_*_null` argument accepts either a single prior or a list:
#' \itemize{
#'   \item **Single prior**: Creates one component in the mixture
#'   \item **List of priors**: Creates multiple components, enabling model-averaging
#'         across different prior specifications (e.g., different effect size scales)
#' }
#'
#' When multiple priors are specified in a list, they are all treated as alternative
#' (or null) hypothesis components for the main inclusion Bayes factor. The inclusion
#' Bayes factor tests the combined evidence for all alternative components against all
#' null components. The detailed summary output shows posterior probabilities and 
#' inclusion Bayes factors for each individual component separately, allowing finer-grained
#' inference about which specific prior specification is most supported by the data.
#'
#' ### Omitting hypothesis components
#'
#' Setting a prior argument to `NULL` or `FALSE` removes that hypothesis component:
#' \itemize{
#'   \item `prior_effect_null = NULL`: No null hypothesis for effect (assumes effect exists)
#'   \item `prior_effect = NULL`: No alternative hypothesis (assumes no effect)
#'   \item `prior_bias_null = NULL`: No "no bias" model (assumes some bias mechanism)
#'   \item `prior_bias = NULL`: No bias model (assumes no bias mechanism)
#' }
#'
#' ### Parameter-specific customization for moderators
#'
#' For moderator priors (`prior_mods`, `prior_mods_null`), you can specify different
#' priors for different predictors using named lists:
#' \preformatted{
#' RoBMA(...,
#'   mods = ~ x1 + x2,
#'   prior_mods_null = list(x1 = NULL)  # No null for x1, default null for x2
#' )
#' }
#' This allows testing whether specific moderators have effects while assuming others
#' are always included or excluded.
#'
#' ## Mixture priors for different parameter types
#'
#' ### Effect size (\eqn{\mu}) and heterogeneity (\eqn{\tau})
#'
#' Effect and heterogeneity priors combine spike-and-slab components:
#' \preformatted{
#' # Default: spike(0) null + normal alternative
#' # Custom: multiple alternatives with different scales
#' RoBMA(...,
#'   prior_effect = list(
#'     prior("normal", list(mean = 0, sd = 0.5)),
#'     prior("normal", list(mean = 0, sd = 1.0))
#'   )
#' )
#' }
#'
#' ### Publication bias
#'
#' Bias priors combine different publication bias adjustment mechanisms:
#' \itemize{
#'   \item No publication bias (null hypothesis)
#'   \item Weight functions: Selection models with different step functions
#'   \item PET/PEESE: Regression-based bias adjustment
#' }
#'
#' The `model_type` argument provides convenient presets:
#' \preformatted{
#' RoBMA(..., model_type = "PSMA")  # Full ensemble (default)
#' RoBMA(..., model_type = "PP")    # Only PET-PEESE models
#' }
#'
#' ### Heterogeneity allocation (\eqn{\rho}) for multilevel models
#'
#' When `study_ids` is specified, \eqn{\rho} controls the within-study vs. between-study
#' variance allocation. By default, only an alternative (Beta(1,1)) is specified,
#' but null hypotheses can be added:
#' \preformatted{
#' RoBMA(...,
#'   study_ids = study,
#'   prior_heterogeneity_allocation_null = prior("spike", list(location = 0.5))
#' )
#' }
#'
#' ### Moderator and scale regression terms
#'
#' When `mods` is specified, the effect prior becomes the intercept prior, and each
#' moderator receives a mixture prior. The **intercept** inherits the effect mixture
#' structure, while **regression coefficients** get separate spike-and-slab mixtures.
#'
#' Similarly, when `scale = ~ ...` is specified for location-scale models, the
#' heterogeneity prior becomes the scale intercept prior, and each scale predictor
#' receives a mixture prior following the same pattern as moderator terms.
#' Note that the scale intercept (i.e., baseline heterogeneity) is always positive and 
#' the multiplicative scale coefficients require non-zero prior distribution on the 
#' intercept.
#' 
#' ## Model weights and prior odds
#'
#' Each component in a mixture prior has an associated prior weight (prior model
#' probability). By default, all components receive equal weights. Custom weights can
#' be specified via the `prior_weights` argument in individual prior objects.
#'
#' The prior odds between null and alternative hypotheses affect the Bayes factor
#' interpretation. With equal prior weights on null and alternative, the posterior
#' odds equal the Bayes factor.
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso
#' \code{\link{prior_specification}} for base prior specification options,
#' \code{\link{RoBMA}} for the main model-averaging function,
#' \code{\link[BayesTools]{prior}} for creating prior distribution objects
#'
#' @examples
#' \dontrun{
#' # Default mixture priors (spike-and-slab for effect and heterogeneity)
#' fit1 <- RoBMA(yi = effect, sei = se, data = dat)
#'
#' # Custom alternative priors (two effect size scales)
#' fit2 <- RoBMA(yi = effect, sei = se, data = dat,
#'   prior_effect = list(
#'     prior("normal", list(mean = 0, sd = 0.35)),
#'     prior("normal", list(mean = 0, sd = 0.70))
#'   )
#' )
#'
#' # No null hypothesis for effect (assume effect exists)
#' fit3 <- RoBMA(yi = effect, sei = se, data = dat,
#'   prior_effect_null = NULL
#' )
#'
#' # Only selection model bias adjustment (no PET-PEESE)
#' fit4 <- RoBMA(yi = effect, sei = se, data = dat,
#'   model_type = "6w"
#' )
#'
#' # Multilevel model with custom rho prior
#' fit5 <- RoBMA(yi = effect, sei = se, study_ids = study, data = dat,
#'   prior_heterogeneity_allocation = prior("beta", list(alpha = 2, beta = 2)),
#'   prior_heterogeneity_allocation_null = prior("spike", list(location = 0.5))
#' )
#'
#' # Meta-regression with different null specifications per moderator
#' fit6 <- RoBMA(yi = effect, sei = se, mods = ~ x1 + x2, data = dat,
#'   prior_mods_null = list(x1 = NULL)  # Always include x1, test x2
#' )
#' }
#'
#' @aliases RoBMA_prior_specification
NULL




# verify that the the measure is available
.check_measure <- function(measure) {

  BayesTools::check_char(measure, "measure", allow_values = c("SMD", "ZCOR", "RR", "OR", "HR", "RD", "IRR", "GEN"))

  return()
}


### prior scale functions
.get_unit_information_sd <- function(data, measure) {

  if (measure %in% c("SMD", "ZCOR", "RR", "OR", "HR", "IRR")) {
    # use the pre-specified UISD for known effect sizes
    prior_unit_information_sd <- .get_unit_information_sd.known(measure)
  } else {
    # compute sd based on the sample sizes for the general effect sizes
    if (all(is.na(data[["outcome"]][["ni"]]))){
      stop("Sample size 'ni' or unit information sd 'unit_information_sd' must be specified to set-up prior distributions without known UISD.", call. = FALSE)
    }
    prior_unit_information_sd <- .get_unit_information_sd.norm(sei = data[["outcome"]][["sei"]], ni = data[["outcome"]][["ni"]])

  }

  return(prior_unit_information_sd)
}

# returns the sd of unit information
# based on Chapter 2.4 in Spiegelhalter, Abrams, and Myles 2004 and Chapter 1 in Grieve 2022
# as reported in Table 2 of Pawel, S., & Held, L. (2025). Closed-form power and sample size calculations for Bayes factors. The American Statistician 9(3), 330-344, 10.1080/00031305.2025.2467919
.get_unit_information_sd.known <- function(measure) {

  BayesTools::check_char(measure, "measure", allow_values = c("SMD", "ZCOR", "RR", "OR", "HR", "IRR"), allow_NA = FALSE)

  return(switch(
    measure,
    "SMD"  = sqrt(2),
    "ZCOR" = 1,
    "RR"   = sqrt(4),
    "OR"   = sqrt(4),
    "HR"   = sqrt(4),
    "IRR"  = sqrt(4)
  ))
}

# computes the unit information based on data
# based on Eq. 6 of Röver, C., Bender, R., Dias, S., Schmid, C. H., Schmidli, H., Sturtz, S., ... & Friede, T. (2021). On weakly informative prior distributions for the heterogeneity parameter in Bayesian random‐effects meta‐analysis. Research Synthesis Methods, 12(4), 448-474. 10.1002/jrsm.1475
.get_unit_information_sd.norm <- function(sei, ni) {

  UISD <- sqrt(sum(ni) / sum(sei^(-2)))

  return(UISD)
}

# computes the unit information based on Poisson (rate) data
# for IRR effect size with GLMM models
.get_unit_information_sd.pois <- function(x1i, x2i, t1i, t2i) {

  # Total events represents the Total Precision (Fisher Information) for the log-rate
  total_events <- sum(x1i + x2i)

  # Total person-time represents the Total Sample Size
  total_time <- sum(t1i + t2i)

  # UISD: sqrt(Total N / Total Precision)
  # Derived from Fisher Information I(alpha) = lambda * t
  # Unit variance approx 1/(lambda * 1_unit_time) -> SD = sqrt(1/lambda) -> sqrt(t/y)
  UISD <- sqrt(total_time / total_events)

  return(UISD)
}


#' @title Estimate Unit Information Standard Deviation
#'
#' @description Estimates the unit information standard deviation (UISD) from
#' sample sizes and standard errors. The UISD is used to scale weakly informative
#' prior distributions for meta-analytic parameters.
#'
#' @param sei a numeric vector of standard errors for each study.
#' @param ni a numeric vector of sample sizes for each study (must have the same
#' length as `sei`).
#'
#' @details
#' The unit information standard deviation is computed following Equation 6 in
#' \insertCite{rover2021weakly;textual}{RoBMA}:
#' \deqn{\text{UISD} = \sqrt{\frac{\sum n_i}{\sum \text{se}_i^{-2}}}}
#' where \eqn{n_i} is the sample size and \eqn{\text{se}_i} is the standard error
#' for each study.
#'
#' This function is useful when you want to compute the UISD once and pass it
#' to multiple analyses via the `prior_unit_information_sd` argument in
#' [brma.norm()] or related functions. This ensures consistent prior scaling
#' across analyses performed on different subsets of the same data.
#'
#' @references
#' \insertAllCited{}
#'
#' @return Returns a single numeric value representing the estimated unit
#' information standard deviation.
#'
#' @examples
#' # Example with simulated data
#' sei <- c(0.2, 0.3, 0.25, 0.15)
#' ni  <- c(50, 30, 40, 80)
#' estimate_unit_information_sd(sei = sei, ni = ni)
#'
#' @export
estimate_unit_information_sd <- function(sei, ni) {

  BayesTools::check_real(sei, "sei", lower = 0, allow_bound = FALSE, allow_NULL = FALSE, allow_NA = FALSE, check_length = FALSE)
  BayesTools::check_real(ni,  "ni",  lower = 0, allow_bound = FALSE, allow_NULL = FALSE, allow_NA = FALSE, check_length = length(sei))

  return(.get_unit_information_sd.norm(sei = sei, ni = ni))
}

### assign prior distributions
.assign_prior.simple                   <- function(
    prior, parameter, measure, data, prior_unit_information_sd,
    prior_informed_field, prior_informed_subfield, rescale_priors) {


  ### always use the user-specified prior if possible
  if (!missing(prior)) {

    # allow using FALSE/NULL arguments to set spike(0) prior distributions
    if (is.null(prior) || isFALSE(prior)) {
      return(BayesTools::prior("spike", parameters = list(0)))
    }

    # apply rescaling if specified
    prior <- .rescale_prior_object(prior, rescale_priors)

    # check the specified prior distribution
    switch(
      parameter,
      "effect"        = .check_prior.effect(prior),
      "heterogeneity" = .check_prior.heterogeneity(prior)
    )

    return(prior)
  }

  ### use informed prior distributions if specified
  if (!missing(prior_informed_field)) {

    if (prior_informed_field == "medicine") {

      # input translation for BayesTools::prior_informed
      if (missing(prior_informed_subfield)) {
        prior_name <- "Cochrane"
      } else {
        prior_name <- prior_informed_subfield
      }

      type <- switch(
        tolower(measure),
        "smd"   = "smd",
        "or"    = "logOR",
        "rr"    = "logRR",
        "rd"    = "RD",
        "hr"    = "logHR"
      )

      # simple informed priors from BayesTools package
      prior <- BayesTools::prior_informed(prior_name, parameter = parameter, type = type)
    }

    # apply rescaling if specified
    prior <- .rescale_prior_object(prior, rescale_priors)
    return(prior)
  }

  ### use prior_unit_information_sd otherwise
  if (missing(prior_unit_information_sd)) {
    # if not passed directly, obtain from data
    prior_unit_information_sd <- .get_unit_information_sd(data = data, measure = measure)
  }

  # get the prior standard deviation based on unit information
  prior_sd <- switch(
    parameter,
    "effect"        = prior_unit_information_sd * RoBMA.get_option("default_UISD.effect"),
    "heterogeneity" = prior_unit_information_sd * RoBMA.get_option("default_UISD.heterogeneity")
  )

  # create to corresponding prior
  prior <- switch(
    parameter,
    "effect"        = BayesTools::prior("normal", parameters = list("mean" = 0, "sd" = prior_sd)),
    "heterogeneity" = BayesTools::prior("normal", parameters = list("mean" = 0, "sd" = prior_sd), truncation = list(0, Inf))
  )

  # apply rescaling if specified
  prior <- .rescale_prior_object(prior, rescale_priors)

  return(prior)
}
.assign_prior.heterogeneity_allocation <- function(prior) {

  # use default settings if prior distribution is not specified
  if (missing(prior)) {
    prior <- BayesTools::prior("beta", parameters = list("alpha" = 1, "beta" = 1))
    return(prior)
  }

  # check the user specified prior distribution
  .check_prior.restricted_01(prior, prior_name = "heterogeneity_allocation")

  return(prior)
}
.assign_prior.baserate                 <- function(prior) {

  # use default settings if prior distribution is not specified
  if (missing(prior)) {
    prior <- BayesTools::prior_factor("beta", parameters = list("alpha" = 1, "beta" = 1), contrast = "independent")
    return(prior)
  }

  # check the user specified prior distribution
  .check_prior.restricted_01(prior, prior_name = "baserate")

  # transform the prior into independent factor prior
  prior <- BayesTools::prior_factor(prior[["distribution"]], parameters = prior[["parameters"]], truncation = prior[["truncation"]], contrast = "independent")

  return(prior)
}
.assign_prior.lograte                  <- function(prior, data) {
  # an uninformative prior distribution is not possible in JAGS (the log-rate is unbounded)
  # therefore, we set independent unit information priors on the log-rate if user unspecified

  if (missing(prior)) {
    return(.get_unit_information_lograte_prior(
      x1i = data[["outcome"]][["x1i"]],
      x2i = data[["outcome"]][["x2i"]],
      t1i = data[["outcome"]][["t1i"]],
      t2i = data[["outcome"]][["t2i"]])
    )
  }

  # transform the prior into independent factor prior
  prior <- BayesTools::prior_factor(prior[["distribution"]], parameters = prior[["parameters"]], truncation = prior[["truncation"]], contrast = "independent")

  return(prior)
}
.assign_prior.bias                     <- function(prior, measure, data, prior_unit_information_sd, bias_type, steps) {

  # assigns either a weight function or PET/PEESE model based on bias_type
  # there is no scaling of the publication bias priors based on rescale_priors
  # (only PEESE prior is rescaled to the match the UISD of the measure)

  if (bias_type == "selmodel") {

    # specify default settings / check the prior distribution
    if (missing(prior) && missing(steps)) {
      prior <- BayesTools::prior_weightfunction("one-sided", parameters = list("steps" = c(0.025), "alpha" = rep(RoBMA.get_option("default_bias_weightfunction.alpha"), 2)))
    } else if (missing(prior) && !missing(steps)) {
      prior <- BayesTools::prior_weightfunction("one-sided", parameters = list("steps" = steps, "alpha" = rep(RoBMA.get_option("default_bias_weightfunction.alpha"), length(steps) + 1)))
    } else if (!BayesTools::is.prior.weightfunction(prior)) {
      stop("'prior_bias' must be a `prior_weightfunction` object", call. = FALSE)
    }

    return(prior)
  }

  if (bias_type == "PET") {

    # specify default settings / check the prior distribution
    if (missing(prior)) {
      # PET prior is invariant to different effect size types
      prior <- BayesTools::prior_PET("cauchy", parameters = list("location" = 0, "scale" = RoBMA.get_option("default_bias_PET.scale")), truncation = list(0, Inf))
    } else if (!BayesTools::is.prior.PET(prior)) {
      stop("'prior_bias' must be a `prior_PET` object", call. = FALSE)
    }

    return(prior)
  }

  if (bias_type == "PEESE") {

    # specify default settings / check the prior distribution
    if (missing(prior)) {
      ### PEESE prior is inversely related to the specified effect size
      # for SMD we use Cauchy with scale 1
      # for other effect sizes, we re-scale by the corresponding UISD
      UISD_SMD     <- .get_unit_information_sd(measure = "SMD")
      UISD_measure <- .get_unit_information_sd(data = data, measure = measure)
      UISD_ratio   <- UISD_SMD / UISD_measure
      prior <- BayesTools::prior_PEESE("cauchy", parameters = list("location" = 0, "scale" = RoBMA.get_option("default_bias_PEESE.scale") * UISD_ratio), truncation = list(0, Inf))
    } else if (!BayesTools::is.prior.PEESE(prior)) {
      stop("'prior_bias' must be a `prior_PEESE` object", call. = FALSE)
    }

    return(prior)
  }

  if (bias_type == "none") {
    # for simplified handling in RoBMA
    return(prior)
  }
}
.assign_prior_list.terms               <- function(
    prior_list, prior_intercept, parameter, measure, data, prior_unit_information_sd,
    prior_informed_field, prior_informed_subfield, rescale_priors,
    set_default_spike = FALSE) {

  # set_default_spike argument is used when calling from .assign_prior_list.terms_mixture
  # to set-up prior distributions under the null hypotheses

  ### for scale regression, the intercept cannot be spike(0) since it's used as log(intercept)
  # in exp( log(intercept) + sum(beta_i * x_i) )
  if (parameter == "scale") {
    if (BayesTools::is.prior.point(prior_intercept) && mean(prior_intercept) == 0) {
      stop("Intercept prior distribution for scale regression cannot be set to 0 since the model is parameterized as tau = exp( log(intercept) + sum(beta_i * x_i) ).", call. = FALSE)
    }
  }

  ### check the user-specified priors
  if (!missing(prior_list)) {
    # check the priors and apply rescaling if specified
    for (i in seq_along(prior_list)) {
      .check_prior.term(prior_list[[i]])
      prior_list[[i]] <- .rescale_prior_object(prior_list[[i]], rescale_priors)
    }
  } else {
    # create an empty placeholder instead
    prior_list <- list()
  }

  ### add intercept to be beginning of the list (parsed previously)
  prior_list <- c(list(intercept = prior_intercept), prior_list)

  ### add default priors for remaining terms if left unspecified
  # the assignment of missing priors is performed within BayesTools::JAGS_formula
  # (it's easier to check and fill in the formula on assignment when it's being parsed)

  # use informed prior distributions if specified
  # otherwise rely on unit information based priors
  if (set_default_spike) {

    # set_default_spike argument is used when calling from .assign_prior_list.terms_mixture
    # to set-up prior distributions under the null hypotheses
    default_prior_continuous <- BayesTools::prior("spike", parameters = list(0))
    default_prior_factor     <- BayesTools::prior_factor("spike", parameters = list(0), contrast = attr(data, "set_contrast_factor_predictors"))

  } else if (!missing(prior_informed_field)) {

    if (prior_informed_field == "medicine") {

      # input translation for BayesTools::prior_informed
      if (missing(prior_informed_subfield)) {
        prior_name <- "Cochrane"
      } else {
        prior_name <- prior_informed_subfield
      }

      type <- switch(
        tolower(measure),
        "smd"   = "smd",
        "or"    = "logOR",
        "rr"    = "logRR",
        "rd"    = "RD",
        "hr"    = "logHR"
      )

      if (parameter == "mods") {
        # treat moderators by scaling them according to 1/2 of the effect
        default_prior_continuous <- BayesTools::prior_informed(prior_name, parameter = "effect", type = type)
        default_prior_continuous <- .rescale_prior_object(default_prior_continuous, RoBMA.get_option("default_informed_priors.mods"))
      } else if (parameter == "scale") {
        # scale parameters use the default scaling (no informed priors exist)
        default_prior_continuous <- BayesTools::prior("normal", parameters = list("mean" = 0, "sd" = .rescale_prior_object(prior, RoBMA.get_option("default_informed_priors.scale"))))
      }

      # transform the continuous prior distributions into prior distributions for factors
      if (attr(data, "set_contrast_factor_predictors") == "treatment") {
        default_prior_factor <- BayesTools::prior_factor(
          distribution = default_prior_continuous[["distribution"]],
          parameters   = default_prior_continuous[["parameters"]],
          contrast     = "treatment"
        )
      } else {
        default_prior_factor <- BayesTools::prior_factor(
          distribution = paste0("m", default_prior_continuous[["distribution"]]),
          parameters   = default_prior_continuous[["parameters"]],
          contrast     = attr(data, "set_contrast_factor_predictors")
        )
      }

    }

  } else {

    ### use prior_unit_information_sd otherwise
    if (missing(prior_unit_information_sd)) {
      # if not passed directly, obtain from data
      prior_unit_information_sd <- .get_unit_information_sd(data = data, measure = measure)
    }

    # get the prior standard deviation based on unit information
    prior_sd <- switch(
      parameter,
      "mods"          = prior_unit_information_sd * RoBMA.get_option("default_UISD.mods"),
      "scale"         = RoBMA.get_option("default_UISD.scale") # moderators for scale do not depend on UISD (multiplicative)
    )

    default_prior_continuous <- BayesTools::prior("normal", parameters = list("mean" = 0, "sd" = prior_sd))
    if (attr(data, "set_contrast_factor_predictors") == "treatment") {
      default_prior_factor <- BayesTools::prior_factor(
        distribution = "normal",
        parameters   = list("mean" = 0, "sd" = prior_sd),
        contrast     = "treatment"
      )
    } else {
      default_prior_factor <- BayesTools::prior_factor(
        distribution = "mnormal",
        parameters   = list("mean" = 0, "sd" = prior_sd),
        contrast     = attr(data, "set_contrast_factor_predictors")
      )
    }
  }

  # apply rescaling if specified
  default_prior_continuous <- .rescale_prior_object(default_prior_continuous, rescale_priors)
  default_prior_factor     <- .rescale_prior_object(default_prior_factor, rescale_priors)

  # list the priors into the default positions
  prior_list[["__default_continuous"]] <- default_prior_continuous
  prior_list[["__default_factor"]]     <- default_prior_factor

  ### use BayesTools::JAGS_formula to obtain the assigned list of priors
  prior_list <- BayesTools::JAGS_formula(
    parameter  = switch (parameter, "mods" = "mu", "scale" = "log_tau"),
    data       = data[[parameter]],
    formula    = attr(data[[parameter]], "formula"),
    prior_list = prior_list
  )[["prior_list"]]
  names(prior_list) <- BayesTools::format_parameter_names(
    parameters         = names(prior_list),
    formula_parameters = switch (parameter, "mods" = "mu", "scale" = "log_tau"),
    formula_prefix     = FALSE
  )

  return(prior_list)
}
### assign mixture prior distributions
.assign_prior.simple_mixture                   <- function(
    prior, prior_null, parameter, measure, data, prior_unit_information_sd,
    prior_informed_field, prior_informed_subfield, rescale_priors) {

  ### check the alternative prior input
  # prior can be either
  # - unspecified (set default as in brma)
  # - NULL/FALSE (omit prior distribution)
  # - single prior distribution
  # - list of prior distributions
  # each of these prior distributions needs to satisfy the same properties as in brma
  # if only a single distribution is specified, place it in a list for a simpler treatment
  if (missing(prior) || BayesTools::is.prior(prior)) {
    priors_alt <- list(.assign_prior.simple(
      prior = prior, parameter = parameter, measure = measure,
      data = data, prior_unit_information_sd = prior_unit_information_sd,
      prior_informed_field = prior_informed_field, prior_informed_subfield = prior_informed_subfield,
      rescale_priors = rescale_priors
    ))
  } else if (is.null(prior) || isFALSE(prior)) {
    priors_alt <- list()
  } else if (is.list(prior) && all(sapply(prior, BayesTools::is.prior))) {
    priors_alt <- prior
    for (i in seq_along(priors_alt)) {
      priors_alt[[i]] <- .assign_prior.simple(
        prior = priors_alt[[i]], parameter = parameter, measure = measure,
        data = data, prior_unit_information_sd = prior_unit_information_sd,
        prior_informed_field = prior_informed_field, prior_informed_subfield = prior_informed_subfield,
        rescale_priors = rescale_priors
      )
    }
  } else {
    stop(gettextf("The 'prior_%1$s' must be either prior distribution or a list of prior distributions.", parameter), call. = FALSE)
  }

  ### check the null prior input
  # prior can be either
  # - empty (set point prior as default)
  # - NULL/FALSE (omit prior distribution)
  # - single prior distribution
  # - list of prior distributions
  # each of these prior distributions needs to satisfy the same properties as in brma
  # if only a single distribution is specified, place it in a list for a simpler treatment
  if (missing(prior_null)) {
    priors_null <- list(BayesTools::prior("spike", parameters = list(0)))
  } else if (is.null(prior_null) || isFALSE(prior_null)) {
    priors_null <- list()
  } else if (BayesTools::is.prior(prior_null)) {
    priors_null <- list(.assign_prior.simple(
      prior = prior_null, parameter = parameter, measure = measure,
      data = data, prior_unit_information_sd = prior_unit_information_sd,
      prior_informed_field = prior_informed_field, prior_informed_subfield = prior_informed_subfield,
      rescale_priors = rescale_priors
    ))
  } else if (is.list(prior_null) && all(sapply(prior_null, BayesTools::is.prior))) {
    priors_null <- prior_null
    for (i in seq_along(priors_null)) {
      priors_null[[i]] <- .assign_prior.simple(
        prior = priors_null[[i]], parameter = parameter, measure = measure,
        data = data, prior_unit_information_sd = prior_unit_information_sd,
        prior_informed_field = prior_informed_field, prior_informed_subfield = prior_informed_subfield,
        rescale_priors = rescale_priors
      )
    }
  } else {
    stop(gettextf("The 'prior_%1$s_null' must be either prior distribution or a list of prior distributions.", parameter), call. = FALSE)
  }

  ### specify the corresponding mixture prior
  prior <- BayesTools::prior_mixture(
    prior_list = c(priors_null, priors_alt),
    is_null    = c(rep(TRUE, length(priors_null)), rep(FALSE, length(priors_alt)))
  )

  return(prior)
}
.assign_prior.bias_mixture                     <- function(
    prior, prior_null, measure, data, prior_unit_information_sd, model_type) {

  ### use precanned publication bias adjustment models
  # if neither `prior` or `prior_null` input is specified use the `model_type`
  # there is no scaling of the publication bias priors based on rescale_priors
  # (only PEESE prior is rescaled to the match the UISD of the measure)

  if (missing(prior) && missing(prior_null)) {

    ### PEESE prior is inversely related to the specified effect size
    # for SMD we use Cauchy with scale 1
    # for other effect sizes, we re-scale by the corresponding UISD
    UISD_SMD <- .get_unit_information_sd(measure = "SMD")
    if (missing(prior_unit_information_sd)) {
      UISD_measure <- .get_unit_information_sd(data = data, measure = measure)
    } else {
      UISD_measure <- prior_unit_information_sd
    }
    UISD_ratio <- UISD_SMD / UISD_measure

    if (model_type == "2w") {
      # prior based on Maier 2022 original RoBMA article
      prior <- BayesTools::prior_mixture(list(
        BayesTools::prior_none(prior_weights = 1),
        BayesTools::prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1),       steps = c(0.05)),             prior_weights = 1/2),
        BayesTools::prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.05, 0.10)),       prior_weights = 1/2)
      ), is_null = c(TRUE, rep(FALSE, 2))
      )
    } else if (model_type == "6w") {
      # subset of selection models used in Bartos 2022 extended RoBMA article
      prior <- BayesTools::prior_mixture(list(
        BayesTools::prior_none(prior_weights = 1),
        BayesTools::prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1),       steps = c(0.05)),             prior_weights = 1/6),
        BayesTools::prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.05, 0.10)),       prior_weights = 1/6),
        BayesTools::prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1),       steps = c(0.05)),             prior_weights = 1/6),
        BayesTools::prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.025, 0.05)),      prior_weights = 1/6),
        BayesTools::prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.05, 0.5)),        prior_weights = 1/6),
        BayesTools::prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1, 1), steps = c(0.025, 0.05, 0.5)), prior_weights = 1/6)
      ), is_null = c(TRUE, rep(FALSE, 6))
      )
    } else if (model_type == "PP") {
      # subset of PET-PEESE models used in Bartos 2022 extended RoBMA article
      prior <- BayesTools::prior_mixture(list(
        BayesTools::prior_none(prior_weights = 1),
        BayesTools::prior_PET(distribution = "Cauchy",   parameters = list(0, RoBMA.get_option("default_bias_PET.scale")),                truncation = list(0, Inf),   prior_weights = 1/2),
        BayesTools::prior_PEESE(distribution = "Cauchy", parameters = list(0, RoBMA.get_option("default_bias_PEESE.scale") * UISD_ratio), truncation = list(0, Inf),   prior_weights = 1/2)
      ), is_null = c(TRUE, rep(FALSE, 2))
      )
    } else if (model_type == "PSMA") {
      # all models used in Bartos 2022 extended RoBMA article
      prior <- BayesTools::prior_mixture(list(
        BayesTools::prior_none(prior_weights = 1),
        BayesTools::prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1),       steps = c(0.05)),             prior_weights = 1/12),
        BayesTools::prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.05, 0.10)),       prior_weights = 1/12),
        BayesTools::prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1),       steps = c(0.05)),             prior_weights = 1/12),
        BayesTools::prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.025, 0.05)),      prior_weights = 1/12),
        BayesTools::prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.05, 0.5)),        prior_weights = 1/12),
        BayesTools::prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1, 1), steps = c(0.025, 0.05, 0.5)), prior_weights = 1/12),
        BayesTools::prior_PET(distribution = "Cauchy",   parameters = list(0, RoBMA.get_option("default_bias_PET.scale")),                truncation = list(0, Inf),   prior_weights = 1/4),
        BayesTools::prior_PEESE(distribution = "Cauchy", parameters = list(0, RoBMA.get_option("default_bias_PEESE.scale") * UISD_ratio), truncation = list(0, Inf),   prior_weights = 1/4)
      ), is_null = c(TRUE, rep(FALSE, 8))
      )
    }

    return(prior)
  }


  ### use the user specified prior distributions if specified

  ### check the alternative prior input
  # prior can be either
  # - unspecified (set alternative part of PSMA)
  # - NULL/FALSE (omit prior distribution)
  # - single prior distribution
  # - list of prior distributions
  # each of these prior distributions needs to satisfy the same properties as in brma
  # if only a single distribution is specified, place it in a list for a simpler treatment
  if (missing(prior)) {
    priors_alt <- list(
      BayesTools::prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1),       steps = c(0.05)),             prior_weights = 1/12),
      BayesTools::prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.05, 0.10)),       prior_weights = 1/12),
      BayesTools::prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1),       steps = c(0.05)),             prior_weights = 1/12),
      BayesTools::prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.025, 0.05)),      prior_weights = 1/12),
      BayesTools::prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.05, 0.5)),        prior_weights = 1/12),
      BayesTools::prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1, 1), steps = c(0.025, 0.05, 0.5)), prior_weights = 1/12),
      BayesTools::prior_PET(distribution = "Cauchy",   parameters = list(0, RoBMA.get_option("default_bias_PET.scale")),                truncation = list(0, Inf),   prior_weights = 1/4),
      BayesTools::prior_PEESE(distribution = "Cauchy", parameters = list(0, RoBMA.get_option("default_bias_PEESE.scale") * UISD_ratio), truncation = list(0, Inf),   prior_weights = 1/4)
    )
  } else if (is.null(prior) || isFALSE(prior)) {
    priors_alt <- list()
  } else if (BayesTools::is.prior(prior)) {
    priors_alt <- list(.assign_prior.bias(
      prior = prior, measure = measure, data = data,
      prior_unit_information_sd = prior_unit_information_sd, bias_type = .get_prior_bias_type(prior)))
  } else if (is.list(prior) && all(sapply(prior, BayesTools::is.prior))) {
    priors_alt <- prior
    for (i in seq_along(priors_alt)) {
      priors_alt[[i]] <- .assign_prior.bias(
        prior = priors_alt[[i]], measure = measure, data = data,
        prior_unit_information_sd = prior_unit_information_sd, bias_type = .get_prior_bias_type(priors_alt[[i]]))
    }
  } else {
    stop("The 'prior_bias' must be either prior distribution or a list of prior distributions.", call. = FALSE)
  }

  ### check the null prior input
  # prior can be either
  # - empty (set point prior as default)
  # - NULL/FALSE (omit prior distribution)
  # - single prior distribution
  # - list of prior distributions
  # each of these prior distributions needs to satisfy the same properties as in brma
  # if only a single distribution is specified, place it in a list for a simpler treatment
  if (missing(prior_null)) {
    priors_null <- list(BayesTools::prior_none(prior_weights = 1))
  } else if (is.null(prior_null) || isFALSE(prior_null)) {
    priors_null <- list()
  } else if (BayesTools::is.prior(prior_null)) {
    priors_null <- list(.assign_prior.bias(
      prior = prior_null, measure = measure, data = data,
      prior_unit_information_sd = prior_unit_information_sd, bias_type = .get_prior_bias_type(prior_null)))
  } else if (is.list(prior_null) && all(sapply(prior_null, BayesTools::is.prior))) {
    priors_null <- prior_null
    for (i in seq_along(priors_null)) {
      priors_null[[i]] <- .assign_prior.bias(
        prior = priors_null[[i]], measure = measure, data = data,
        prior_unit_information_sd = prior_unit_information_sd, bias_type = .get_prior_bias_type(priors_null[[i]]))
    }
  } else {
    stop("The 'prior_bias_null' must be either prior distribution or a list of prior distributions.", call. = FALSE)
  }

  ### specify the corresponding mixture prior
  prior <- BayesTools::prior_mixture(
    prior_list = c(priors_null, priors_alt),
    is_null    = c(rep(TRUE, length(priors_null)), rep(FALSE, length(priors_alt)))
  )

  return(prior)
}
.assign_prior.heterogeneity_allocation_mixture <- function(prior, prior_null) {

  ### check the alternative prior input
  # prior can be either
  # - unspecified (set default as in brma)
  # - NULL/FALSE (omit prior distribution)
  # - single prior distribution
  # - list of prior distributions
  # each of these prior distributions needs to satisfy the same properties as in brma
  # if only a single distribution is specified, place it in a list for a simpler treatment
  if (missing(prior) || BayesTools::is.prior(prior)) {
    priors_alt <- list(.assign_prior.heterogeneity_allocation(prior = prior))
  } else if (is.null(prior) || isFALSE(prior)) {
    priors_alt <- list()
  } else if (is.list(prior) && all(sapply(prior, BayesTools::is.prior))) {
    priors_alt <- prior
    for (i in seq_along(priors_alt)) {
      priors_alt[[i]] <- .assign_prior.heterogeneity_allocation(prior = priors_alt[[i]])
    }
  } else {
    stop("The 'prior_heterogeneity_allocation' must be either prior distribution or a list of prior distributions.", call. = FALSE)
  }

  ### check the null prior input
  # prior can be either
  # - empty (omit prior as default)
  # - NULL/FALSE (omit prior distribution)
  # - single prior distribution
  # - list of prior distributions
  # each of these prior distributions needs to satisfy the same properties as in brma
  # if only a single distribution is specified, place it in a list for a simpler treatment
  if (missing(prior_null) || is.null(prior_null) || isFALSE(prior_null)) {
    priors_null <- list()
  } else if (BayesTools::is.prior(prior_null)) {
    priors_null <- list(.assign_prior.heterogeneity_allocation(prior = prior_null))
  } else if (is.list(prior_null) && all(sapply(prior_null, BayesTools::is.prior))) {
    priors_null <- prior_null
    for (i in seq_along(priors_null)) {
      priors_null[[i]] <- .assign_prior.heterogeneity_allocation(prior = priors_null[[i]])
    }
  } else {
    stop("The 'prior_heterogeneity_allocation_null' must be either prior distribution or a list of prior distributions.", call. = FALSE)
  }

  ### specify the corresponding mixture prior
  prior <- BayesTools::prior_mixture(
    prior_list = c(priors_null, priors_alt),
    is_null    = c(rep(TRUE, length(priors_null)), rep(FALSE, length(priors_alt)))
  )

  return(prior)
}
.assign_prior_list.terms_mixture               <- function(
    prior_list, prior_list_null, prior_intercept, parameter, measure, data, prior_unit_information_sd,
    prior_informed_field, prior_informed_subfield, rescale_priors) {

  # in contrast to effect/heterogeneity prior list settings terms currently only allow a single prior for alternative
  # and a single prior for the null hypothesis for each term

  # one issue when re-using .assign_prior_list.terms is that we need to forward prior for the intercept which is
  # an already parsed mixture prior - when merging mixture priors for terms, we ought to not duplicate the intercept

  # we want to allow NULL/FALSE to omit prior distribution for a given component, .assign_prior_list.terms
  # fills all missing prior distributions with default priors
  # as such, we store names of parameters that should be omited and remove them manually
  if (!missing(prior_list) && length(prior_list) > 0 && is.list(prior_list)) {
    par_alt_omit <- names(prior_list)[sapply(prior_list, function(x) is.null(x) || isFALSE(x))]
    # filter out NULL/FALSE values before passing to .assign_prior_list.terms
    prior_list <- prior_list[!sapply(prior_list, function(x) is.null(x) || isFALSE(x))]
  } else {
    par_alt_omit <- NULL
  }
  if (!missing(prior_list_null) && length(prior_list_null) > 0 && is.list(prior_list_null)) {
    par_null_omit <- names(prior_list_null)[sapply(prior_list_null, function(x) is.null(x) || isFALSE(x))]
    # filter out NULL/FALSE values before passing to .assign_prior_list.terms
    prior_list_null <- prior_list_null[!sapply(prior_list_null, function(x) is.null(x) || isFALSE(x))]
  } else {
    par_null_omit <- NULL
  }


  if (parameter == "scale") {
    # in the case of scale moderation, spike at zero prior (i.e., fixed-effect model) is not possible since the intercept
    # is on a log scale
    # go over the prior and remove any spike(0) priors & throw a warning
    for (i in rev(seq_along(prior_intercept))) {
      if (BayesTools::is.prior.point(prior_intercept[[i]]) && mean(prior_intercept[[i]]) == 0) {
        prior_intercept[[i]] <- NULL
        warning("Spike(0) prior distribution for between-study heterogeneity tau was removed since scale regression requires non-zero between-study heterogeneity.", immediate. = TRUE)
      }
    }
  }

  prior_list_alt <- .assign_prior_list.terms(
    prior_list = prior_list, prior_intercept = prior_intercept, parameter = parameter,
    measure = measure, data = data,  prior_unit_information_sd = prior_unit_information_sd,
    prior_informed_field = prior_informed_field, prior_informed_subfield = prior_informed_subfield,
    rescale_priors = rescale_priors)

  prior_list_null <- .assign_prior_list.terms(
    prior_list = prior_list_null, prior_intercept = prior_intercept, parameter = parameter,
    measure = measure, data = data,  prior_unit_information_sd = prior_unit_information_sd,
    prior_informed_field = prior_informed_field, prior_informed_subfield = prior_informed_subfield,
    rescale_priors = rescale_priors)

  ### initiate the list of prior mixtures with the intercept
  prior_list <- list(intercept = prior_intercept)
  for (par in names(prior_list_alt)) {

    # skip the already assigned intercept
    if (par == "intercept") {
      next
    }

    # create temporal holder for the parameter
    temp_prior_list  <- list()
    temp_is_null     <- NULL

    if (!par %in% par_null_omit) {
      temp_prior_list <- c(temp_prior_list, list(prior_list_null[[par]]))
      temp_is_null    <- c(temp_is_null, TRUE)
    }

    if (!par %in% par_alt_omit) {
      temp_prior_list <- c(temp_prior_list, list(prior_list_alt[[par]]))
      temp_is_null    <- c(temp_is_null, FALSE)
    }

    if (length(temp_prior_list) == 0) {
      stop(sprintf("At least one prior distribution needs to be defined for '%1$s' parameter in the '%2$s' model.", par, parameter))
    }

    prior_list[[par]] <- BayesTools::prior_mixture(
      prior_list = temp_prior_list,
      is_null    = temp_is_null
    )
  }

  return(prior_list)
}
### check specified prior distributions whether they follow basic requirements
.check_prior.effect                   <- function(prior) {

  # check object type
  if (!BayesTools::is.prior(prior))
    stop("The 'prior_effect' argument must be a prior distribution.", call. = FALSE)
  if (!(BayesTools::is.prior.simple(prior) || !BayesTools::is.prior.point(prior)))
    stop("The 'prior_effect' argument must be a either a simple or point prior distribution.", call. = FALSE)

  return()
}
.check_prior.heterogeneity            <- function(prior) {

  # check object type
  if (!BayesTools::is.prior(prior))
    stop("The 'prior_heterogeneity' argument must be a prior distribution.", call. = FALSE)
  if (!(BayesTools::is.prior.simple(prior) || !BayesTools::is.prior.point(prior)))
    stop("The 'prior_heterogeneity' argument must be a either a simple or point prior distribution.", call. = FALSE)

  # check range restriction
  if (prior[["truncation"]][["lower"]] < 0) {
    warning("Lower truncation point for 'prior_heterogeneity' prior distribution must be at least 0. The prior distribution was modified.",
            immediate. = TRUE, call. = FALSE)
    prior[["truncation"]][["lower"]] <- 0
  }

  return()
}
.check_prior.restricted_01            <- function(prior, prior_name) {

  # check object type
  if (!BayesTools::is.prior(prior))
    stop(sprintf("The '%s' argument must be a prior distribution.", prior_name), call. = FALSE)
  if (!(BayesTools::is.prior.simple(prior) || !BayesTools::is.prior.point(prior)))
    stop(sprintf("The '%s' argument must be a either a simple or point prior distribution.", prior_name), call. = FALSE)

  # check range restriction
  if (prior[["truncation"]][["lower"]] < 0) {
    warning(sprintf("Lower truncation point for '%s' prior distribution must be at least 0. The prior distribution was modified.", prior_name),
            immediate. = TRUE, call. = FALSE)
    prior[["truncation"]][["lower"]] <- 0
  }

  if (prior[["truncation"]][["upper"]] > 1) {
    warning(sprintf("Upper truncation point for '%s' prior distribution must be at most 1. The prior distribution was modified.", prior_name),
            immediate. = TRUE, call. = FALSE)
    prior[["truncation"]][["upper"]] <- 1
  }

  return()
}
.check_prior.term                     <- function(prior) {

  # check object type
  if (!BayesTools::is.prior(prior))
    stop("The 'prior_mods' / 'prior_scale' arguments must be a prior distribution.", call. = FALSE)
  if (!(BayesTools::is.prior.simple(prior) || BayesTools::is.prior.point(prior) || BayesTools::is.prior.factor(prior)))
    stop("The 'prior_mods' / 'prior_scale' arguments for predictors must be a either a simple, factor, or point prior distribution.", call. = FALSE)

  return()
}
.check_prior_specification_conflict   <- function(prior_unit_information_sd, prior_informed_field) {

  # check that both UISD and informed priors are not specified simultaneously
  if (!missing(prior_unit_information_sd) && !missing(prior_informed_field)) {
    stop(paste0(
      "Both 'prior_unit_information_sd' and 'prior_informed_field' were specified. ",
      "These arguments are mutually exclusive: 'prior_unit_information_sd' is used to scale default priors, ",
      "while 'prior_informed_field' uses pre-specified informed prior distributions that do not rely on UISD scaling."
    ), call. = FALSE)
  }

  return()
}
.get_prior_bias_type                  <- function(prior) {
  if (BayesTools::is.prior.none(prior)) {
    priors_alt <- "none"
  } else if (BayesTools::is.prior.weightfunction(prior)) {
    bias_type <- "selmodel"
  } else if (BayesTools::is.prior.PET(prior)) {
    bias_type <- "PET"
  } else if (BayesTools::is.prior.PEESE(prior)) {
    bias_type <- "PEESE"
  } else {
    stop("The publication bias prior was not recognized", call. = FALSE)
  }
  return(bias_type)
}

# helper function for defining unit information prior on lograte for IRR models
.get_unit_information_lograte_prior   <- function(x1i, x2i, t1i, t2i) {

  # compute the unit information standard deviation from Poisson data
  sigma_prior <- .get_unit_information_sd.pois(x1i, x2i, t1i, t2i)

  # compute the pooled crude log-rate for the mean
  total_events <- sum(x1i + x2i)
  total_time   <- sum(t1i + t2i)
  mu_prior     <- log(total_events / total_time)

  # define the prior object
  prior <- BayesTools::prior_factor("normal", parameters = list(mean = mu_prior, sd = sigma_prior), contrast = "independent")

  return(prior)
}

# priors rescaling function
.rescale_prior_object <- function(prior, rescale_priors) {
  # re-scaling prior distributions

  # no need if re-scaling is not specified
  if (missing(rescale_priors) || rescale_priors == 1)
    return(prior)

  # no need for point / no prior distributions
  if (is.prior.point(prior) || is.prior.none(prior))
    return(prior)

  # only specific prior distributions can be rescaled
  # (UPGRADE: once BayesTools incorporates direct scaling, extend to all kinds of priors)
  can_be_rescaled <- c("normal", "mnormal", "cauchy", "mcauchy", "t", "mt", "invgamma")
  if (!is.element(prior[["distribution"]], can_be_rescaled))
    stop(sprintf(
      paste0("The '%1$s' prior distribution cannot be rescaled using the 'rescale_priors' argument. ",
             "Only one of the following distributions can be rescaled: %2$s"),
      prior[["distribution"]],
      paste0("'", can_be_rescaled, "'", sep = ", ")
    ), call. = FALSE)

  if (prior[["distribution"]] %in% c("normal", "mnormal")) {
    prior[["parameters"]][["sd"]] <- prior[["parameters"]][["sd"]] * rescale_priors
  } else if (prior[["distribution"]] %in% c("cauchy", "mcauchy", "t", "mt", "invgamma")) {
    prior[["parameters"]][["scale"]] <- prior[["parameters"]][["scale"]] * rescale_priors
  }

  return(prior)
}

### prior specification functions
.check_and_list_priors.brma <- function(
    prior_effect, prior_heterogeneity,
    prior_mods, prior_scale,
    prior_heterogeneity_allocation, prior_bias,
    prior_baserate, prior_lograte,
    rescale_priors,
    prior_unit_information_sd,
    prior_informed_field, prior_informed_subfield,
    data, measure,
    bias_type, steps) {

  ### check input
  if (!missing(rescale_priors))
    BayesTools::check_real(rescale_priors, "rescale_priors", lower = 0, allow_bound = FALSE, allow_NA = FALSE)
  if (!missing(prior_unit_information_sd))
    BayesTools::check_real(prior_unit_information_sd, "prior_unit_information_sd", lower = 0, allow_bound = FALSE, allow_NA = FALSE)
  if (!missing(prior_informed_field))
    BayesTools::check_char(prior_informed_field, "prior_informed_field", allow_values = c("medicine"), allow_NA = FALSE)
  if (!missing(prior_informed_subfield))
    BayesTools::check_char(prior_informed_subfield, "prior_informed_subfield", allow_NA = FALSE)
  .check_prior_specification_conflict(prior_unit_information_sd, prior_informed_field)
  if (!missing(steps))
    BayesTools::check_real(steps, "steps", lower = 0, upper = 1, allow_bound = FALSE, check_length = FALSE, allow_NA = FALSE)

  ### set prior distributions
  measure       <- .data_measure(data)
  prior_outcome <- list()

  prior_outcome[["mu"]] <- .assign_prior.simple(
      prior = prior_effect, parameter = "effect", measure = measure,
      data = data, prior_unit_information_sd = prior_unit_information_sd,
      prior_informed_field = prior_informed_field, prior_informed_subfield = prior_informed_subfield,
      rescale_priors = rescale_priors
    )
  prior_outcome[["tau"]] <- .assign_prior.simple(
    prior = prior_heterogeneity, parameter = "heterogeneity", measure = measure,
    data = data, prior_unit_information_sd = prior_unit_information_sd,
    prior_informed_field = prior_informed_field, prior_informed_subfield = prior_informed_subfield,
    rescale_priors = rescale_priors
  )

  # only for multilevel models
  if (.is_data_multilevel(data)) {
    prior_outcome[["rho"]]   <- .assign_prior.heterogeneity_allocation(prior = prior_heterogeneity_allocation)
    # gamma is an auxiliary parameter for non-centered random-effects parameterization (for study-level)
    prior_outcome[["gamma"]] <- prior_factor("normal", parameters = list("mean" = 0, "sd" = 1), contrast = "independent")
  }
  # only for glmm models of OR
  if (.data_outcome_type(data) == "bin") {
    prior_outcome[["pi"]]    <- .assign_prior.baserate(prior = prior_baserate)
    # theta is an auxiliary parameter for non-centered random-effects parameterization (for estimate-level)
    # (these are marginalized in the normal model)
    prior_outcome[["theta"]] <- prior_factor("normal", parameters = list("mean" = 0, "sd" = 1), contrast = "independent")
  }
  # only for glmm models of IRR
  if (.data_outcome_type(data) == "pois") {
    prior_outcome[["phi"]] <- .assign_prior.lograte(prior = prior_lograte, data)
    # theta is an auxiliary parameter for non-centered random-effects parameterization (for estimate-level)
    # (these are marginalized in the normal model)
    prior_outcome[["theta"]] <- prior_factor("normal", parameters = list("mean" = 0, "sd" = 1), contrast = "independent")
  }
  # only for selection / PET / PEESE models
  if (!missing(bias_type)) {
    prior_outcome[["bias"]] <- .assign_prior.bias(
      prior = prior_bias, measure = measure, data = data, prior_unit_information_sd = prior_unit_information_sd,
      bias_type = bias_type, steps = steps)
  }


  if (.is_data_mods(data)) {
    prior_mods <- .assign_prior_list.terms(
      prior_list = prior_mods, prior_intercept = prior_outcome[["mu"]], parameter = "mods",
      measure = measure, data = data,  prior_unit_information_sd = prior_unit_information_sd,
      prior_informed_field = prior_informed_field, prior_informed_subfield = prior_informed_subfield,
      rescale_priors = rescale_priors)
    prior_outcome[["mu"]] <- NULL # the prior is forwarded through "mods" to the intercept
  } else {
    prior_mods <- NULL
  }

  if (.is_data_scale(data)) {
    prior_scale <- .assign_prior_list.terms(
      prior_list = prior_scale, prior_intercept = prior_outcome[["tau"]], parameter = "scale",
      measure = measure, data = data,  prior_unit_information_sd = prior_unit_information_sd,
      prior_informed_field = prior_informed_field, prior_informed_subfield = prior_informed_subfield,
      rescale_priors = rescale_priors)
    prior_outcome[["tau"]] <- NULL # the prior is forwarded through "scale" to the log(intercept)
  } else {
    prior_scale <- NULL
  }

  return(list(
    outcome  = prior_outcome,
    mods     = prior_mods,
    scale    = prior_scale
  ))
}


.check_and_list_priors.RoBMA <- function(
    prior_effect, prior_heterogeneity,
    prior_mods, prior_scale,
    prior_heterogeneity_allocation, prior_bias,
    prior_effect_null, prior_heterogeneity_null,
    prior_mods_null, prior_scale_null,
    prior_heterogeneity_allocation_null, prior_bias_null,
    prior_baserate, prior_lograte,
    rescale_priors,
    prior_unit_information_sd,
    prior_informed_field, prior_informed_subfield,
    data, model_type) {

  ### check input
  if (!missing(rescale_priors))
    BayesTools::check_real(rescale_priors, "rescale_priors", lower = 0, allow_bound = FALSE, allow_NA = FALSE)
  if (!missing(prior_unit_information_sd))
    BayesTools::check_real(prior_unit_information_sd, "prior_unit_information_sd", lower = 0, allow_bound = FALSE, allow_NA = FALSE)
  if (!missing(prior_informed_field))
    BayesTools::check_char(prior_informed_field, "prior_informed_field", allow_values = c("medicine"), allow_NA = FALSE)
  if (!missing(prior_informed_subfield))
    BayesTools::check_char(prior_informed_subfield, "prior_informed_subfield", allow_NA = FALSE)
  .check_prior_specification_conflict(prior_unit_information_sd, prior_informed_field)
  BayesTools::check_char(model_type, "model_type", allow_values = c("2w", "6w", "PP", "PSMA"))

  ### set prior distributions
  measure       <- .data_measure(data)
  prior_outcome <- list()

  prior_outcome[["mu"]] <- .assign_prior.simple_mixture(
    prior = prior_effect, prior_null = prior_effect_null, parameter = "effect", measure = measure,
    data = data, prior_unit_information_sd = prior_unit_information_sd,
    prior_informed_field = prior_informed_field, prior_informed_subfield = prior_informed_subfield,
    rescale_priors = rescale_priors
  )
  prior_outcome[["tau"]] <- .assign_prior.simple_mixture(
    prior = prior_heterogeneity, prior_null = prior_heterogeneity_null, parameter = "heterogeneity", measure = measure,
    data = data, prior_unit_information_sd = prior_unit_information_sd,
    prior_informed_field = prior_informed_field, prior_informed_subfield = prior_informed_subfield,
    rescale_priors = rescale_priors
  )
  prior_outcome[["bias"]] <- .assign_prior.bias_mixture(
    prior = prior_bias, prior_null = prior_bias_null, measure = measure,
    data = data, prior_unit_information_sd = prior_unit_information_sd, model_type = model_type
  )

  # only for multilevel models
  if (.is_data_multilevel(data)) {
    prior_outcome[["rho"]]   <- .assign_prior.heterogeneity_allocation_mixture(
      prior = prior_heterogeneity_allocation, prior_null = prior_heterogeneity_allocation_null)
    # gamma is an auxiliary parameter for non-centered random-effects parameterization (for study-level)
    prior_outcome[["gamma"]] <- prior_factor("normal", parameters = list("mean" = 0, "sd" = 1), contrast = "independent")
  }
  # only for glmm models of OR
  if (.data_outcome_type(data) == "bin") {
    prior_outcome[["pi"]]    <- .assign_prior.baserate(prior = prior_baserate)
    # theta is an auxiliary parameter for non-centered random-effects parameterization (for estimate-level)
    # (these are marginalized in the normal model)
    prior_outcome[["theta"]] <- prior_factor("normal", parameters = list("mean" = 0, "sd" = 1), contrast = "independent")
  }
  # only for glmm models of IRR
  if (.data_outcome_type(data) == "pois") {
    prior_outcome[["phi"]] <- .assign_prior.lograte(prior = prior_lograte, data)
    # theta is an auxiliary parameter for non-centered random-effects parameterization (for estimate-level)
    # (these are marginalized in the normal model)
    prior_outcome[["theta"]] <- prior_factor("normal", parameters = list("mean" = 0, "sd" = 1), contrast = "independent")
  }


  if (.is_data_mods(data)) {
    prior_mods <- .assign_prior_list.terms_mixture(
      prior_list = prior_mods, prior_list_null = prior_mods_null, prior_intercept = prior_outcome[["mu"]], parameter = "mods",
      measure = measure, data = data,  prior_unit_information_sd = prior_unit_information_sd,
      prior_informed_field = prior_informed_field, prior_informed_subfield = prior_informed_subfield,
      rescale_priors = rescale_priors)
    prior_outcome[["mu"]] <- NULL # the prior is forwarded through "mods" to the intercept
  } else {
    prior_mods <- NULL
  }

  if (.is_data_scale(data)) {
    prior_scale <- .assign_prior_list.terms_mixture(
      prior_list = prior_scale, prior_list_null = prior_scale_null, prior_intercept = prior_outcome[["tau"]], parameter = "scale",
      measure = measure, data = data,  prior_unit_information_sd = prior_unit_information_sd,
      prior_informed_field = prior_informed_field, prior_informed_subfield = prior_informed_subfield,
      rescale_priors = rescale_priors)
    prior_outcome[["tau"]] <- NULL # the prior is forwarded through "scale" to the log(intercept)
  } else {
    prior_scale <- NULL
  }

  return(list(
    outcome  = prior_outcome,
    mods     = prior_mods,
    scale    = prior_scale
  ))
}
