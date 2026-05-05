#' @title Prior specification
#' @name prior_specification
#'
#' @description
#' Prior distributions can be supplied for all model parameters. When omitted,
#' the fitting functions construct defaults from the effect-size measure, sample
#' sizes, or informed-prior settings when available.
#'
#' This typically includes the pooled effect \eqn{\mu} and between-study heterogeneity
#' \eqn{\tau}. In the case of meta-regression, the pooled effect \eqn{\mu} corresponds to
#' the intercept, and additional prior distributions for the regression coefficients
#' are required. In the case of a location-scale model, the between-study heterogeneity
#' corresponds to the intercept of the scale regression, and additional prior
#' distributions for the scale regression coefficients are required.
#'
#' @param prior_effect prior distribution for the effect size (\eqn{\mu}) parameter
#' (the intercept). If omitted, a default prior is constructed. In single-model
#' functions, explicit `NULL` or `FALSE` sets a spike at zero.
#' @param prior_heterogeneity prior distribution for the heterogeneity (\eqn{\tau})
#' parameter. If omitted, a default prior is constructed. In single-model
#' functions, explicit `NULL` or `FALSE` sets a spike at zero.
#' @param prior_mods prior distribution for the moderators (\eqn{\beta}) parameters.
#' A single prior applies to all terms; a named list can specify term-specific
#' priors. If omitted or `NULL`, default priors are used.
#' @param prior_scale prior distribution for the scale (\eqn{\delta}) parameters.
#' A single prior applies to all terms; a named list can specify term-specific
#' priors. If omitted or `NULL`, default priors are used.
#' @param prior_heterogeneity_allocation prior distribution for the fraction of
#' heterogeneity allocated to the cluster-level component in multilevel models
#' (\eqn{\rho}). If omitted or `NULL`, defaults to `Beta(1, 1)`.
#' @param prior_baserate prior distribution for the estimate-specific midpoint
#' base-rate probability in binomial GLMM models. If omitted or `NULL`, defaults
#' to independent `Beta(1, 1)` priors.
#' @param prior_lograte prior distribution for the estimate-specific midpoint
#' log-rate in Poisson GLMM models. If omitted or `NULL`, a data-based
#' unit-information normal prior is used independently for each estimate.
#' @param prior_unit_information_sd numeric. The unit information standard deviation (\eqn{\sigma_{unit}}).
#' Cannot be used together with `prior_informed_field`.
#' @param rescale_priors numeric. A scaling factor for supported prior distributions.
#' Point and none priors are unchanged. For constructors with publication-bias
#' prior distributions, `rescale_priors` does not rescale them except for the
#' default PEESE prior's UISD adjustment. Defaults to 1.
#' @param standardize_continuous_predictors logical. Whether to standardize continuous predictors.
#' Defaults to `TRUE`.
#' @param set_contrast_factor_predictors character. How to set contrast for factor predictors.
#' Defaults are constructor-specific and shown in each function usage; single-model
#' constructors use `"treatment"`, while model-averaging constructors use `"meandif"`.
#' @param prior_informed_field character. The field of the informed prior distributions.
#' Omit to use the standard default prior specification; explicit `NULL` is invalid.
#' @param prior_informed_subfield character. The subfield of the informed prior distributions.
#' Omit to use the field-specific default, such as `"Cochrane"` for
#' `prior_informed_field = "medicine"`; explicit `NULL` is invalid.
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
#' (section (1)). It cannot be combined with `prior_informed_field`.
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
#' The `rescale_priors` argument allows rescaling supported prior distributions by a
#' multiplicative factor. For example, `rescale_priors = 2` doubles the standard
#' deviations/scales of normal, Cauchy, t, and inverse-gamma prior distributions,
#' making them more diffuse. Point and none priors are unchanged. For publication-bias
#' prior distributions, see \code{\link{publication_bias_prior_specification}}.
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link{publication_bias_prior_specification}},
#' \code{\link[BayesTools]{prior}}, \code{\link{RoBMA.options}},
#' \code{\link{brma}}
#' @aliases prior_specification
NULL


#' @title Publication-bias prior specification
#' @name publication_bias_prior_specification
#'
#' @description
#' Publication-bias prior distributions are specified separately from the
#' meta-analytic prior distributions documented in
#' \code{\link{prior_specification}}. They define selection-model weights,
#' PET/PEESE regression coefficients, or the publication-bias model space in
#' \code{\link{RoBMA}}.
#'
#' @param prior_bias prior distribution for publication-bias adjustment. For
#' \code{\link{bselmodel}}, this is usually a weightfunction prior created by
#' \code{\link{prior_weightfunction}}. For \code{\link{bPET}}, use
#' \code{\link{prior_PET}}. For \code{\link{bPEESE}}, use
#' \code{\link{prior_PEESE}}. For \code{\link{RoBMA}}, this can be a single
#' publication-bias prior distribution or a list of publication-bias prior
#' distributions. In the single-model bias-adjustment constructors, omitted or
#' \code{NULL} uses the corresponding default prior distribution.
#' @param prior_bias_null prior distribution(s) for null publication-bias
#' component(s) in \code{\link{RoBMA}}, usually \code{\link{prior_none}()}.
#' A single prior distribution object creates one null component, a list creates
#' multiple null components, and \code{NULL} or \code{FALSE} omits the null
#' publication-bias component.
#' @param model_type character string specifying predefined publication-bias
#' model ensembles for \code{\link{RoBMA}}. One of \code{"PSMA"},
#' \code{"6w"}, \code{"2w"}, or \code{"PP"}.
#' @param steps numeric vector of one-sided p-value cut points for the default
#' \code{\link{bselmodel}} weightfunction prior. If \code{prior_bias} is supplied,
#' the prior distribution carries its own p-value cut points.
#'
#' @details
#' ## Single-model bias-adjustment priors
#'
#' \code{\link{bselmodel}}, \code{\link{bPET}}, and \code{\link{bPEESE}} fit one
#' publication-bias adjustment at a time. The \code{prior_bias} argument must
#' match the fitted bias-adjustment type.
#'
#' \tabular{lll}{
#' Constructor \tab Prior constructor \tab Default prior distribution \cr
#' \code{bselmodel()} \tab \code{\link{prior_weightfunction}} \tab one-sided cumulative weightfunction with \code{steps = 0.025} \cr
#' \code{bPET()} \tab \code{\link{prior_PET}} \tab positive Cauchy centered at 0 \cr
#' \code{bPEESE()} \tab \code{\link{prior_PEESE}} \tab positive Cauchy centered at 0, with UISD-adjusted scale
#' }
#'
#' The default PET prior distribution uses
#' \code{RoBMA.get_option("default_bias_PET.scale")}. The default PEESE prior
#' distribution uses \code{RoBMA.get_option("default_bias_PEESE.scale")} after
#' rescaling to the analyzed effect-size measure. The default weightfunction
#' prior distribution uses \code{RoBMA.get_option("default_bias_weightfunction.alpha")}
#' for each cumulative-Dirichlet alpha parameter.
#'
#' ## Model-averaged publication-bias priors
#'
#' \code{\link{RoBMA}} averages over publication-bias components. By default,
#' \code{prior_bias_null} is \code{\link{prior_none}()} and \code{model_type}
#' determines the alternative components:
#'
#' \tabular{ll}{
#' \code{"PSMA"} \tab six weight functions, PET, and PEESE \cr
#' \code{"6w"} \tab six weight functions \cr
#' \code{"2w"} \tab two two-sided weight functions \cr
#' \code{"PP"} \tab PET and PEESE
#' }
#'
#' Custom \code{prior_bias} replaces the preset alternative components. If
#' \code{prior_bias} is omitted, \code{model_type} still supplies the default
#' alternative components even when \code{prior_bias_null} is customized or
#' omitted. Setting \code{prior_bias = NULL} or \code{prior_bias = FALSE}
#' omits all alternative bias-adjustment models. Setting
#' \code{prior_bias_null = NULL} or \code{prior_bias_null = FALSE} omits the
#' no-bias component.
#'
#' Each publication-bias component has a prior model weight. User-specified
#' components receive equal weights unless their prior objects set
#' \code{prior_weights}. The \code{"2w"} and \code{"6w"} presets split the
#' alternative bias mass equally across their weight functions, \code{"PP"}
#' splits it equally across PET and PEESE, and \code{"PSMA"} assigns half of the
#' alternative bias mass to the six weight functions combined and one quarter
#' each to PET and PEESE.
#'
#' Publication-bias prior distributions are not rescaled by \code{rescale_priors}.
#' The exception is the default PEESE prior distribution, whose scale is adjusted
#' to the effect-size measure via the unit-information standard deviation.
#'
#' @seealso
#' \code{\link{prior_weightfunction}}, \code{\link{prior_PET}},
#' \code{\link{prior_PEESE}}, \code{\link{prior_none}},
#' \code{\link{RoBMA_prior_specification}},
#' \code{\link{prior_specification}}
#'
#' @aliases publication_bias_prior_specification bias_prior_specification
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
#' @param prior_effect prior distribution(s) for the alternative effect
#' component(s).
#' @param prior_heterogeneity prior distribution(s) for the alternative
#' heterogeneity component(s).
#' @param prior_mods prior distribution(s) for alternative moderator
#' components. A single prior applies to all terms; a named list can specify
#' term-specific components.
#' @param prior_scale prior distribution(s) for alternative scale-regression
#' components. A single prior applies to all terms; a named list can specify
#' term-specific components.
#' @param prior_heterogeneity_allocation prior distribution(s) for the
#' alternative cluster-level heterogeneity allocation component(s).
#' @param prior_bias prior distribution(s) for alternative publication-bias
#' component(s), such as weight functions, PET, or PEESE.
#'
#' Alternative prior arguments can be:
#' \itemize{
#'   \item A single prior distribution object (creates a mixture with one alternative)
#'   \item A list of prior distributions (creates a mixture with multiple alternatives)
#'   \item `NULL` or `FALSE` (omits the alternative hypothesis component)
#' }
#' See \code{\link{publication_bias_prior_specification}} for details on
#' specifying publication-bias priors and \code{\link{prior_specification}} for
#' details on specifying meta-analytic parameter priors.
#'
#' @param prior_effect_null prior distribution(s) for the null effect
#' component(s).
#' @param prior_heterogeneity_null prior distribution(s) for the null
#' heterogeneity component(s).
#' @param prior_mods_null prior distribution(s) for null moderator components.
#' A single prior applies to all terms; a named list can specify term-specific
#' components.
#' @param prior_scale_null prior distribution(s) for null scale-regression
#' components. A single prior applies to all terms; a named list can specify
#' term-specific components.
#' @param prior_heterogeneity_allocation_null prior distribution(s) for the
#' null cluster-level heterogeneity allocation component(s).
#' @param prior_bias_null prior distribution(s) for null publication-bias
#' component(s), usually `prior_none()`. See
#' \code{\link{publication_bias_prior_specification}}.
#'
#' Null prior arguments can be:
#' \itemize{
#'   \item A single prior distribution object (creates a mixture with one null)
#'   \item A list of prior distributions (creates a mixture with multiple nulls)
#'   \item `NULL` or `FALSE` (omits the null hypothesis component)
#' }
#' Defaults to a point mass (spike) at zero for effect and heterogeneity parameters.
#'
#' @param model_type character string specifying predefined publication-bias model
#' ensembles. One of:
#' \itemize{
#'   \item `"PSMA"` (default): Full RoBMA-PSMA ensemble with 6 weight functions + PET + PEESE
#'   \item `"6w"`: Six weight function models
#'   \item `"2w"`: Two weight function models 
#'   \item `"PP"`: PET-PEESE models only
#' }
#' Custom `prior_bias` replaces the preset alternative bias components. If
#' `prior_bias` is omitted, `model_type` determines the default alternatives even
#' when `prior_bias_null` is customized or omitted.
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
#' For top-level mixture components, setting a prior argument to `NULL` or `FALSE`
#' removes that hypothesis component:
#' \itemize{
#'   \item `prior_effect_null = NULL`: No null hypothesis for effect (assumes effect exists)
#'   \item `prior_effect = NULL`: No alternative hypothesis (assumes no effect)
#'   \item `prior_bias_null = NULL` or `FALSE`: No "no bias" model (assumes some bias mechanism)
#'   \item `prior_bias = NULL` or `FALSE`: No bias model (assumes no bias mechanism)
#' }
#' For moderator and scale regression terms, an omitted or whole-argument `NULL`
#' `prior_mods` / `prior_scale` means "use defaults". To omit a component for a
#' specific term, use a named list entry such as `prior_mods_null = list(x1 = NULL)`.
#'
#' ### Parameter-specific customization for moderators
#'
#' For moderator priors (`prior_mods`, `prior_mods_null`), you can specify different
#' priors for different predictors using named lists. A single prior object applies
#' to all moderator terms:
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
#' See \code{\link{publication_bias_prior_specification}} for the bias-prior
#' constructors and preset model spaces.
#'
#' ### Heterogeneity allocation (\eqn{\rho}) for multilevel models
#'
#' When `cluster` is specified, \eqn{\rho} controls the cluster-level variance
#' allocation. The decomposition is
#' \eqn{\tau_{between} = \tau\sqrt{\rho}} and
#' \eqn{\tau_{within} = \tau\sqrt{1 - \rho}}. By default, only an alternative
#' (Beta(1,1)) is specified, but null hypotheses can be added:
#' \preformatted{
#' RoBMA(...,
#'   cluster = study,
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
#' probability). User-specified components receive equal weights unless their
#' prior objects set `prior_weights`. Publication-bias presets use fixed weights:
#' `"2w"` and `"6w"` split the alternative bias mass equally across their
#' weight functions, `"PP"` splits it equally across PET and PEESE, and `"PSMA"`
#' assigns half of the alternative bias mass to the six weight functions
#' combined and one quarter each to PET and PEESE. Custom weights can be
#' specified via the `prior_weights` argument in individual prior objects.
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
#' \code{\link{publication_bias_prior_specification}} for publication-bias priors,
#' \code{\link{RoBMA}} for the main model-averaging function,
#' \code{\link[BayesTools]{prior}} for creating prior distribution objects
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'
#'   # Default mixture priors (spike-and-slab for effect and heterogeneity)
#'   fit1 <- RoBMA(
#'     yi      = yi,
#'     vi      = vi,
#'     data    = dat.lehmann2018,
#'     measure = "SMD",
#'     seed    = 1,
#'     silent  = TRUE
#'   )
#'
#'   # Custom alternative priors (two effect size scales)
#'   fit2 <- RoBMA(
#'     yi           = yi,
#'     vi           = vi,
#'     data         = dat.lehmann2018,
#'     measure      = "SMD",
#'     prior_effect = list(
#'       prior("normal", list(mean = 0, sd = 0.35)),
#'       prior("normal", list(mean = 0, sd = 0.70))
#'     ),
#'     seed         = 1,
#'     silent       = TRUE
#'   )
#'
#'   # No null hypothesis for effect (assume effect exists)
#'   fit3 <- RoBMA(
#'     yi                = yi,
#'     vi                = vi,
#'     data              = dat.lehmann2018,
#'     measure           = "SMD",
#'     prior_effect_null = NULL,
#'     seed              = 1,
#'     silent            = TRUE
#'   )
#'
#'   # Only selection model bias adjustment (no PET-PEESE)
#'   fit4 <- RoBMA(
#'     yi         = yi,
#'     vi         = vi,
#'     data       = dat.lehmann2018,
#'     measure    = "SMD",
#'     model_type = "6w",
#'     seed       = 1,
#'     silent     = TRUE
#'   )
#'
#'   # Meta-regression with different null specifications per moderator
#'   fit5 <- RoBMA(
#'     yi              = yi,
#'     vi              = vi,
#'     mods            = ~ Preregistered,
#'     data            = dat.lehmann2018,
#'     measure         = "SMD",
#'     prior_mods_null = list(Preregistered = NULL),
#'     seed            = 1,
#'     silent          = TRUE
#'   )
#' }
#' }
#'
#' @aliases RoBMA_prior_specification
NULL




# stop if public fitting constructors did not receive an explicit measure
.stop_missing_measure <- function(caller) {

  stop(
    paste0(
      caller, " requires explicit 'measure'. Use 'GEN' only for general ",
      "effect sizes; otherwise specify one of 'SMD', 'ZCOR', 'RR', 'OR', ",
      "'HR', 'RD', or 'IRR'."
    ),
    call. = FALSE
  )
}

# verify that the the measure is available
.check_measure <- function(measure, class = "norm") {

  BayesTools::check_char(measure, "measure")

  if (identical(class, "glmm")) {
    if (!measure %in% c("OR", "IRR")) {
      stop(
        "GLMM models support only 'OR' and 'IRR' measures.",
        call. = FALSE
      )
    }
    return()
  }

  BayesTools::check_char(
    measure,
    "measure",
    allow_values = c("SMD", "ZCOR", "RR", "OR", "HR", "RD", "IRR", "GEN")
  )

  return()
}


### prior scale functions
.get_unit_information_sd <- function(data, measure) {

  if (measure %in% c("SMD", "ZCOR", "RR", "OR", "HR", "IRR")) {
    # use the pre-specified UISD for known effect sizes
    prior_unit_information_sd <- .get_unit_information_sd.known(measure)
  } else {
    # compute sd based on the sample sizes for the general effect sizes
    ni <- data[["outcome"]][["ni"]]
    if (all(is.na(ni))) {
      stop("Sample size 'ni' or unit information sd 'unit_information_sd' must be specified to set-up prior distributions without known UISD.", call. = FALSE)
    }
    if (any(is.na(ni))) {
      stop("Sample size 'ni' must be complete to estimate unit information sd for measures without known UISD.", call. = FALSE)
    }
    prior_unit_information_sd <- .get_unit_information_sd.norm(sei = data[["outcome"]][["sei"]], ni = ni)

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
#' @param sei a complete numeric vector of strictly positive standard errors
#' for each study.
#' @param ni a complete numeric vector of strictly positive sample sizes for
#' each study, with the same length as `sei`.
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
#' [brma()] or related functions. This ensures consistent prior scaling
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
    prior <- switch(
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
      if (is.null(type)) {
        stop("Informed medicine priors are only available for measures 'SMD', 'OR', 'RR', 'RD', and 'HR'.", call. = FALSE)
      }

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
  if (missing(prior) || is.null(prior)) {
    prior <- BayesTools::prior("beta", parameters = list("alpha" = 1, "beta" = 1))
    return(prior)
  }

  # check the user specified prior distribution
  prior <- .check_prior.restricted_01(prior, prior_name = "heterogeneity_allocation")

  return(prior)
}
.assign_prior.baserate                 <- function(prior) {

  # use default settings if prior distribution is not specified
  if (missing(prior) || is.null(prior)) {
    prior <- BayesTools::prior_factor("beta", parameters = list("alpha" = 1, "beta" = 1), contrast = "independent")
    return(prior)
  }

  # check the user specified prior distribution
  prior <- .check_prior.restricted_01(prior, prior_name = "baserate")

  # transform the prior into independent factor prior
  prior <- BayesTools::prior_factor(prior[["distribution"]], parameters = prior[["parameters"]], truncation = prior[["truncation"]], contrast = "independent")

  return(prior)
}
.assign_prior.lograte                  <- function(prior, data) {
  # the log-rate is unbounded, so use a proper data-based default prior
  # rather than an improper flat nuisance-parameter prior

  if (missing(prior) || is.null(prior)) {
    return(.get_unit_information_lograte_prior(
      x1i = data[["outcome"]][["x1i"]],
      x2i = data[["outcome"]][["x2i"]],
      t1i = data[["outcome"]][["t1i"]],
      t2i = data[["outcome"]][["t2i"]])
    )
  }

  # check the user specified prior distribution
  prior <- .check_prior.simple_or_point(prior, prior_name = "prior_lograte")

  # transform the prior into independent factor prior
  prior <- BayesTools::prior_factor(prior[["distribution"]], parameters = prior[["parameters"]], truncation = prior[["truncation"]], contrast = "independent")

  return(prior)
}
.get_PEESE_UISD_ratio                  <- function(measure, data, prior_unit_information_sd) {

  UISD_SMD <- .get_unit_information_sd.known("SMD")
  if (missing(prior_unit_information_sd)) {
    UISD_measure <- .get_unit_information_sd(data = data, measure = measure)
  } else {
    UISD_measure <- prior_unit_information_sd
  }

  return(UISD_SMD / UISD_measure)
}
.assign_prior.bias                     <- function(prior, measure, data, prior_unit_information_sd, bias_type, steps) {

  # assigns either a weight function or PET/PEESE model based on bias_type
  # there is no scaling of the publication bias priors based on rescale_priors
  # (only PEESE prior is rescaled to the match the UISD of the measure)

  if (bias_type == "selmodel") {

    # specify default settings / check the prior distribution
    if (missing(prior) || is.null(prior)) {
      if (missing(steps)) {
        steps <- c(0.025)
      }
      prior <- BayesTools::prior_weightfunction(
        "one-sided",
        steps   = steps,
        weights = BayesTools::wf_cumulative(
          alpha = rep(RoBMA.get_option("default_bias_weightfunction.alpha"), length(steps) + 1)
        )
      )
    } else if (!.prior_is_selection_kernel(prior)) {
      stop(
        "'prior_bias' must be a `prior_weightfunction` or `prior_bias` object",
        call. = FALSE
      )
    } else if (.prior_has_phacking(prior)) {
      .selection_stop_phacking_deferred()
    }

    return(prior)
  }

  if (bias_type == "PET") {

    # specify default settings / check the prior distribution
    if (missing(prior) || is.null(prior)) {
      # PET prior is invariant to different effect size types
      prior <- BayesTools::prior_PET("cauchy", parameters = list("location" = 0, "scale" = RoBMA.get_option("default_bias_PET.scale")), truncation = list(0, Inf))
    } else if (!BayesTools::is.prior.PET(prior)) {
      stop("'prior_bias' must be a `prior_PET` object", call. = FALSE)
    }

    return(prior)
  }

  if (bias_type == "PEESE") {

    # specify default settings / check the prior distribution
    if (missing(prior) || is.null(prior)) {
      ### PEESE prior is inversely related to the specified effect size
      # for SMD we use Cauchy with scale 1
      # for other effect sizes, we re-scale by the corresponding UISD
      UISD_ratio <- .get_PEESE_UISD_ratio(
        measure = measure, data = data, prior_unit_information_sd = prior_unit_information_sd)
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
.term_default_continuous_prior         <- function(prior) {

  if (!BayesTools::is.prior.factor(prior)) {
    return(prior)
  }

  distribution <- prior[["distribution"]]
  if (distribution %in% c("mnormal", "mcauchy", "mt")) {
    distribution <- substring(distribution, 2)
  }

  return(BayesTools::prior(
    distribution = distribution,
    parameters   = prior[["parameters"]],
    truncation   = prior[["truncation"]]
  ))
}
.term_default_factor_prior             <- function(prior, contrast) {

  if (BayesTools::is.prior.factor(prior)) {
    return(prior)
  }

  distribution <- prior[["distribution"]]
  if (contrast != "treatment" && !distribution %in% c("spike", "mnormal", "mcauchy", "mt")) {
    distribution <- paste0("m", distribution)
  }

  return(BayesTools::prior_factor(
    distribution = distribution,
    parameters   = prior[["parameters"]],
    truncation   = prior[["truncation"]],
    contrast     = contrast
  ))
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
  if (parameter == "mods" &&
      attr(stats::terms(attr(data[[parameter]], "formula")), "intercept") == 0) {
    prior_intercept <- BayesTools::prior("spike", parameters = list(0))
  }

  ### check the user-specified priors
  user_default_prior <- NULL
  if (!missing(prior_list) && BayesTools::is.prior(prior_list)) {
    .check_prior.term(prior_list)
    user_default_prior <- prior_list
    prior_list         <- list()
  } else if (!missing(prior_list) && is.null(prior_list)) {
    prior_list <- list()
  } else if (!missing(prior_list)) {
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
      if (is.null(type)) {
        stop("Informed medicine priors are only available for measures 'SMD', 'OR', 'RR', 'RD', and 'HR'.", call. = FALSE)
      }

      if (parameter == "mods") {
        # treat moderators by scaling them according to 1/2 of the effect
        default_prior_continuous <- BayesTools::prior_informed(prior_name, parameter = "effect", type = type)
        default_prior_continuous <- .rescale_prior_object(default_prior_continuous, RoBMA.get_option("default_informed_priors.mods"))
      } else if (parameter == "scale") {
        # scale parameters use the default scaling (no informed priors exist)
        default_prior_continuous <- BayesTools::prior("normal", parameters = list("mean" = 0, "sd" = RoBMA.get_option("default_informed_priors.scale")))
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

  if (!is.null(user_default_prior)) {
    default_prior_continuous <- .term_default_continuous_prior(user_default_prior)
    default_prior_factor     <- .term_default_factor_prior(
      prior    = user_default_prior,
      contrast = attr(data, "set_contrast_factor_predictors")
    )
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

# Remove invalid fixed-effect heterogeneity components for scale regression.
.remove_zero_scale_intercept_prior <- function(prior_intercept) {

  is_zero_point <- function(prior) {
    BayesTools::is.prior.point(prior) && mean(prior) == 0
  }

  if (!BayesTools::is.prior.mixture(prior_intercept)) {
    if (is_zero_point(prior_intercept)) {
      stop(
        "At least one non-zero between-study heterogeneity prior must be defined for scale regression.",
        call. = FALSE
      )
    }
    return(prior_intercept)
  }

  keep <- !vapply(prior_intercept, is_zero_point, logical(1))

  if (all(keep))
    return(prior_intercept)

  warning(
    "Spike(0) prior distribution for between-study heterogeneity tau was removed since scale regression requires non-zero between-study heterogeneity.",
    immediate. = TRUE
  )

  if (!any(keep)) {
    stop(
      "At least one non-zero between-study heterogeneity prior must be defined for scale regression.",
      call. = FALSE
    )
  }

  components <- attr(prior_intercept, "components")
  prior_list <- prior_intercept[keep]
  attributes(prior_list) <- NULL

  return(BayesTools::prior_mixture(
    prior_list = prior_list,
    components = components[keep]
  ))
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

  if (length(priors_null) + length(priors_alt) == 0) {
    stop(gettextf("At least one prior distribution needs to be defined for '%1$s'.", parameter), call. = FALSE)
  }

  ### specify the corresponding mixture prior
  prior <- BayesTools::prior_mixture(
    prior_list = c(priors_null, priors_alt),
    is_null    = c(rep(TRUE, length(priors_null)), rep(FALSE, length(priors_alt)))
  )

  return(prior)
}
.default_prior.bias_alt                         <- function(model_type, measure, data, prior_unit_information_sd) {

  if (model_type == "2w") {
    return(list(
      BayesTools::prior_weightfunction("two-sided", c(0.05),       BayesTools::wf_cumulative(c(1, 1)),    prior_weights = 1/2),
      BayesTools::prior_weightfunction("two-sided", c(0.05, 0.10), BayesTools::wf_cumulative(c(1, 1, 1)), prior_weights = 1/2)
    ))
  }

  if (model_type == "6w") {
    return(list(
      BayesTools::prior_weightfunction("two-sided", c(0.05),             BayesTools::wf_cumulative(c(1, 1)),       prior_weights = 1/6),
      BayesTools::prior_weightfunction("two-sided", c(0.05, 0.10),       BayesTools::wf_cumulative(c(1, 1, 1)),    prior_weights = 1/6),
      BayesTools::prior_weightfunction("one-sided", c(0.05),             BayesTools::wf_cumulative(c(1, 1)),       prior_weights = 1/6),
      BayesTools::prior_weightfunction("one-sided", c(0.025, 0.05),      BayesTools::wf_cumulative(c(1, 1, 1)),    prior_weights = 1/6),
      BayesTools::prior_weightfunction("one-sided", c(0.05, 0.5),        BayesTools::wf_cumulative(c(1, 1, 1)),    prior_weights = 1/6),
      BayesTools::prior_weightfunction("one-sided", c(0.025, 0.05, 0.5), BayesTools::wf_cumulative(c(1, 1, 1, 1)), prior_weights = 1/6)
    ))
  }

  UISD_ratio <- .get_PEESE_UISD_ratio(
    measure = measure, data = data, prior_unit_information_sd = prior_unit_information_sd)

  if (model_type == "PP") {
    return(list(
      BayesTools::prior_PET(distribution = "Cauchy",   parameters = list(0, RoBMA.get_option("default_bias_PET.scale")),                truncation = list(0, Inf), prior_weights = 1/2),
      BayesTools::prior_PEESE(distribution = "Cauchy", parameters = list(0, RoBMA.get_option("default_bias_PEESE.scale") * UISD_ratio), truncation = list(0, Inf), prior_weights = 1/2)
    ))
  }

  return(list(
    BayesTools::prior_weightfunction("two-sided", c(0.05),             BayesTools::wf_cumulative(c(1, 1)),       prior_weights = 1/12),
    BayesTools::prior_weightfunction("two-sided", c(0.05, 0.10),       BayesTools::wf_cumulative(c(1, 1, 1)),    prior_weights = 1/12),
    BayesTools::prior_weightfunction("one-sided", c(0.05),             BayesTools::wf_cumulative(c(1, 1)),       prior_weights = 1/12),
    BayesTools::prior_weightfunction("one-sided", c(0.025, 0.05),      BayesTools::wf_cumulative(c(1, 1, 1)),    prior_weights = 1/12),
    BayesTools::prior_weightfunction("one-sided", c(0.05, 0.5),        BayesTools::wf_cumulative(c(1, 1, 1)),    prior_weights = 1/12),
    BayesTools::prior_weightfunction("one-sided", c(0.025, 0.05, 0.5), BayesTools::wf_cumulative(c(1, 1, 1, 1)), prior_weights = 1/12),
    BayesTools::prior_PET(distribution = "Cauchy",   parameters = list(0, RoBMA.get_option("default_bias_PET.scale")),                truncation = list(0, Inf), prior_weights = 1/4),
    BayesTools::prior_PEESE(distribution = "Cauchy", parameters = list(0, RoBMA.get_option("default_bias_PEESE.scale") * UISD_ratio), truncation = list(0, Inf), prior_weights = 1/4)
  ))
}
.assign_prior.bias_mixture                     <- function(
    prior, prior_null, measure, data, prior_unit_information_sd, model_type) {

  ### use precanned publication bias adjustment models
  # if neither `prior` or `prior_null` input is specified use the `model_type`
  # there is no scaling of the publication bias priors based on rescale_priors
  # (only PEESE prior is rescaled to the match the UISD of the measure)

  if (missing(prior) && missing(prior_null)) {
    priors_null <- list(BayesTools::prior_none(prior_weights = 1))
    priors_alt  <- .default_prior.bias_alt(
      model_type = model_type, measure = measure, data = data,
      prior_unit_information_sd = prior_unit_information_sd)
    prior <- BayesTools::prior_mixture(
      prior_list = c(priors_null, priors_alt),
      is_null    = c(TRUE, rep(FALSE, length(priors_alt)))
    )

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
    priors_alt <- .default_prior.bias_alt(
      model_type = model_type, measure = measure, data = data,
      prior_unit_information_sd = prior_unit_information_sd)
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

  if (length(priors_null) + length(priors_alt) == 0) {
    stop("At least one prior distribution needs to be defined for 'bias'.", call. = FALSE)
  }

  if (length(priors_alt) == 0 &&
      length(priors_null) == 1 &&
      BayesTools::is.prior.none(priors_null[[1]])) {
    return(priors_null[[1]])
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

  if (length(priors_null) + length(priors_alt) == 0) {
    stop("At least one prior distribution needs to be defined for 'heterogeneity_allocation'.", call. = FALSE)
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
  if (!missing(prior_list) && !BayesTools::is.prior(prior_list) && length(prior_list) > 0 && is.list(prior_list)) {
    par_alt_omit <- names(prior_list)[sapply(prior_list, function(x) is.null(x) || isFALSE(x))]
    # filter out NULL/FALSE values before passing to .assign_prior_list.terms
    prior_list <- prior_list[!sapply(prior_list, function(x) is.null(x) || isFALSE(x))]
  } else {
    par_alt_omit <- NULL
  }
  if (!missing(prior_list_null) && !BayesTools::is.prior(prior_list_null) && length(prior_list_null) > 0 && is.list(prior_list_null)) {
    par_null_omit <- names(prior_list_null)[sapply(prior_list_null, function(x) is.null(x) || isFALSE(x))]
    # filter out NULL/FALSE values before passing to .assign_prior_list.terms
    prior_list_null <- prior_list_null[!sapply(prior_list_null, function(x) is.null(x) || isFALSE(x))]
  } else {
    par_null_omit <- NULL
  }


  if (parameter == "scale") {
    # in the case of scale moderation, spike at zero prior (i.e., fixed-effect model) is not possible since the intercept
    # is on a log scale
    # remove spike(0) priors and rebuild mixture metadata
    prior_intercept <- .remove_zero_scale_intercept_prior(prior_intercept)
  }
  if (parameter == "mods" &&
      attr(stats::terms(attr(data[[parameter]], "formula")), "intercept") == 0) {
    prior_intercept <- BayesTools::prior("spike", parameters = list(0))
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
    rescale_priors = rescale_priors, set_default_spike = TRUE)

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

  prior_list <- .repair_formula_prior_list(
    prior_list = prior_list,
    parameter  = switch(parameter, "mods" = "mu", "scale" = "log_tau")
  )

  return(prior_list)
}

# Restore BayesTools formula metadata lost when formula prior mixtures are rebuilt.
.repair_formula_prior_list <- function(prior_list, parameter) {

  if (is.null(prior_list)) {
    return(prior_list)
  }

  for (i in seq_along(prior_list)) {
    prior_list[[i]] <- .repair_formula_prior(
      prior    = prior_list[[i]],
      parameter = parameter
    )
  }

  return(prior_list)
}
.repair_formula_prior <- function(prior, parameter) {

  attr(prior, "parameter") <- parameter

  formula_attributes <- c(
    "interaction", "interaction_terms", "levels", "level_names",
    "term_components", "factor_terms", "factor_contrasts",
    "factor_design", "factor_cell_names", "multiply_by"
  )

  if (.has_formula_prior_components(prior)) {
    source_attributes <- NULL
    for (i in seq_along(prior)) {
      if (!is.null(attr(prior[[i]], "parameter"))) {
        source_attributes <- attributes(prior[[i]])
        break
      }
    }

    if (!is.null(source_attributes)) {
      for (attribute in formula_attributes) {
        if (is.null(attr(prior, attribute)) && !is.null(source_attributes[[attribute]])) {
          attr(prior, attribute) <- source_attributes[[attribute]]
        }
      }
    }

    for (i in seq_along(prior)) {
      attr(prior[[i]], "parameter") <- parameter
      for (attribute in formula_attributes) {
        if (is.null(attr(prior[[i]], attribute)) && !is.null(attr(prior, attribute))) {
          attr(prior[[i]], attribute) <- attr(prior, attribute)
        }
      }
    }
  }

  return(prior)
}
.has_formula_prior_components <- function(prior) {

  return(BayesTools::is.prior.mixture(prior) ||
    BayesTools::is.prior.spike_and_slab(prior))
}
### check specified prior distributions whether they follow basic requirements
.check_prior.effect                   <- function(prior) {

  prior <- .check_prior.simple_or_point(prior, prior_name = "prior_effect")

  return(prior)
}
.check_prior.heterogeneity            <- function(prior) {

  prior <- .check_prior.simple_or_point(prior, prior_name = "prior_heterogeneity")

  # check range restriction
  if (BayesTools::is.prior.point(prior) && mean(prior) < 0) {
    stop("Point prior distribution for 'prior_heterogeneity' must be at least 0.", call. = FALSE)
  }
  if (prior[["truncation"]][["lower"]] < 0) {
    warning("Lower truncation point for 'prior_heterogeneity' prior distribution must be at least 0. The prior distribution was modified.",
            immediate. = TRUE, call. = FALSE)
    prior[["truncation"]][["lower"]] <- 0
  }

  return(prior)
}
.check_prior.simple_or_point          <- function(prior, prior_name) {

  # check object type
  if (!BayesTools::is.prior(prior))
    stop(sprintf("The '%s' argument must be a prior distribution.", prior_name), call. = FALSE)
  if (!(BayesTools::is.prior.simple(prior) || BayesTools::is.prior.point(prior)))
    stop(sprintf("The '%s' argument must be a either a simple or point prior distribution.", prior_name), call. = FALSE)

  return(prior)
}
.check_prior.restricted_01            <- function(prior, prior_name) {

  prior <- .check_prior.simple_or_point(prior, prior_name)

  # check range restriction
  if (BayesTools::is.prior.point(prior) && (mean(prior) < 0 || mean(prior) > 1)) {
    stop(sprintf("Point prior distribution for '%s' must be within [0, 1].", prior_name), call. = FALSE)
  }
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

  return(prior)
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
    bias_type <- "none"
  } else if (.prior_is_selection_kernel(prior)) {
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

  if (anyNA(x1i) || anyNA(x2i) || anyNA(t1i) || anyNA(t2i)) {
    stop("Default 'prior_lograte' requires complete Poisson event counts and person-times.", call. = FALSE)
  }

  total_events <- sum(x1i + x2i)
  total_time   <- sum(t1i + t2i)
  if (!is.finite(total_events) || total_events <= 0) {
    stop("Default 'prior_lograte' requires at least one observed Poisson event. Specify 'prior_lograte' manually for all-zero event counts.", call. = FALSE)
  }
  if (!is.finite(total_time) || total_time <= 0) {
    stop("Default 'prior_lograte' requires positive finite total person-time.", call. = FALSE)
  }

  # compute the unit information standard deviation from Poisson data
  sigma_prior <- .get_unit_information_sd.pois(x1i, x2i, t1i, t2i)

  # compute the pooled crude log-rate for the mean
  mu_prior     <- log(total_events / total_time)
  if (!is.finite(mu_prior) || !is.finite(sigma_prior)) {
    stop("Default 'prior_lograte' could not be computed from the supplied Poisson data.", call. = FALSE)
  }

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
  if (BayesTools::is.prior.point(prior) || BayesTools::is.prior.none(prior))
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
    # gamma is an auxiliary parameter for non-centered random-effects parameterization (cluster-level)
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
  if (!missing(model_type))
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
  # only for multilevel models
  if (.is_data_multilevel(data)) {
    prior_outcome[["rho"]]   <- .assign_prior.heterogeneity_allocation_mixture(
      prior = prior_heterogeneity_allocation, prior_null = prior_heterogeneity_allocation_null)
    # gamma is an auxiliary parameter for non-centered random-effects parameterization (cluster-level)
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
  # only for RoBMA
  if (!missing(model_type)) {
    prior_outcome[["bias"]] <- .assign_prior.bias_mixture(
      prior = prior_bias, prior_null = prior_bias_null, measure = measure,
      data = data, prior_unit_information_sd = prior_unit_information_sd, model_type = model_type
    )
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
