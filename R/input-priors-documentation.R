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
#' functions, explicit `NULL` or `FALSE` sets a spike at zero. For
#' \code{\link{brma.mv}} random-effect formulas, this can also be a partial or
#' complete \code{\link{prior_random}} object; see
#' \code{\link{random_effect_prior_specification}}.
#' @param prior_mods prior distribution for the moderators (\eqn{\beta}) parameters.
#' A single prior applies to all terms; a named list can specify term-specific
#' priors. If omitted or `NULL`, default priors are used. This argument can be
#' specified only when `mods` is specified.
#' @param prior_scale prior distribution for the scale (\eqn{\delta}) parameters.
#' A single prior applies to all terms; a named list can specify term-specific
#' priors. If omitted or `NULL`, default priors are used. This argument can be
#' specified only when `scale` is specified.
#' @param prior_heterogeneity_allocation prior distribution for the fraction of
#' heterogeneity allocated to the cluster-level component in multilevel models
#' (\eqn{\rho}). If omitted or `NULL`, defaults to `Beta(1, 1)`. This argument
#' can be specified only when `cluster` is specified.
#' @param prior_baserate prior distribution for the estimate-specific midpoint
#' base-rate probability in binomial GLMM models (`measure = "OR"`). If omitted
#' or `NULL`, defaults to independent `Beta(1, 1)` priors. Only beta and point
#' priors are supported; beta priors may be truncated within `[0, 1]`, and point
#' priors must lie strictly inside `(0, 1)`.
#' @param prior_lograte prior distribution for the estimate-specific midpoint
#' log-rate in Poisson GLMM models (`measure = "IRR"`). If omitted or `NULL`, a
#' normal prior centered at the pooled crude log-rate is used independently for
#' each estimate. Its standard deviation defaults to `1` and can be changed via
#' `RoBMA.options(default_lograte.sd = ...)`. Rescaling all person-times shifts
#' the prior mean by the corresponding log conversion while leaving its
#' standard deviation unchanged. Only normal and point priors are supported;
#' normal priors may be truncated.
#' @param prior_unit_information_sd numeric. The unit information standard deviation (\eqn{\sigma_{unit}}).
#' Cannot be used together with `prior_informed_field`.
#' @param rescale_priors numeric. A scaling factor for supported prior distributions.
#' Point and none priors are unchanged. For constructors with publication-bias
#' prior distributions, `rescale_priors` does not rescale them except for the
#' default PEESE prior's UISD adjustment. Defaults to 1.
#' @param standardize_continuous_predictors logical. Whether to standardize continuous predictors.
#' Defaults to `TRUE`.
#' @param set_contrast_factor_predictors character. How to set contrast for factor predictors.
#' Defaults are constructor-specific and shown in each function usage. Single-model
#' constructors use `"treatment"`; model-averaging constructors use `"meandif"`.
#' Accepted values are `"treatment"`, `"meandif"`, `"orthonormal"`, and
#' `"independent"`. Independent coding creates one coefficient per factor
#' level. Random-coefficient blocks reuse the resolved fixed-factor contrast
#' unless `random_block(contrasts = ...)` overrides it. A warning is issued
#' when default or inherited independent contrasts make a fixed or random
#' design overspecified.
#' @param prior_informed_field character. The field of the informed prior distributions.
#' Omit to use the standard default prior specification; explicit `NULL` is invalid.
#' @param prior_informed_subfield character. The subfield of the informed prior distributions.
#' Can only be specified together with `prior_informed_field`. Omit to use the
#' field-specific default, such as `"Cochrane"` for `prior_informed_field =
#' "medicine"`; explicit `NULL` is invalid.
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
#' For `brma.mv()` random-effect formulas, see
#' \code{\link{random_effect_prior_specification}}.
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
#' \code{\link{random_effect_prior_specification}},
#' \code{\link[BayesTools]{prior}}, \code{\link{RoBMA.options}},
#' \code{\link{brma}}
#' @aliases prior_specification
NULL


#' @title Random-effect prior specification
#' @name random_effect_prior_specification
#'
#' @description
#' Random-effect prior distributions are used by \code{\link{brma.mv}} when a
#' `random` formula is supplied. The formula and the prior have deliberately
#' separate roles. The formula defines the random-effect design, grouping
#' variables, and covariance family. The `prior_heterogeneity` argument defines
#' the prior scale and can additionally customize correlations, variance
#' allocation, monitoring, prediction, and computational parameterization.
#'
#' RoBMA constructs a complete random-effect prior when
#' `prior_heterogeneity` is omitted. An ordinary positive prior distribution
#' changes only the base random-effect standard-deviation prior. A
#' \code{\link{prior_random}} object can either customize selected non-scale
#' settings while inheriting RoBMA's scale defaults, or explicitly own the
#' complete scale architecture.
#'
#' @param random an optional formula or list of formulas specifying
#' BayesTools random-effect terms for \code{\link{brma.mv}}. See
#' \code{\link{random_effect_formula_tags}} for the formula syntax, component
#' naming rules, and public semantic parameter names.
#' @param prior_heterogeneity prior distribution for residual heterogeneity.
#' In random-formula \code{\link{brma.mv}} models, omission constructs all
#' defaults; an ordinary positive prior distribution replaces the base or
#' aggregate random-effect SD prior while retaining automatic allocation and correlation
#' defaults; and a \code{\link{prior_random}} object supplies partial or
#' complete random-effect settings as described below.
#' @param scale an optional formula or named list of formulas specifying
#' location-scale predictors. With random-effect formulas, `scale` models the
#' row-wise total random-effect standard deviation consumed by a random-effect
#' component.
#' @param prior_scale prior distribution for scale-regression coefficients.
#' A single prior applies to all terms; a named list can specify term-specific
#' priors. If omitted or `NULL`, default priors are used.
#' @param prior_unit_information_sd numeric. The unit information standard
#' deviation used to construct the default base heterogeneity prior distribution.
#' @param rescale_priors numeric. A scaling factor for supported prior
#' distributions. It rescales the base heterogeneity prior distribution before
#' RoBMA constructs the random-effect prior specification.
#'
#' @details
#' ## Division of responsibility
#'
#' The random-effect setup has three layers:
#' \tabular{lll}{
#' Input \tab Controls \tab Does not redefine \cr
#' `random` \tab blocks, grouping, coefficient/index basis, covariance family \tab prior scale or prior family \cr
#' ordinary `prior_heterogeneity` \tab base or aggregate random-effect SD prior \tab formula structure or automatic allocation \cr
#' `prior_random(...)` \tab block overrides, covariance priors, allocation, and policies \tab the structure parsed from `random`
#' }
#'
#' Plain `(expr | group)` syntax and `us()` / `un()` define unstructured
#' random-coefficient blocks. `||` and `diag()` define independent
#' random-coefficient blocks, while `id()` gives independent coefficients one
#' shared SD. In these families, `1`, `0`, and `-1` control the random intercept;
#' continuous slopes, factor slopes, and interactions are supported. Factor
#' coding comes from the resolved contrast metadata, not from intercept syntax.
#' Use `random_block(contrasts = ...)` when a block-specific factor basis is
#' required.
#'
#' The `cs()`, `hcs()`, `ar1()` / `ar()`, `har()`, and `car()` families instead
#' own a structure-specific index basis. Their left side supplies index
#' variables rather than a coefficient formula, so these structures reject
#' intercept controls and block contrast overrides. See
#' \code{\link{random_effect_formula_tags}} for the complete formula and
#' variable-type contract.
#'
#' Nested grouping syntax, such as `random = ~ 1 | study / esid`, creates
#' multiple blocks. A formula list creates top-level components:
#'
#' \preformatted{
#' random = list(
#'   study   = ~ 1 | study,
#'   outcome = ~ 1 | outcome
#' )
#' }
#' A bare formula or unnamed one-entry list has no redundant public component
#' prefix. An explicitly named one-entry list retains its name. In lists with
#' two or more entries, missing names become `component 1`, `component 2`, and
#' so on. Component names are used by summaries and by component-specific
#' `scale` models. Block overrides in `prior_random()` must match the resolved
#' block names, which can be inspected with
#' `print_prior(object, parameter = "random")`.
#'
#' ## Default standard-deviation prior
#'
#' RoBMA first constructs the same base heterogeneity SD prior described in
#' \code{\link{prior_specification}}. Under the standard UISD defaults,
#' \deqn{\sigma_{\mathrm{base}} \sim
#'       \mathrm{Normal}^{+}(0, \mathrm{UISD}/4).}
#' The scale multiplier can be changed through
#' `RoBMA.options(default_UISD.heterogeneity = ...)`. The same `measure`, `ni`,
#' `prior_unit_information_sd`, informed-prior, and `rescale_priors` rules used
#' for ordinary heterogeneity apply before the prior is assigned to the
#' random-effect structure.
#'
#' Supplying an ordinary positive prior through `prior_heterogeneity` replaces
#' this base prior only. RoBMA still determines whether it applies directly to
#' one SD, to `sd_total` for a total-variance allocation, or to `sd_common` for
#' a mean-variance allocation.
#'
#' ## Structure-specific defaults
#'
#' Let \eqn{K} denote the resolved number of random-coefficient columns for a
#' coefficient structure or the number of index levels for an index structure.
#' The formula family determines the number of marginal SDs and the default
#' correlation prior:
#' \tabular{llll}{
#' Structure \tab Basis \tab Marginal SDs \tab Default correlation prior \cr
#' `id()` \tab coefficients \tab one shared SD \tab none \cr
#' `diag()` or `||` \tab coefficients \tab one per column \tab none \cr
#' plain `|`, `us()` / `un()` \tab coefficients \tab one per column \tab `LKJ(1)` matrix \cr
#' `cs()` \tab index levels \tab one shared SD \tab `Uniform(-1 / (K - 1), 1)` on raw `cor` \cr
#' `hcs()` \tab index levels \tab one per level \tab `Uniform(-1 / (K - 1), 1)` on raw `cor` \cr
#' `ar1()` / `ar()` \tab ordered index \tab one shared SD \tab `Uniform(-1, 1)` on raw `cor` \cr
#' `har()` \tab ordered index \tab one per level \tab `Uniform(-1, 1)` on raw `cor` \cr
#' `car()` \tab numeric coordinates \tab one shared SD \tab `Uniform(0, 1)` on raw `cor`
#' }
#'
#' Correlation priors are added after the formula basis and \eqn{K} are
#' resolved, and only when the block has more than one random coefficient. The
#' CS/HCS lower bound is the exact positive-definiteness bound, so the default
#' includes all admissible negative correlations. `LKJ(1)` is uniform over
#' correlation matrices and does not favor positive correlations. In
#' particular, `hcs(index | group)` has one common pairwise correlation, while
#' `us(0 + index | group)` can estimate a different correlation for every pair
#' of resolved coefficient columns.
#'
#' ## Automatic variance allocation
#'
#' RoBMA uses symmetric `Dirichlet(1, ..., 1)` priors whenever a variance scale
#' must be divided. The allocation is hierarchical when the formula contains
#' multiple top-level components, nested blocks, or heterogeneous blocks:
#' \itemize{
#'   \item A single homogeneous block receives the base SD prior directly.
#'   \item A split across blocks or components uses
#'   `scale = "total_variance"`. For child weight \eqn{w_j},
#'   \deqn{\sigma_j = \sigma_{\mathrm{total}}\sqrt{w_j}, \qquad
#'         \sum_j \sigma_j^2 = \sigma_{\mathrm{total}}^2.}
#'   The public aggregate quantities are `sd_total` and `var_total`; components
#'   are `var_prop(...)`.
#'   \item A heterogeneous block, such as DIAG, US, HCS, or HAR with
#'   \eqn{K > 1}, adds an SD-component allocation with
#'   `scale = "mean_variance"`:
#'   \deqn{\sigma_j = \sigma_{\mathrm{common}}\sqrt{K w_j}, \qquad
#'         K^{-1}\sum_j \sigma_j^2 = \sigma_{\mathrm{common}}^2.}
#'   The public aggregate quantities are `sd_common` and `var_common`;
#'   components are `var_ratio(...)` and `sd_ratio(...)`.
#' }
#'
#' Consequently, the base prior always controls a clearly defined total or
#' average variance scale; it is not independently copied to every SD when that
#' would make the implied total variance grow with the number of components.
#'
#' ## Public random-effect parameter names
#'
#' BayesTools stores concrete coordinates, semantic quantities, and aliases in
#' one versioned parameter map. RoBMA summaries, plots, density estimation, and
#' hypotheses use its semantic catalog view.
#' Canonical random-effect names have the form
#' `(formula) owner: quantity(arguments)`; the formula prefix can be omitted.
#' A bare formula or unnamed one-entry list also suppresses its redundant owner,
#' while an explicitly named one-entry list or a multi-component model retains
#' owners for disambiguation. Parentheses contain coefficient or parameter
#' names and square brackets contain factor or index levels. Examples include
#' `sd(intercept)`, `cor(group[sensitivity],group[specificity])`,
#' `study: sd(intercept)`, `sd_total`, and `var_prop(study)`.
#' A custom `random_variance_allocation()` always requires a stable internal
#' `name`; its `display_name` and `component_names` separately control these
#' public owner and component labels.
#'
#' Formula-random correlations are always public `cor` quantities. Compact
#' scalar `rho` coordinates and LKJ primitives are internal backend details,
#' not aliases. This convention is distinct from the released legacy
#' `cluster` interface, where `rho` remains the public variance-allocation
#' parameter described in \code{\link{prior_specification}}.
#'
#' ## Ways to customize random-effect priors
#'
#' There are three increasingly explicit customization levels:
#' \enumerate{
#'   \item Supply an ordinary positive prior distribution to
#'   `prior_heterogeneity`. This changes the base SD prior but preserves all
#'   automatic allocations and default correlation priors.
#'   \item Supply a partial `prior_random()` object containing only correlation,
#'   contrast, monitoring, new-level prediction, or parameterization settings.
#'   RoBMA supplies its base SD prior and every required variance allocation.
#'   \item Supply a `prior_random()` object containing scale architecture. The
#'   object is complete and RoBMA does not add or merge scale defaults.
#' }
#'
#' A `prior_random()` object is partial only when all scale-defining fields are
#' absent. The following fields make it complete: top-level or block-specific
#' `sd`; `random_covariance(sd = ...)`; block `sd_source`; block `terms`
#' containing term-specific SD overrides; or `allocation`. Once any of these is
#' present, the user owns the entire scale architecture and must provide every
#' required SD or allocation source.
#'
#' For example, the following partial specification changes only the factor
#' basis of block `study`:
#' \preformatted{
#' random_prior <- prior_random(
#'   study = random_block(contrasts = c(group = "independent"))
#' )
#' }
#' Combined with `random = ~ us(0 + group | study)`, this creates one
#' correlated random coefficient per `group` level. RoBMA supplies the
#' UISD-scaled `sd_common` and mean-variance Dirichlet allocation; BayesTools adds
#' the default `LKJ(1)` correlation prior after resolving the coefficient count.
#' No empty `random_covariance()` or manual allocation is needed.
#'
#' A correlation-only customization is partial as well. Use
#' `random_covariance(cor = prior_lkj(...))` for US/UN and
#' `random_covariance(cor = ..., cor_scale = ...)` for CS/HCS, AR1/HAR, or CAR.
#' An explicitly supplied scalar `cor` prior defaults to the Fisher-z scale;
#' `cor_scale = "cor"` instead places it directly on the raw correlation, and
#' `cor_scale = "logit"` uses a logit transform of the structure's admissible
#' interval. The `structure` in the formula remains authoritative, and a
#' conflicting covariance override is rejected.
#'
#' Complete specifications can define global SD defaults, block-specific
#' \code{\link{random_block}} settings, term-specific SDs, external SD sources,
#' and allocations built with \code{\link{random_variance_allocation}}. Use a
#' complete specification when the desired scale hierarchy differs from
#' RoBMA's automatic total/mean-variance allocation; do not add an SD merely to
#' customize a correlation or contrast.
#'
#' The `parameterization` argument and block-specific overrides select
#' `"noncentered"`, `"centered"`, or `"auto"` backend representations without
#' changing the specified standard-deviation, allocation, or correlation prior
#' distributions. Explicit centering requires strictly positive backend-owned
#' standard-deviation coordinates; BayesTools rejects unsupported centered
#' combinations and lets `"auto"` fall back to the noncentered representation.
#' Monitoring and new-level prediction settings similarly affect saved
#' quantities and prediction policy, not the statistical prior. Keep latent
#' effects monitored when bridge sampling or conditional prediction is needed.
#' See \code{\link{prior_random}} for every constructor argument and complete
#' examples.
#'
#' ## Location-scale random-effect models
#'
#' When `scale` is used together with `random`, the scale model targets the
#' row-wise total random-effect standard deviation. With one top-level
#' random-effect component, a single `scale` formula is unambiguous. With
#' multiple top-level random-effect components, `scale` must be a named list
#' whose names match the random-effect components. The random-effect allocation
#' rules then split the row-wise standard-deviation source rather than a scalar
#' total SD parameter. The `prior_scale` argument controls the scale-regression
#' coefficients. A custom `prior_random()` containing explicit SD architecture
#' must agree with these row-shaped scale sources; incompatible sources are
#' rejected rather than silently replaced.
#'
#' ## Inspecting the resolved priors
#'
#' The completed prior depends on the resolved formula columns and index levels.
#' Use `only_priors = TRUE` to perform data, formula, and prior resolution
#' without MCMC, then inspect the result with
#' `print_prior(object, parameter = "random")`. The printed object shows the
#' resolved block names, SD sources, allocation hierarchy, correlation priors,
#' contrasts, and policies that will be used for fitting. Its mathematical
#' `sigma_total` and `sigma_common` labels describe the prior architecture;
#' fitted public selectors use `sd_total` and `sd_common`, respectively.
#'
#' @examples
#' \dontrun{
#' random_prior <- prior_random(
#'   study = random_block(contrasts = c(group = "independent"))
#' )
#'
#' priors <- brma.mv(
#'   yi = yi, V = V, ni = ni,
#'   mods = ~ 0 + group,
#'   random = ~ us(0 + group | study),
#'   prior_mods = group_prior,
#'   prior_heterogeneity = random_prior,
#'   data = dat, measure = "GEN",
#'   only_priors = TRUE
#' )
#' print_prior(priors, parameter = "random")
#' }
#'
#' @seealso
#' \code{\link{prior_specification}},
#' \code{\link{random_effect_formula_tags}},
#' \code{\link{brma.mv}},
#' \code{\link{prior_random}},
#' \code{\link{random_block}},
#' \code{\link{random_covariance}},
#' \code{\link{random_variance_allocation}}
#'
#' @aliases random_effect_prior_specification random_prior_specification
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
#' term-specific components. This argument can be specified only when `mods` is
#' specified.
#' @param prior_scale prior distribution(s) for alternative scale-regression
#' components. A single prior applies to all terms; a named list can specify
#' term-specific components. This argument can be specified only when `scale`
#' is specified.
#' @param prior_heterogeneity_allocation prior distribution(s) for the
#' alternative cluster-level heterogeneity allocation component(s). This
#' argument can be specified only when `cluster` is specified.
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
#' components. This argument can be specified only when `mods` is specified.
#' @param prior_scale_null prior distribution(s) for null scale-regression
#' components. A single prior applies to all terms; a named list can specify
#' term-specific components. This argument can be specified only when `scale`
#' is specified.
#' @param prior_heterogeneity_allocation_null prior distribution(s) for the
#' null cluster-level heterogeneity allocation component(s). This argument can
#' be specified only when `cluster` is specified.
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
# based on Eq. 6 of RÃ¶ver, C., Bender, R., Dias, S., Schmid, C. H., Schmidli, H., Sturtz, S., ... & Friede, T. (2021). On weakly informative prior distributions for the heterogeneity parameter in Bayesian randomâ€effects metaâ€analysis. Research Synthesis Methods, 12(4), 448-474. 10.1002/jrsm.1475
.get_unit_information_sd.norm <- function(sei, ni) {

  UISD <- sqrt(sum(ni) / sum(sei^(-2)))

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
