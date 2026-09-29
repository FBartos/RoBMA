## version 4.1.0 (IN PROGRESS)
### Breaking changes
These changes affect code and saved objects written for RoBMA 4.0.0.
- installation and saved objects:
  - requires R >= 4.3.0 (was 4.0.0), BayesTools >= 0.3.1 (was 0.3.0), and
    loo >= 2.10.0, and imports Matrix (>= 1.6-2), ps, and statmod.
  - requires JAGS 4.x (>= 4.3.1, < 5.0.0). The RoBMA JAGS module implements
    the JAGS 4 module interface and does not compile against JAGS 5.
    Installation stops with a message naming the reported version when
    `configure` finds another major version through pkg-config or the
    `JAGS_MAJOR` declared by the selected headers. On Windows the build
    selects the newest installed `JAGS-4.*` when `JAGS_ROOT` is unset and
    stops when `JAGS_ROOT`, `JAGS_VERSION`, `JAGS_MAJOR_VERSION`, or the
    headers name another major version.
  - models fitted with RoBMA 4.0.0 must be refitted: their fits lack the
    parameter map and the fitted-object contract of BayesTools 0.3.1.
    `summary()`, `print()`, `plot()`, `as_draws()`, `hypothesis()`,
    `update()`, and the other functions that read the fitted parameters stop
    on them with an error of class `BayesTools_refit_required`. RoBMA's own
    refit requests have the class `RoBMA_refit_required` with the parent
    class `BayesTools_refit_required`, so callers match
    `BayesTools_refit_required` for every stale fit (objects without
    posterior samples, such as `only_priors = TRUE` objects, are asked to be
    fitted instead).
- removed and reordered arguments and methods:
  - removes the `normal_approximation` argument of `marginal_means()`.
    Marginal-means Bayes factors are Savage-Dickey density ratios with a
    kernel density estimate of the posterior ordinate (see the qCMDE/IWMDE
    features for the other estimators). `n_samples` moves from the fourth to
    the third position; pass it and the later arguments by name.
  - removes the `logLik()` method of `brma` objects; `AIC()` and `BIC()`
    stop with an explanation. Use `log_lik()` for pointwise posterior
    log-likelihood draws.
  - `update()` stops when `seed` is supplied: extended chains continue their
    stored JAGS random-number state and cannot be reseeded (BayesTools'
    `JAGS_extend()` no longer takes `seed`). Omit `seed`.
  - `RoBMA()` and `bselmodel()` take the new `selection` and
    `selection_control` arguments before `sample`; pass `sample` and the
    later arguments by name.
  - the generic `as_draws()` of `brma_samples` objects returns a
    `draws_array`, so chain dimensions are explicit; `as_draws_matrix()`
    remains available and keeps its chain-count metadata.
  - cluster-unit LOO and WAIC of binomial and Poisson GLMMs
    (`unit = "cluster"`) stop until certified nested adaptive quadrature is
    available; use `unit = "estimate"`.
- stricter inputs:
  - `BMA.norm()`, `RoBMA()`, `BMA.glmm()`, `BMA.mv()`, and `RoBMA.mv()` with
    their default mean-difference factor contrasts, and any fit with
    explicitly set mean-difference or orthonormal contrasts (e.g.,
    `brma(..., set_contrast_factor_predictors = "meandif")`), stop when a
    moderator or scale formula interacts a factor with a predictor whose own
    term is missing, such as `g:x` in `mods = ~ g + g:x` or `mods = ~ g / x`.
    Such a term codes `g` by level indicators and has one coefficient per
    level; with BayesTools 0.3.0 its mean-difference prior was fitted as
    independent priors on these level coefficients and summarized as
    differences from the mean. Specify an independent prior for that term
    (`prior_mods = list("g:x" = prior_factor(..., contrast = "independent"))`),
    include the missing term (`mods = ~ g * x`), or set
    `set_contrast_factor_predictors = "treatment"`.
  - explicitly supplied `prior_mods`, `prior_scale`, and
    `prior_heterogeneity_allocation` stop when the model has no moderators,
    no scale formula, or no multilevel structure; 4.0.0 discarded them
    silently.
  - binomial GLMM base-rate priors (`prior_baserate`) must be point or beta
    priors and Poisson GLMM log-rate priors (`prior_lograte`) point or
    normal priors, including their truncations, so every accepted nuisance
    prior has a certified post-fit likelihood route.
  - nonfinite `weights` are rejected when the data are validated; positive
    fractional weights and the existing missing-value handling are kept.
- defaults:
  - the default Poisson GLMM `prior_lograte` is an independent
    `Normal(log(pooled crude rate), 1)` prior instead of a prior whose
    standard deviation depended on the data and the exposure-time unit. The
    new `RoBMA.options(default_lograte.sd = )` sets its standard deviation.
    Default Poisson GLMM fits change, and the prior's informativeness no
    longer depends on the exposure-time unit.
  - `BMA.glmm()` uses mean-difference factor contrasts by default
    (`set_contrast_factor_predictors = "meandif"`, was `"treatment"`), as
    the other model-averaging constructors do; moderator coefficients of
    calls without an explicit contrast change their interpretation.
  - `RoBMA.options(silent = )` defaults to `TRUE` (was `FALSE`) and is the
    default of `silent` in every constructor: `brma()`, `bselmodel()`,
    `bPET()`, `bPEESE()`, and `brma.glmm()` no longer print the JAGS output
    by default, and `BMA()`, `RoBMA()`, and `BMA.glmm()` follow the option
    (their `silent = TRUE` default ignored it).
  - `plot()` of `zplot()` results defaults to `plot_extrapolation = FALSE`:
    the plot shows the observed z-statistics and the fitted selected density
    only, and draws the extrapolated curve, which rescales the display by an
    inferred selection convention, only with `plot_extrapolation = TRUE`.
    `lines()` keeps `extrapolate = FALSE`.
  - convergence checks and autofit (`set_convergence_checks()`,
    `set_autofit_control()`, `update()`) check the product-space model
    indicators by default (new `check_indicators = TRUE`) through
    BayesTools' label-invariant state-occupancy diagnostics; 4.0.0 left them
    out. `monitor` restricts the checked parameters without dropping
    eligible indicators, and `allow_not_assessable = FALSE` fails a sampled
    parameter whose draws never change. Declared constants (point priors,
    reference and fixed weight-function bins) are reported as structural
    constants.
- results:
  - fits with a given `seed` give different draws than 4.0.0: BayesTools
    derives each chain's JAGS seed from `seed` through R's random-number
    generator instead of `seed + chain`, and samples two-bin cumulative
    weight-function priors through their exact Beta marginal. Seeded fitting
    no longer resets R's global random-number stream to a state determined
    by `seed` (BayesTools restores the caller's `.Random.seed` and
    `RNGkind()`); code that relied on that reset must set its own seed.
  - `predict()` has one two-axis contract for ordinary, multilevel,
    multivariate, and GLMM models: `type` selects fixed terms, latent
    effects, or observed responses, and the new `conditioning_depth` selects
    marginal (`"marginal"`, the default), fitted-cluster (`"cluster"`), or
    fitted-estimate (`"estimate"`) prediction. Marginal prediction is the
    default for implicit fitted designs and for equivalent explicit
    `newdata`; estimate depth includes the posterior uncertainty of the
    fitted latent effects, while `blup()` and `fitted()` keep conditional
    means. This changes the released same-data `type = "estimate"` and the
    clustered and GLMM `type = "response"` predictions. Marginal and cluster
    GLMM responses draw new nuisance base rates or rates from their priors;
    estimate depth keeps the fitted posterior nuisance rates.
  - `zplot()` targets the marginal, empirical-design posterior predictive
    distribution projected into z-space by default (new
    `conditioning_depth = "marginal"`); `conditioning_depth = "cluster"`
    gives the target that 4.0.0 used for multilevel models, and
    `"estimate"` adds same-effect posterior prediction that integrates the
    posterior uncertainty of each fitted latent effect. Fitted selection and
    inverse-probability extrapolation apply at the selected depth.
  - `funnel()` draws plug-in contours: location and bias parameters at their
    within-model posterior means and heterogeneity at its root mean square
    SD, `sqrt(E(tau^2))`, with the contours of the complete joint models of
    model-averaged fits mixed by their posterior model probabilities before
    inversion. The posterior-averaged contours that 4.0.0 drew for
    model-averaged publication-bias fits are drawn by the new `bfunnel()`.
    Outcome mode needs an intercept-only model with one common marginal
    heterogeneity distribution; other models use LOO-PIT residual mode.
    Moderator-dependent PET/PEESE curves and outcome-mode funnels of models
    with location or scale predictors stop, as a single unconditional curve
    is not a fitted estimand there. GLMM outcome-mode contours are marked as
    descriptive normal effect-size approximations.
  - `cooks.distance()` reports the unscaled squared posterior Mahalanobis
    distance, following metafor's chi-square-based convention; values were
    previously divided by the fixed-effect model rank and increase by that
    rank when it exceeds one.
  - `vif()` averages the posterior coefficient covariance over the draws
    instead of plugging the posterior mean heterogeneity into one covariance
    matrix.
  - estimate-unit LOO and WAIC of Gaussian models are deletion-conditioned
    scores, p(y_i | y_{-i}, theta), with every local Gaussian random effect
    integrated: an estimate of a multilevel model is scored given the other
    estimates of its cluster, where 4.0.0 conditioned on the sampled cluster
    effects. Sampled, marginalized, and mixed random-effect
    parameterizations share one estimand; LOO-PIT and studentized residuals
    use the matching conditional distribution.
  - marginal-means Bayes factors follow the Savage-Dickey rules of
    BayesTools 0.3.1 (see also Fixes): a null hypothesis outside the
    posterior draws gives a finite Bayes factor instead of `Inf`; and a
    marginal mean structurally fixed at the null hypothesis, or whose prior
    density at the null is a point mass, zero (a null outside the prior
    support), infinite, undefined, or without an exact value, has an `NA`
    Bayes factor with the reason among the table's warnings, while the other
    marginal means are computed (4.0.0 reported the density ratio, with a
    warning that it is likely invalid for a mean fixed at the null). A point
    hypothesis on such a mean with `hypothesis()` stops with the BayesTools
    class (e.g., `BayesTools_zero_ordinate` or
    `BayesTools_inexact_ordinate`).
  - `ranef()` of multilevel models returns the cluster-level effects with one
    column per cluster by default (`u_cluster[<cluster>]`); 4.0.0 returned
    one column per estimate (`u_cluster[<cluster>|<estimate>]`), which
    `expand = TRUE` keeps.
- labels, tables, and posterior draws:
  - `as_draws()`, `as_draws_array()`, `as_draws_df()`, `as_draws_list()`,
    `as_draws_matrix()`, and `as_draws_rvars()` of fitted models return the
    model-parameter schema by default: backend-private latent effects,
    random-effect simulation and Cholesky factors, and prior-parameterization
    variables are omitted (4.0.0 returned every monitored column), and
    parameters with point priors are exact constant columns. The new
    `include_auxiliary = TRUE` adds the backend variables.
  - labels are rendered by BayesTools from the parameter catalog: square
    brackets after a factor term hold a level label (`mu_g[5]`, `g[5]`, and
    `(mu) g[5]` mean the level labelled "5"), and contrast coefficients that
    are not a level (mean-difference, orthonormal, and ordered coefficients)
    are labelled `g{j}` in untransformed output.
  - summary tables of `BMA()`/`BMA.norm()` and `RoBMA()` fits label the
    intercept of a scale formula `exp(intercept)` (was `intercept`), as those
    of the other constructors already did, which makes its SD scale
    explicit; the table note describes exponentiated scale slopes as
    multipliers of the SD.
  - a location intercept fixed at zero is omitted from meta-regression
    summaries, stored coefficients, marginal means, plots, density
    estimates, and hypotheses whenever moderators are present;
    intercept-only models keep it.

### Features
- multivariate and multilevel models with known sampling covariance:
  - adds `brma.mv()`, normal-likelihood meta-analysis with a known sampling
    covariance matrix `V` (latent, whitened, and block-MVN backends through
    `known_v_parameterization`) and BayesTools formula random effects
    through `random` (e.g., `random = ~ 1 | study / esid`, `us()`, `cs()`,
    `hcs()`, AR1, HAR, and CAR structures, and plain random slopes such as
    `(1 + x | study)` with `us()` semantics), with metafor-style `R` and
    `Rscale` for a known random-effect group covariance of random-intercept
    blocks.
  - adds `BMA.mv()` for product-space model averaging with the complete
    `brma.mv()` workflow. Independent random-component gates multiply their
    allocated slab variances without changing the fitted component scales;
    `tau_total` and `tau2_total` report the realized gated aggregate,
    including the all-off zero branch, and `tau2_prop(...)` the active
    component shares conditional on positive total heterogeneity.
  - adds `RoBMA.mv()`, which combines the `BMA.mv()` product space with
    unadjusted, selection-model, PET, and PEESE branches, and `bPET.mv()`,
    `bPEESE.mv()`, and `bselmodel.mv()`. PET uses `sqrt(diag(V))` and PEESE
    `diag(V)` as row-level bias predictors while keeping the full known
    sampling covariance in the likelihood.
  - all `.mv` constructors share the summaries, predictions (with the new
    `V_new` argument of `predict()`, required for explicit known-`V`
    response predictions), LOO/WAIC, diagnostics, plots, hypotheses,
    posterior conversion, and `update()` of the other RoBMA models; bridge
    sampling marginalizes their Gaussian location random effects exactly. An
    omitted or `NULL` `random` is a fixed-effect model, and random terms
    must be declared before heterogeneity priors or scale formulas are
    supplied.
  - `V` can be a dense, diagonal, or block-list matrix, or a
    `known_v_factor()` declaration `diag(d) + UU'`. Ordinary covariance
    matrices, such as `metafor::vcalc()` output, are accepted directly:
    block-constant and other exact `diag(d) + UU'` structures are
    recovered to working precision for the exact selection routes, and
    other blocks keep the general route. Numerically symmetric matrices are
    made exactly symmetric once at input; exactly singular blocks that
    reach working precision (e.g., `vcalc(rho = 1)`) are accepted, and
    positive-semidefinite `V` is supported when priors or marginalized
    random effects regularize every null direction. Negative eigenvalues,
    correlations above one, and singular designs without covering model
    structure fail input validation.
  - named `scale` lists target concrete random-effect blocks, including an
    unambiguous terminal nested grouping name such as `esid` in
    `random = ~ 1 | study / esid`; untargeted blocks keep their ordinary SD
    priors.
- random effects:
  - re-exports BayesTools' `prior_random()` interface with `random_block()`,
    `random_covariance()`, `random_monitor()`, `random_new_levels()`,
    `random_sd_source()`, `random_variance_allocation()`,
    `allocation_ref()`, and `prior_lkj()`, documented under
    `?prior_random`. Partial specifications are completed with RoBMA's
    UISD-scaled SD and variance-allocation defaults; correlation defaults
    are BayesTools' (`LKJ(1)` for US/UN and uniform priors over the complete
    admissible range for CS/HCS, AR1/HAR, and CAR).
    `random_block(parameterization = "mean_centered")` samples study means
    of eligible retained scalar random intercepts.
  - random-effect quantities have semantic names in summaries, plots,
    density estimates, hypotheses, and `as_draws()`: `tau`, `tau2`, `rho`,
    `tau_total`, `tau2_total`, `tau_common`, `tau2_common`, `tau2_prop`,
    `tau_mult`, and `tau2_mult`. A sole random intercept prints as `tau` or
    `study: tau`, non-intercept coefficients as `study: tau(x)`; a bare
    formula or unnamed one-entry list omits the owner prefix (`rho(...)`,
    `tau_common`), lists with several unnamed components are
    `component 1`, `component 2`, ..., and the shared correlation of `cs()`
    and `hcs()` blocks is also selected by the pairwise aliases of its
    levels (`rho(outcome[a],outcome[b])`). Owner-free shorthand is accepted
    when it resolves uniquely.
  - `print_prior()` prints the structured random-effect prior once,
    `summary()` reports the quantities aligned with the prior
    specification, and `summary_heterogeneity()` (new `component` argument)
    the aggregate variances, allocation-derived component SDs and
    variances, multipliers, and correlations in their owning component.
    Original-scale correlations `rho(...)` are summarized, tested, and
    plotted over their defined draws (a correlation is undefined where one
    of its SDs is zero), with a table note on the number of defined draws;
    the exported `as_draws()` column has `NA` in those draws.
  - random-effect inclusion is reported in the Component Inclusion table
    (rows `Heterogeneity: ...`, before Publication Bias) and by
    `interpret()`, and the inclusion gates are exported by `as_draws()` as
    columns named after their SD with an `_indicator` suffix (e.g.,
    `study: tau_indicator`, and `tau_total_indicator` for the component of an
    allocation-SD prior).
  - `ranef()` gains `component`, `simplify`, and metafor-compatible `expand`
    (one column per unique grouping-level contribution by default,
    `expand = TRUE` for observation-aligned output); for `brma.mv()`
    random-formula models it returns a flat list keyed by canonical block
    names, with nested blocks in metafor's outer-to-inner order.
    `pooled_heterogeneity()` gains `component`.
  - `plot()` of random-effect quantities (with `component = "random"`)
    draws exact BayesTools prior densities: the LKJ marginal of fitted-scale
    correlations, Beta marginals of Dirichlet allocation fractions, and the
    allocated SD, SD-multiplier, and variance-multiplier priors. Quantities
    without an exact prior density (original-scale correlations that
    combine the correlation with the SDs of a scaled block) are drawn
    without the prior curve and with a warning of class
    `BayesTools_prior_curve_unavailable`. SD and variance multipliers accept
    `transform = "LOG"`.
- selection models:
  - `bselmodel()`, `RoBMA()`, `bselmodel.mv()`, and `RoBMA.mv()` fit
    conditional Gaussian vector selection models specified with
    `selection = selection_model(...)` (re-exported from BayesTools; also
    `prior_weightfunction(model = )`): `estimate_random_effects`,
    `other_random_effects`, and `known_sampling_variance` are each
    conditioned on or integrated (`"condition"`/`"integrate"`, by default
    integrate, condition, integrate), weights follow `weight_rule`
    (`"product"`, the default, or `"best"`), and `group` binds publication
    groups. Explicit weight-function priors keep their own specification.
  - best-p-value weights use the smallest actual p-value within each
    publication event (unequal standard errors, one- and two-sided
    geometry, nonmonotone weights, and zero-weight bins). Publication groups
    come from an explicit prior column or a supported constructor cluster,
    never from covariance blocks or random-effect groups; product models
    ignore publication grouping. Selection thresholds always use the
    original sampling standard errors, and the sampling covariance defines
    the model whatever its representation.
  - one execution plan serves fitting, bridge sampling, posterior
    likelihoods, densities, predictions, LOO, and z-plots: independent
    factors are analytic, supported low-rank factors use deterministic
    quadrature (nested, tensor, or Smolyak rules) with a controlled QMC
    fallback, and every integral reports its error diagnostics.
    `set_selection_likelihood_control()` (arguments `selection_control` of
    the constructors, and `density_control$integration_control` of
    densities and hypotheses) sets the QMC controls. Models that condition
    on every source with positive weights use their exact Gaussian
    likelihood, with a warning when fitting starts. Non-unit observation
    `weights` require integrated sampling variation and independent product
    factors.
  - selection normalizers are cached under one budget
    (`RoBMA.options(selection.cache_max_bytes = )`, by default the smaller
    of 4 GiB and a quarter of the available RAM, shared by all chains of a
    fit); `selection.cache_retain = TRUE` keeps the caches with the fitted
    object (released by `remove_selection_cache()`), and
    `selection.sampler = "coarse_corrected"` proposes slice moves from
    coarse normalizers (`selection.coarse_grid`,
    `selection.coarse_max_rules`) and corrects each proposal against the
    full target. `selection_cache_info()` and `selection_sampler_info()`
    report the counters.
  - adds `selection_sensitivity_diagnostics()` and
    `add_selection_sensitivity_diagnostics()`, which compare the explicit
    conditioning models at the same weighting rule and publication
    partition on request; stored comparisons are refreshed when chains are
    extended.
  - `zplot()` of multivariate models, estimate-conditioned models, and
    best-rule models uses the normalized pre-selection Gaussian marginal as
    the bias-adjusted reference and reports its significance probability as
    the EDR; missing-study counts are unavailable, because relative weights
    do not identify absolute publication probabilities. Independent and
    conditional-cluster univariate product models keep the conventional
    inverse-weight diagnostic when estimate effects are integrated. Its
    numerical
    integration takes `integration_control`, and `plot()` and `lines()` of
    `zplot()` results take `integration_control`, `parallel`, and `cores`
    (1,000 posterior draws by default for `brma.mv()` models).
  - `funnel()` and `bfunnel()` draw outcome-mode contours of correlated
    known-`V` selection models under `weight_rule = "product"`. Selected
    funnel contours and regression-plot sampling intervals stop for
    selection configurations that the scalar calculation does not support;
    bias-adjusted contours (`sampling_bias = FALSE`) remain available.
- hypothesis tests:
  - adds `hypothesis()` for fitted models and for `marginal_means()`
    results, with the aliases `bf_hypothesis()` and `BF_hypothesis()`, and
    `hypothesis_quantities()`, which lists the names `hypothesis()` accepts
    and the tests it runs on each quantity (columns `alias`, `parameter`,
    `component`, `term`, `bracket`, `point_test`, `direction_test`,
    `contrast_test`, `point_test_methods`, `contrast_test_methods`, and
    `reason`). The `alias` column lists only names that `hypothesis()`
    resolves to the row's quantity with the row's `component` (a name shared
    by several quantities of a component, such as the `intercept` of several
    scale formulas, is omitted); the alias of a factor term also names its
    levels (`<alias>[<level>]`) with that `component`.
  - statements are point, interval, and directional hypotheses, with
    explicit comparisons (`"mu = 0 vs mu > 0"`) and several statements per
    call that may target different parameters. They reference parameters,
    factor levels by their labels (`g[a]`), linear combinations of the
    levels of one factor term (`g[a] = g[b]`, `2 * g[a] = 0.1`), unquoted
    formula interactions (`alloc:ablat[random]`), and the display labels of
    the summary tables (`(mu) intercept`, `(log_tau) intercept`,
    `exp(intercept)`, `(mu) g[a]`, `(mu) study: sd(intercept)`); constants
    may stand on the left (`0 > mu`). `component` (`"mods"`/`"location"`,
    `"scale"`, or `"random"`) selects between parameters of several model
    components; with it, a level of a factor alias shared by the location and
    scale formulas (e.g., `Preregistered[Pre-Registered]` for
    `mods = ~ Preregistered, scale = ~ Preregistered`) resolves within that
    component. Level-alias labels of mean-difference factors that
    `summary()` prints (e.g., `g1[dif: 5]`) are not accepted; use the level
    name (`g1[5]`).
  - every statement is planned once: the plan records the tested target
    (its weights on the fitted coefficients, its BayesTools prior density,
    and its declared atoms), classifies each point value with
    `BayesTools::prior_ordinate_status()`, and states which of KDE, qCMDE,
    IWMDE, and the normal approximation can evaluate it. `hypothesis()`, its
    marginal-means method, and `hypothesis_quantities()` execute or render
    these plans only, so a listed test is run and a refused statement stops
    with the same reason on every route.
  - point hypotheses need an exact, finite, positive prior ordinate at the
    null. They are available for the levels of factors of every contrast
    (treatment, independent, mean-difference, orthonormal, and ordered),
    including levels and level contrasts of multivariate t factor priors
    (`prior_factor("mt", ...)`, an exact univariate t) and, in
    model-averaged fits, levels off their null atom (the exact ordinate of
    the mixture of the factor priors; linear combinations of levels need
    `conditional = TRUE`); for the original-scale intercept of a
    log-intercept scale regression with normal slope priors (other slope
    priors are refused with `BayesTools_inexact_ordinate`); and for gated
    random-effect SDs, variances, totals, and variance proportions away from
    their point masses (the exact continuous prior density). A variance
    point hypothesis is evaluated through its SD, so `tau2 = v^2` and
    `tau = v` give the same Bayes factor. Point nulls on an exact support
    boundary, such as `tau2_prop(study) = 0`, use the one-sided prior
    ordinate when BayesTools classifies it as exact, finite, and positive,
    for KDE (with its boundary-reflected kernel), qCMDE, and IWMDE.
  - statements on a factor level fixed by the contrast (the treatment
    reference level) are refused with one reason whatever the density
    method; contrasts with the reference level remain available. Ordered
    levels whose Dirichlet share gives an infinite prior ordinate at the null
    (e.g., `g[b] = g[a]` for the first of several increments) stop with
    `BayesTools_infinite_ordinate`. Allocation-derived component SDs at zero
    are refused as nonregular product boundaries with a pointer to
    `tau2_prop(...) = 0` or `tau2_mult(...) = 0`, and point tests whose
    transformation has a singular Jacobian at the requested value fail with
    a structural explanation. SD components of nested allocations, which
    have no exact prior density, take region prior probabilities from prior
    draws and refuse point tests; original-scale correlations of scaled
    `us()` blocks support region hypotheses (e.g.,
    `"rho(intercept,x) > 0"`) over their defined draws.
  - refusals are classed conditions. Tests that cannot be computed have the
    class `RoBMA_hypothesis_unavailable` with `RoBMA_hypothesis_fixed` (a
    quantity fixed by the model, such as the treatment reference level),
    `RoBMA_hypothesis_target` (no supported target, e.g., publication-bias
    parameters, weight-function coordinates such as `omega[0,0.025]`,
    inclusion indicators, `inclusion(study)`, or latent cluster effects), or
    `RoBMA_hypothesis_method` (the density method is unavailable for the
    target), and BayesTools refusals keep their classes (e.g.,
    `BayesTools_point_mass_at_null`, `BayesTools_infinite_ordinate`,
    `BayesTools_inexact_ordinate`, `BayesTools_linear_target_unavailable`).
    qCMDE/IWMDE refusals whose cause `plot()` refuses too also have the
    classes of that refusal (see qCMDE/IWMDE below). The fitted model's
    qCMDE/IWMDE capability is checked before the quantity's own method
    refusals, such as "qCMDE/IWMDE point hypotheses for factor parameters
    must specify a level" for a point statement on a whole factor term.
  - statements that have to be restated, or an argument that has to
    change, stop with `RoBMA_hypothesis_statement` without
    `RoBMA_hypothesis_unavailable`: unsupported point expressions (e.g.,
    `2 * mu = 0`), nonlinear expressions of factor levels, one quantity
    named by several names in one statement (e.g., `mu > 0 & intercept < 1`
    or `g[a] > mu_g[a]`; use one name), factor contrast coefficients such
    as `g{1}` (a selector of a coefficient that is a level, such as `g{1}`
    of a treatment factor, also has the class
    `BayesTools_selector_unavailable`, which names the level form), unknown
    names and levels (also `BayesTools_parameter_not_found` and
    `BayesTools_parameter_resolution_error`, with the fields `alias` and
    `available`), statements without parameter references (also
    `BayesTools_hypothesis_no_parameters` and
    `BayesTools_parameter_resolution_error`, on fitted models and marginal
    means alike), ambiguous references
    (`RoBMA_hypothesis_ambiguous`; set `component`, or `parameter` for
    marginal means), and references to parameters of another component
    (`RoBMA_component_mismatch`). Every other refusal of a statement's
    references by BayesTools (`BayesTools_parameter_resolution_error`) also
    stops with `RoBMA_hypothesis_statement` before its classes and with its
    fields, on fitted models and marginal means, also when BayesTools raises
    it only when it evaluates the statement (e.g., the whole factor term
    next to one of its levels, `"g[a] > mu_g"`). With an explicit
    `component`, a name known only outside it, or a display label next to
    another component's name, is a component mismatch, also in statements
    with factor levels (e.g., `"log_tau_g[a] > mu_g[a]"` with
    `component = "mods"`). The exception is the alias of a factor term
    of both the location and the scale formula next to a name unknown
    within `component`, such as another component's name or a display
    label (e.g., `"g > log_tau_intercept"` with `component = "mods"`): the
    statement stops with `RoBMA_hypothesis_statement` and the BayesTools
    refusal of its first unresolved reference, the ambiguity of the alias
    (`BayesTools_parameter_ambiguous`) or the unknown name
    (`BayesTools_parameter_not_found`, e.g., for
    `"g[a] > log_tau_intercept"`).
  - `hypothesis()` warns when a parameter with null and alternative
    components is tested on the full product-space ensemble because
    `conditional` was omitted; set `conditional = FALSE` to keep that test
    without the warning, or `conditional = TRUE` to test only the models in
    which the parameter is active.
  - hypothesis tables name repeated rows by statement (`g1 (1)`,
    `g1 (3)`), keep the plain label for single statements, and print the
    public labels of the tested quantities.
- posterior densities and point ordinates (qCMDE and IWMDE):
  - adds likelihood-aware qCMDE and IWMDE estimates of posterior densities
    and point ordinates through `density_method` (`"KDE"`, `"qCMDE"`, or
    `"IWMDE"`) and `density_control` in `plot()`, `lines()`,
    `hypothesis()`, `marginal_means()`, and marginal-means plots and
    hypotheses. qCMDE is the default of `hypothesis()` on fitted models;
    plots and marginal means use KDE unless requested. `marginal_means()`
    gains `density_method`, `density_control`, and the selectors
    `parameter`, `type`, and `levels`, which restrict the precomputed
    densities (and, for `parameter` and `levels`, the conditional ordinates
    of the inclusion Bayes factors); marginal-means plots and hypotheses use
    the stored method unless overridden, and `bf = FALSE` skips the
    ordinates.
  - `density_control` takes `samples` (a fixed simple-random sample of
    posterior rows, drawn without replacement from the eligible continuous
    rows), `n_points`, `normalization_points`, `normalization_prob`,
    `target_relative_mcse`, `display_grid`, and `integration_control`, and
    rejects other settings. qCMDE places each row's normalization range at
    that row's conditional quantiles at `normalization_prob` (which must be
    below 1), in closed form where a Gaussian likelihood kernel meets a
    normal prior, from a Gaussian tail bound with another proper prior, and
    otherwise by extending the range until the endpoint tail estimate meets
    the target, and normalizes each row on nested uniform grids.
  - adds `density_diagnostics()` with the classes `RoBMA_density_diagnostics`
    and `RoBMA_density_ordinate_error`: row-sampling uncertainty and
    selected-row MCMC error are reported separately, qCMDE reports its
    largest row truncation (exact, bound, or estimate) and the ordinate
    error bound t/(1 - t), density curves are graded over their 5%-95%
    bulk, and failed point computations are kept with their reason. The
    diagnostics never change the estimator's row sample. Estimation fails
    when a selected row lacks a valid normalizer, proposal density, or
    finite contribution. Results without attached diagnostics (KDE
    ordinates, tests without a point hypothesis) stop with the class
    `RoBMA_density_diagnostics_unavailable`.
  - supported targets are fitted coefficients and factor levels (on the
    displayed coefficient scale, through the exact fitted-to-original
    transform), marginal means, and the random-component quantities of
    `brma.mv()` family models; known-`V` formula random effects are
    marginalized, allocation endpoints use their exact Dirichlet densities,
    GLMM ordinates condition on the sampled estimate-level effects and
    nuisance rates, and each estimate carries RoBMA's schema and version
    provenance.
  - qCMDE/IWMDE requests that the fitted model does not support stop before
    any density is estimated with an error of class
    `RoBMA_density_method_unavailable` and a class naming the cause:
    `RoBMA_density_method_random_unknown_v` (`brma.mv()` random-formula
    models without known `V`), `RoBMA_density_method_scale_components`
    (component-specific scale formulas, a named `scale` list),
    `RoBMA_density_method_glmm` (IWMDE for binomial and Poisson GLMMs),
    `RoBMA_density_method_conditional_random` (conditional random-effect
    plots and statements), `RoBMA_density_method_random_target` (a
    random-effect quantity without a supported scalar coordinate), or
    `RoBMA_density_method_original_scale` (an original-scale coefficient
    whose fitted map is nonlinear, such as the exponentiated intercept of a
    scale formula with standardized predictors; use
    `standardized_coefficients = TRUE`). `hypothesis()` adds these classes
    after its own, and `hypothesis_quantities()` reports the same reasons.
- predictions, diagnostics, and model comparison:
  - adds `log_lik()` for pointwise posterior log-likelihood draws, as used
    by LOO and WAIC.
  - re-exports `loo::loo_model_weights()` and adds a `brma` method that
    computes stacking or pseudo-BMA weights from compatible stored LOO
    results, also for fits with different numbers of posterior draws.
    `loo_compare()` tables of `brma` fits have the print class
    `compare.loo.brma` and keep the numeric `compare.loo` matrix.
  - `add_marglik()` gains `parallel`, `cores` (defaulting to the fitted
    settings), `repetitions`, `method`, `maxiter`, and `silent`; repeated
    bridge estimates are kept in Bayes-factor and posterior-probability
    comparisons (Bayes factors from the median log marginal likelihood,
    posterior probabilities with one row per repetition).
  - fully fixed point-prior models are supported for `brma()`,
    publication-bias, GLMM, multilevel, and `brma.mv()` models: fixed
    parameters stay in the posterior summaries, marginal likelihoods of
    zero-dimensional models are exact (with no bridge object, which
    `bridge_sampler()` states), and constant log-likelihood columns get
    exact uniform LOO importance ratios with zero Pareto-k.
  - `pooled_effect()` adds prediction intervals (`PI` columns) for one new
    true effect at the average design, with the heterogeneity of
    `pooled_heterogeneity()`, drawn without advancing the caller's
    random-number stream.
  - `hatvalues()`, `residuals()`, `rstandard()`, `qqnorm()`, and `vif()`
    gain `max_samples`, and `dfbetas()` and `covratio()` gain `component`
    and `parameter`.
  - `update()` accepts named `NULL` values to clear nullable
    convergence/autofit thresholds and monitors; omitted controls keep their
    fitted values.
- plots:
  - adds `bfunnel()`, which averages the posterior-draw sampling CDFs
    before inversion, and a `metafor::forest()` method for `brma` objects
    with `as_metafor_forest()` for the forest-plot data (prediction
    intervals target the pooled average design, or one new true effect for
    explicit `newdata`; shade/dist styles use the posterior predictive
    draws).
  - adds `lines()` methods for fitted models and marginal means, which
    overlay posterior (and prior) densities on an existing plot.
  - `plot()`, `plot_prior()`, `print_prior()`, and `plot_diagnostic()` gain
    `component` (`"auto"`, `"mods"`/`"location"`, `"scale"`, or
    `"random"`), and posterior plots select factor levels by label
    (`g[level]`), including prior overlays and qCMDE densities. These
    functions refuse quantities of the fitted model that are not model
    parameters (inclusion indicators, `bias_indicator`, weight-function
    coordinates such as `omega[1]`, `inclusion(...)`, and latent cluster
    effects such as `gamma[1]`) with a message naming the quantity and the
    classes `RoBMA_hypothesis_target` and `RoBMA_hypothesis_unavailable`,
    which `hypothesis()` gives them too, and a `parameter_mods` or
    `parameter_scale` selection of another component with
    `RoBMA_component_mismatch`.
  - `plot_prior(standardized_coefficients = FALSE)` stops for the prior of a
    term with several fitted coordinates (a factor with three or more
    levels or with independent contrasts, and the interactions of such a
    factor) whose original-scale coordinates the standardization of
    continuous predictors changes: a factor or factor-by-factor interaction
    next to its interaction with a standardized predictor, and that
    interaction term. The error has the classes
    `RoBMA_density_method_original_scale` and
    `RoBMA_density_method_unavailable` and points to
    `standardized_coefficients = TRUE`; it replaces the standardized-scale
    prior that was plotted silently. Terms whose coordinates are unchanged
    (factors and factor-by-factor interactions without an interaction with
    a standardized predictor) plot as before, and single-coordinate terms
    (such as a two-level factor with treatment, mean-difference, or
    orthonormal contrasts) are transformed as before.
  - `transform = "EXP"` applies to individual log-scale ratio
    meta-regression coefficients and to scale-regression coefficients (as
    multiplicative changes in heterogeneity), and `transform = "LOG"` to
    positive heterogeneity intercepts; KDE, qCMDE, and IWMDE curves use the
    exact change-of-variable Jacobian.
  - `regplot()` accepts moderators that appear only in a scale formula
    (a constant location effect with heterogeneity-dependent prediction and
    sampling intervals).
  - `radial()` supports `brma.mv()` models (with row-marginal
    heterogeneity); its bands for non-Gaussian and publication-bias models
    are documented as descriptive reference bands.
- summaries and output:
  - `as.data.frame()` and `data.frame()` convert the printable results
    (model and heterogeneity summaries, posterior, prediction, and pooled
    results, marginal means, model weights, z-plots, VIFs, and influence
    diagnostics) to one long data frame with leading `component` and
    `parameter` columns and `CI_`/`PI_` interval columns.
  - multivariate summaries follow the univariate layout (one table without
    predictors; common, regression, and scale sections with them), `BMA.mv()`
    and `brma.mv()` share one prior-parameterization summary, and models
    with several scale formulas label their scale rows, inclusion rows, and
    `summary_models()` rows by the targeted SD (`(Study: tau) x`).
- other additions:
  - adds the `Hoogeveen2023` data set: 106 analyst-team effect-size
    estimates of the relation between religiosity and self-reported
    well-being.
  - `set_contrast_factor_predictors = "independent"` defines one fixed
    coefficient per factor level, and formula random effects inherit that
    basis unless `random_block(contrasts = )` overrides it; fixed and random
    formulas warn before fitting when automatically independent factor
    columns overspecify their design.
  - adds the options `RoBMA.options(native_threads = )`, the thread budget
    of post-fit native kernels (by default one thread for models fitted
    serially and the fit's parallel settings otherwise), and
    `RoBMA.options(jags.worker_output = )`, which captures the stdout and
    stderr of parallel JAGS workers during fitting and extension.

### Fixes
- marginal-means Bayes factors use the exact Gaussian kernel sum at the null
  hypothesis (bandwidth `bw.nrd0()` of the continuous draws) as the
  posterior ordinate, through BayesTools 0.3.1. The 512-point grid estimate
  of 4.0.0 overestimated the posterior density at the null of long-tailed
  posteriors up to several-fold and so understated their Bayes factors.
- `predict()` stops and names a variable of the moderator or scale formula
  (or of a random-effect formula) that `newdata` lacks. 4.0.0 filled missing
  outcome-named columns (`yi`, `sei`, `ai`, `ci`, `n1i`, `n2i`, ...) with
  zeros for its data parser, so a formula that uses such a column as a
  predictor (e.g., `mods = ~ sei`) was silently evaluated at zero.
- a response formula without an intercept (`yi = yi ~ 0` or `yi ~ -1`)
  fixes the location intercept at zero, as `mods = ~ 0` does; 4.0.0 dropped
  the `0` and estimated the intercept.
- prior-only objects (`only_priors = TRUE`) evaluate moderator and scale
  formulas with standardized continuous predictors, as fitting does. These
  evaluations used the raw predictor values, so scale formulas gave wrong
  heterogeneity (in a test model, row-specific `tau` was 0.58 to 2.07 times
  the correct value) and moderator formulas gave wrong predictions.
- `summary()` of an `only_priors = TRUE` object returns the resolved priors,
  as for `only_data = TRUE`, instead of failing.
- weighted residual diagnostics use the original outcome variances, keeping
  the likelihood weights for fitting and PSIS deletion; internally
  standardized residuals combine the fitted projection with the original
  covariance, as metafor does with custom weights.
- plug-in funnel contours, model-averaged funnel contours, and deleted-row
  location-scale influence heterogeneity use the within-draw (within-model)
  root mean square heterogeneity, averaging integrated and conditioned
  heterogeneity on the variance scale.
- selected-normal funnel and regression-plot CDF tails are normalized
  jointly, so valid low-standard-error contours no longer fail through
  separately rounded normalization masses.
- `regplot(mod = "vi")` and `regplot(mod = "sei")` keep the precision of the
  prediction grid, derive the matching standard error or variance, and use
  pointwise precision in the sampling intervals.
- `pooled_heterogeneity()` is evaluated at the average expanded scale and
  random design; `summary_heterogeneity()` keeps the within-draw RMS over
  the observed design.
- cluster-level selection-model likelihood integrals are certified by
  successive 7-, 15-, and 31-node quadrature rules and fail explicitly when
  the requested accuracy is not established, instead of using heuristic
  curvature fallbacks.
- estimate-unit GLMM log-likelihoods (LOO, WAIC, `log_lik()`) use native
  joint adaptive Gauss-Hermite quadrature with a convergence-checked
  prior-CDF fallback instead of a fixed grid, including sparse binomial and
  Poisson columns, and boundary-matched quadrature for boundary-concentrated
  and endpoint-truncated beta base-rate priors, so valid sparse GLMM fits no
  longer fail in `log_lik()`, LOO, and WAIC. Zero and infinite conditional
  Poisson rates no longer produce `NaN` likelihoods.
- GLMM predictions accept the four-cell binomial count representation used
  in fitting (inconsistent duplicate totals are rejected), return raw
  binomial and Poisson responses as count samples without effect-size
  transformation metadata, and accept scalar JAGS nuisance-parameter names
  of models with one estimate.
- cached LOO/WAIC results and marginal likelihoods are checked against the
  current data, cluster membership, and fitted target before reuse, with
  versioned outcome fingerprints independent of R's serialization format,
  and ranked comparison tables keep the model names.
- `plot_weightfunction()` marks the observed p-values in the convention of
  the weight function: two-sided weight functions show two-sided p-values,
  where 4.0.0 marked one-sided p-values on their axis. Diagnostic plots
  respect their limits and clipping, and radial plots draw high-confidence
  intervals.
- transformed posterior and marginal-means plots take `xlim` on the
  displayed scale, so `transform = "EXP"` prior curves extend toward zero
  instead of starting at an exponentiated fitted-scale limit; `ylim2` and
  `ylab2` control the secondary point-mass axis, and independent factor
  coefficients get distinct default colors.
- `zplot()` honors explicit axes and cutpoints and keeps the histogram
  method's delegation.
- `ranef()` of ordinary multilevel models derives the cluster and estimate
  components from one joint Gaussian BLUP solve.
- `loo_weights()` explains, when several models are passed, that it returns
  single-model PSIS weights and points to `loo_model_weights()`.
- an empty analysis data set names the step that emptied it (an exhaustive
  `subset` or dropped missing values), and an unknown `mods` variable is
  named as an unknown random-effect variable is.
- fitted objects no longer keep the calling environment of their formulas,
  which could store unrelated workspace objects in saved fits.
- the RoBMA JAGS module is unregistered before its native objects are
  destroyed at process exit, preventing intermittent access violations on
  Windows; native kernels validate their counts, dimensions, and selection
  inputs before arithmetic and release their resources on R errors and
  interrupts; the build uses the C++ OpenMP compiler and linker flags and
  tracks all native header dependencies.

### Performance
- saved fits are smaller: through BayesTools 0.3.1, fits no longer carry
  runjags' compiled rjags model, which held another copy of the model and
  its data; `update()` recompiles the model from the stored chain states, as
  for a reloaded fit.
- post-fit native kernels (selected-normal likelihoods, normalizers, CDFs,
  moments, z-curve densities, and the EDR summary) evaluate their posterior
  rows in parallel on the thread budget of `RoBMA.options(native_threads =
  )` with identical values at any thread count, choosing the thread count
  from the work of a row; standard-normal tails inside selection kernels
  are evaluated in one vectorized pass (agreeing with `pnorm()` to 1e-13).
- normal multilevel cluster likelihoods use a native exact
  diagonal-plus-rank-one kernel shared by log-likelihood, LOO, and density
  computations, with a stable positive-sum identity for the quadratic form;
  bridge sampling of `cluster` models integrates the standardized study
  effects exactly.
- `funnel()`, `bfunnel()`, and `regplot()` contours use one native weighted
  posterior-mixture quantile engine (a machine-relative Brent solve for
  continuous mixtures and exact generalized-inverse bisection with atoms),
  and `zplot()` reuses invariant affine-density terms and evaluates fitted
  and extrapolated curves together.
- `add_loo()` builds the centered log-likelihood matrix once for relative
  effective sample sizes and PSIS, and influence diagnostics compute only
  the requested leave-one-out moments.
- repeated plot and summary calls reuse the fitted parameter-map cache, and
  qCMDE/IWMDE densities, selection-model post-processing, and
  `brma.mv()` likelihoods use exact batched, low-rank, sparse, and Markov
  covariance routes (with dense fallbacks) without changing their targets.

### Documentation
- documents the qCMDE/IWMDE estimating equations, the Savage-Dickey nesting
  requirement, method selection, tuning controls, reliability diagnostics,
  unsupported targets, the sensitivity workflow, and the limits of the
  normal approximation, and separates the posterior-row sampling error and
  the selected-row MCMC error; prior-ordinate and IWMDE proposal-weight
  estimation uncertainty remain outside `BF_error`.
- adds the vignettes "Posterior Densities and Parameter Tests" and
  "Multivariate and Multilevel Meta-Analysis".

## version 4.0.0
### Breaking changes
- rewrites the package around the unified `brma` class hierarchy. Single-model fits now use `brma()`, `brma.glmm()`, `bselmodel()`, `bPET()`, and `bPEESE()`; model-averaged fits use `BMA()`, `BMA.glmm()`, and `RoBMA()`.
- removes the legacy `RoBMA.reg()`, `NoBMA()`, `NoBMA.reg()`, `BiBMA()`, and `BiBMA.reg()` constructors. Use `mods`, `scale`, and `cluster` in the new constructors, `BMA()` for no-bias normal-likelihood model averaging, and `BMA.glmm()` for GLMM model averaging.
- replaces old input aliases such as `d`, `r`, `logOR`, `OR`, `z`, `y`, `se`, `v`, `n`, `study_names`, `study_ids`, `weight`, and `transformation` with `yi`, `vi`/`sei`, `ni`, `slab`, `cluster`, `weights`, `measure`, `output_measure`, and `transform`.
- removes legacy helper APIs including `combine_data()`, `check_setup()`, `extract_posterior()`, `marginal_summary()`, `marginal_plot()`, `plot_models()`, `adjusted_effect()`, `as_zcurve()`, and the old z-curve plotting methods.
- normal-likelihood fitting functions now require an explicit `measure` for fitted models. Use `measure = "GEN"` for generic effect sizes without a known unit-information scale.
- `update()` for `brma` objects now focuses on extending MCMC samples, updating labels, and refreshing cached quantities, not changing model structure.
- `set_convergence_checks()` no longer accepts the old `remove_failed` and `balance_probability` arguments.

### Features
- adds `brma()` / `brma.norm()` for single normal-likelihood Bayesian meta-analysis, including random-effects, meta-regression, multilevel, and location-scale models.
- adds `brma.glmm()` for binomial-normal and Poisson-normal GLMM meta-analysis from raw two-arm counts (`measure = "OR"` and `"IRR"`).
- adds single-model publication-bias constructors `bselmodel()`, `bPET()`, and `bPEESE()`.
- adds `BMA()` / `BMA.norm()` for Bayesian model averaging without publication-bias adjustment.
- adds `BMA.glmm()` for Bayesian model averaging of GLMM meta-analyses without publication-bias adjustment.
- rewrites `RoBMA()` as a product-space model-averaged ensemble over effect, heterogeneity, moderator, scale, and publication-bias components.
- adds formula/data-frame input handling for effect sizes, moderators, scale predictors, clusters, labels, subsets, likelihood weights, and raw GLMM counts.
- adds default prior construction from standardized effect-size measures, estimated or manually supplied unit-information standard deviations, and informed empirical priors.
- adds `prior_weightfunction()`, `wf_cumulative()`, `wf_fixed()`, and `wf_independent()` for BayesTools-backed selection-weightfunction priors.
- adds `prior_PET()`, `prior_PEESE()`, `prior_none()`, `prior_factor()`, `prior_informed()`, and BayesTools contrast helpers as package-level prior utilities.
- adds `posterior` package interfaces via `as_draws()`, `as_draws_array()`, `as_draws_df()`, `as_draws_list()`, `as_draws_matrix()`, and `as_draws_rvars()` for fitted models and `brma_samples`.
- adds the `brma_samples` posterior-sample class with print, summary, matrix, and `posterior` conversion methods.
- adds `predict.brma()` for posterior predictions of fixed terms, cluster effects, latent true effects, observed responses, and scale terms, with `newdata`, `conditional`, `bias_adjusted`, `output_measure`, and `transform` support.
- adds convenience wrappers `fitted()`, `pooled_effect()`, `pooled_heterogeneity()`, `blup()`, `true_effects()`, and `ranef()` for `brma` objects.
- adds model-comparison helpers `add_loo()`, `loo()`, `loo_compare()`, `loo_weights()`, `check_loo()`, `add_waic()`, and `waic()` using the `loo` package.
- adds bridge-sampling marginal likelihood support for single-model `brma` fits via `add_marglik()`, `bridge_sampler()`, `logml()`, `bf()`, `bayes_factor()`, and `post_prob()`.
- adds residual and influence diagnostics: `residuals()`, `rstandard()`, `rstudent()` / `LOO-PIT`, `hatvalues()`, `influence()`, `dfbetas()`, `dffits()`, `cooks.distance()`, `covratio()`, and `vif()`.
- adds plotting methods for `brma` objects: posterior/prior plots, `funnel()`, `regplot()`, `qqnorm()`, `radial()` / `galbraith()`, MCMC diagnostic plots, weightfunction plots, and PET-PEESE plots.
- adds `marginal_means()` with summary and plotting methods for moderator models.
- adds `summary_models()` for marginal and individual model-weight summaries of product-space `RoBMA`, `BMA`, and `BMA.glmm` objects.
- adds `interpret()` for concise textual interpretation of fitted `brma` and model-averaged objects.
- renames the zplot diagnostic API to `as_zplot()` and adds the direct plotting wrapper `zplot()`, with `plot()`, `hist()`, `lines()`, `summary()`, and print methods for zplot objects.
- adds `RoBMA.options()` and `RoBMA.get_option()` package options for defaults such as core count, automatic LOO/WAIC/marginal-likelihood computation, prior scaling defaults, and selection-bias defaults.

### Changes
- renames the multilevel clustering argument to `cluster`.
- renames study labels to `slab`, matching `metafor` naming.
- renames likelihood weights to `weights` for supported constructors and applies them consistently to posterior fitting, log-likelihoods, LOO, WAIC, and diagnostics; `brma.mv()` currently rejects likelihood weights.
- uses `measure`, `output_measure`, and `transform` for effect-size scale handling. Supported conversions include `SMD`, `COR`, `ZCOR`, and `OR`; `transform = "EXP"` exponentiates log ratio measures for display.
- standardizes continuous predictors by default and transforms reported coefficients back to the original scale unless standardized coefficients are requested.
- uses treatment contrasts by default for single-model constructors and mean-difference contrasts by default for model-averaged constructors.
- changes `predict.brma()` default to `type = "terms"`. GLMM `type = "response"` predictions return continuity-corrected effect-size estimators by default via `as_measure = TRUE`.
- separates output `unit` from `conditioning_depth` for residuals, fitted values, LOO, WAIC, and related diagnostics.
- supports estimate-level LOO/WAIC targets for `brma.mv()` and estimate-/cluster-level LOO/WAIC targets for multilevel `brma()` models, with target metadata to prevent invalid comparisons.
- keeps bridge-sampling marginal likelihoods for single-model `brma` objects; product-space `RoBMA`, `BMA`, and `BMA.glmm` objects rely on product-space only.
- routes selection-weightfunction priors through the BayesTools selection backend and selected-normal kernel, removing legacy weighted-normal mapping paths.
- uses `bias_indicator` and branch-aware selected-normal contexts for RoBMA publication-bias mixtures instead of inferring selection branches from `omega`.
- increases zplot default posterior thinning controls to `10000` samples and accepts `Inf` where full posterior evaluation is requested.
- adds `max_samples` controls to expensive funnel, regplot, and zplot summaries.
- updates the package startup message to point users to `vignette("v00-introduction", package = "RoBMA")`.
- requires BayesTools 0.3.0 for forward API and selection-backend support.
- adds `bridgesampling`, `loo`, `MASS`, and `parallel` as imports and `posterior` as a suggested package.

### Fixes
- fixes loading and runtime checks for the RoBMA JAGS module and native R routines.

### Performance and internals
- moves fitting to JAGS product-space models with mixture-prior indicators for model averaging.
- replaces legacy weighted-normal and multivariate-normal native code with selected-normal kernels shared by JAGS and R-native calls.
- adds native selected-normal routines for log likelihoods, normalizers, CDFs, moments, RNG, weighted summaries, funnel contours, regplot intervals, and zplot densities/threshold summaries.
- adds native estimate-unit GLMM marginal log-likelihood helpers for binomial and Poisson models.
- caches selected-normal normalizers and uses telescoping selection probabilities with log-space fallbacks for better numerical stability.
- relocates selected-normal C++ code to `src/selnorm/` and updates `Makevars*`, native registration, cleanup rules, and JAGS distribution registration.
- removes unused native matrix/LAPACK helper sources and older source-level transformation helpers.

### Documentation and tests
- reorganizes vignettes into numbered workflows covering introduction, prior distributions, baseline Bayesian meta-analysis, feature coverage, metafor parity, model averaging, RoBMA, multilevel models, medicine examples, and zplot diagnostics.
- regenerates roxygen documentation for the new constructors, priors, predictions, summaries, diagnostics, plots, model-comparison methods, and datasets.
- refreshes the README and pkgdown site for the 4.0.0 API.
- adds cached model fits under numbered vignette/model directories.
- refactors tests into ordered input, fitting, prediction, plotting, diagnostics, model-comparison, selected-normal kernel, and vignette-cache coverage.
- adds regression tests for selected-normal telescope probabilities, native/R fallback parity, posterior-row alignment, GLMM response conversion, LOO/WAIC targets, bridge sampling, and visual outputs.

## version 3.6.1
### Features
- `Explanation` vignette that helps navigate users through the vignettes
- two vignettes demonstrating robust Bayesian meta-analysis and meta-regressions
- `summary()` function now provides publication bias model type summary (`type = "models"`) for models fitted using `algorithm = "ss"`
- improves control over zplot diagnostics (i.e., specifying col, border, etc for the individual elements)

## version 3.6
### Features
- `funnel()` plot to visualize residuals vs the expected sampling distribution for `RoBMA()` and `RoBMA.reg()` models when using the `algorithm = "ss"`
- `residuals()` method for `RoBMA()` and `RoBMA.reg()` models when using the `algorithm = "ss"`
- `as_zplot()` function to transform meta-analytic models into a zplot object, only available for `RoBMA()` and `RoBMA.reg()` fitted using the `algorithm = "ss"`
- `plot()`, `summary()`, and `print()` functions for the `as_zplot` objects

## version 3.5.1
### Features
- `summary()` function now supports a `standardized_coefficients` argument to report either standardized (default) or raw meta-regression coefficients
- `extract()` function to extract the posterior samples of the model parameters
- `true_effects()` function to summarize the true effect size estimates of `RoBMA()` and `RoBMA.reg()` models when using the `algorithm = "ss"`
- `predict()` method for `RoBMA()` and `RoBMA.reg()` models when using the `algorithm = "ss"`

### Fixes
- fitting a meta-regression using predictors with missing values result in a clear error message

### Changes
- improving the speed of unit tests

## version 3.5
### Features
- approximate and computationally feasibly 3lvl selection models via the `RoBMA()` and `RoBMA.reg()` functions with the `cluster` argument when using `algorithm = "ss"`
- 3lvl binomial-normal models for binary data via the `BiBMA` and `BiBMA.reg` functions with the `cluster` argument when using `algorithm = "ss"`
- `pooled_effect()` function to compute the pooled effect size from the `RoBMA.reg`, `NoBMA.reg`, and `BiBMA.reg` models
- `adjusted_effect()` function to compute the adjusted effect size from the `RoBMA.reg`, `NoBMA.reg`, and `BiBMA.reg` models
- enables `summary_heterogeneity()` for BiBMA models

### Fixes
- passing and checks of the `cluster` and `study_labels` arguments
- PEESE prior distribution now scale as 1/scale instead of 1/scale^2 with the `rescale_priors` argument  
- the conditional prediction interval based on `summary_heterogeneity()` is now conditional on the presence of the effect
- additional minor prior handling fixes (i.e., missing marginal estimates when only alternative prior distributions were specified etc)
- diagnostics with mixture baseline priors when using `algorithm = "ss"`
- `summary_heterogeneity()` with only a single study does not produce relative heterogeneity instead of crashing

## version 3.4
### Features
- adding binomial-normal meta-regression models for binary data via the `BiBMA.reg` function
- the spike and slab algorithm for faster model estimation via the `algorithm = "ss"` argument for BiBMA models
- default prior distributions for all parameters of BiBMA models are now set via the `set_default_binomial_priors()` function

## version 3.3
### Features
- the spike and slab algorithm for faster model estimation via the `algorithm = "ss"` argument (see a new vignette for more details)
- refactoring of the JAGS C++ code of weighted distributions and exporting of the lpdfs into JAGS (maintenance)
- weights_mix JAGS prior distribution to sample a mixture of weight functions directly

### Fixes
- incorrectly omitting models with more than one predictor when computing conditional marginal summary

## version 3.2.1
### Features
- default prior distributions for all parameters are now set via the `set_default_priors()` function
- `rescale_priors` argument allows to conveniently re-scale the prior distributions for the effect, heterogeneity, and bias simultaneously

## version 3.2
### Features
- `summary_heterogeneity()` function to summarize the heterogeneity of the RoBMA models (prediction interval, tau, tau^2, I^2, and H^2)
- `check_RoBMA_convergence()` function to check the convergence of the RoBMA models
- adds informed prior distributions for binary and time-to-event outcomes via BayesTools 0.2.17

### Fixes
- checking and fixing the number of available cores upon loading the package (hopefully fixes some parallelization issues)
- `update()` function re-evaluates convergence checks of individual models (https://github.com/FBartos/RoBMA/issues/34) 
- typos and minor issues in the vignettes


## version 3.1
### Features
- binomial-normal models for binary data via the `BiBMA` function
- `NoBMA` and `NoBMA.reg()` functions as wrappers around `RoBMA` `RoBMA.reg()` functions for simpler specification of publication bias unadjusted Bayesian model-averaged meta-analysis
- adding odds ratios output transformation` 
- extending (instead of a complete refitting) of models via the `update.RoBMA()` function (only non-converged models by default or all by setting `extend_all = TRUE`)

### Fixes
- handling of non-converged models

## version 3.0.1
### Fixes (thanks to Don & Rens)
- compilation issues with Clang (https://github.com/FBartos/RoBMA/issues/28)
- lapack path specifications (https://github.com/FBartos/RoBMA/issues/24)

## version 3.0
### Features
- meta-regression with `RoBMA.reg()` function
- posterior marginal summary and plots for the `RoBMA.reg` models with `summary_marginal()` and `plot_marginal()` functions
- new vignette on hierarchical Bayesian model-averaged meta-analysis
- new vignette on robust Bayesian model-averaged meta-regression
- adding vignette from AMPPS tutorial
- faster implementation of JAGS multivariate normal distribution (based on the BUGS JAGS module)
- incorporating `weight` argument in the `RoBMA` and `combine_data` functions in order to pass `custom` likelihood weights
- ability to use inverse square weights in the weighted meta-analysis by setting a `weighted_type = "inverse_sqrt"` argument 

### Changes
- reworked interface for the hierarchical models. Prior distributions are now specified via the `priors_hierarchical` and `priors_hierarchical_null` arguments instead of `priors_rho` and `priors_rho_null`. The model summary now shows `Hierarchical` component summary.

## version 2.3.2
### Fixes
- suppressing start-up message 
- cleaning up imports

## version 2.3.1
### Fixes
- fixing weighted meta-analysis parameterization 

## version 2.3
### Features
- weighted meta-analysis by specifying `cluster` argument in `RoBMA()` and setting `weighted = TRUE`. The likelihood contribution of estimates from each study is down-weighted proportionally to the number of estimates in that study. Note that this experimental feature is supposed to provide a conservative alternative for estimating RoBMA in cases with multiple estimates from a study where the multivariate option is not computationally feasible.

## version 2.2.3
### Fixes
- updating the Makevars to install with R 4.2 and JAGS 4.3.1

## version 2.2.2
### Fixes
- updating the C++ to compile on M1 Mac

## version 2.2.1
### Changes
- message about the effect size scale of parameter estimates is always shown
- compatibility with BayesTools 0.2.0+

## version 2.2
### Features
- three-level meta-analysis by specifying `cluster` argument in `RoBMA`. However, note that this is (1) an experimental feature and (2) the computational expense of fitting selection models with clustering is extreme. As of now, it is almost impossible to have more than 2-3 estimates clustered within a single study).

## version 2.1.2
### Fixes
- adding Windows ucrt patch (thanks to Tomas Kalibera)
- adding BayesTools version check

## version 2.1.1
### Fixes
- incorrectly formatted citations in vignettes and capitalization

### Features
- adding `informed_prior()` function (from the BayesTools package) that allows specification of various informed prior distributions from the field of medicine and psychology
- adding a vignette reproducing the example of dentine sensitivity with the informed Bayesian model-averaged meta-analysis from Bartoš et al., 2021 ([open-access](https://onlinelibrary.wiley.com/doi/10.1002/sim.9170)),
- further reductions of fitted object size when setting `save = "min"`

## version 2.1
### Fixes
- more informative error message when the JAGS module fails to load
- correcting wrong PEESE transformation for the individual models summaries (issue #12)
- fixing error message for missing conditional PET-PEESE
- fixing incorrect lower bound check for log(OR)

### Features
- adding `interpret()` function (issue #11)
- adding effect size transformation via `output_scale` argument to `plot()` and `plot_models()` functions
- better handling of effect size transformations and scaling - BayesTools style back-end functions with Jacobian transformations

## version 2.0
Please notice that this is a major release that breaks backwards compatibility.

### Changes
 - naming of the arguments specifying prior distributions for the different parameters/components of the models changed (`priors_mu` -> `priors_effect`, `priors_tau` -> `priors_heterogeneity`, and `priors_omega` -> `priors_bias`),
 - prior distributions for specifying weight functions now use a dedicated function (`prior(distribution = "two.sided", parameters = ...)` -> `prior_weightfunction(distribution = "two.sided", parameters = ...)`),
 - new dedicated function for specifying no publication bias adjustment component / no heterogeneity component (`prior_none()`),
 - new dedicated functions for specifying models with the PET and PEESE publication bias adjustments (`prior_PET(distribution = "Cauchy", parameters = ...)` and `prior_PEESE(distribution = "Cauchy", parameters = ...)`),
 - new default prior distribution specification for the publication bias adjustment part of the models (corresponding to the RoBMA-PSMA model from Bartoš et al., 2021 [manuscript](https://doi.org/10.1002/jrsm.1594)),
 - new `model_type` argument allowing to specify different "pre-canned" models (`"PSMA"` = RoBMA-PSMA, `"PP"` = RoBMA-PP, `"2w"` = corresponding to Maier et al., in press , [manuscript](https://doi.org/10.1037/met0000405)),
 - `combine_data` function allows combination of different effect sizes / variability measures into a common effect size measure (also used from within the `RoBMA` function),
 - better and improved automatic fitting procedure now enabled by default (can be turned of with `autofit = FALSE`)
 - prior distributions can be specified on the different scale than the supplied effect sizes (the package fits the model on Fisher's z scale and back transforms the results back to the scale that was used for prior distributions specification, Cohen's d by default, but both of them can be overwritten with the `prior_scale` and `transformation` arguments),
 - new prior distributions, e.g., beta or fixed weight functions,
 - estimates from individual models are now plotted with the `plot_models()` function and the forest plot can be obtained with the `forest()` function,
 - the posterior distribution plots for the individual weights are no able supported, however, the weightfunction and the PET-PEESE publication bias adjustments can be visualized with the `plot.RoBMA()` function and `parameter = "weightfunction"` and `parameter = "PET-PEESE"`.

## version 1.2.1
### Fixes
- check_setup function not working at all

## version 1.2.0
### Changes
- the studies's true effects are now marginalized out of the random effects models and are no longer estimated (see Appendix A of our [manuscript](https://doi.org/10.1037/met0000405) for more details). As a results, arguments referring to the true effects are now disabled.
- all models are now being estimated using the likelihood of effect sizes (instead of test-statistics as usually defined). We reproduced the simulation study that we used to evaluate the method performance and it achieved identical results (up to MCMC error, before marginalizing out the true effects). A big advantage of using the normal likelihood for effect sizes is a considerable speed up of the whole estimation process.
- as a results of these two changes, the results of the models will differ to those of pre 1.2.0 version

### Fixes
- autofit being turn on if any control argument was specified

## version 1.1.2
### Fixes
- vdiffr not being used conditionally in unit tests

## version 1.1.1
### Fixes
- inability to fit a model without specifying a seed
- inability to produce individual model plots due to incompatibility with the newer versions of ggplot2  

## version 1.1.0
### Features
- parallel within and between model fitting using the parallel package with 'parallel = TRUE' argument

## version 1.0.5
### Fixes:
- models being fitted automatically until reaching R-hat lower than 1.05 without setting max_rhat and autofit control parameters
- bug preventing to draw a bivariate plot of mu and tau
- range for parameter estimates from individual models no containing 0 (or 1 in case of OR measured effect sizes)
- inability to fit a model with only null mu distributions if correlation or OR measured effect sizes were specified
- ordering of the estimated and observed effects when both of them are requested simultaneously
- formatting of this file (NEWS.md)

### Improvements:
- priors plot: parameter specification, default plotting range, clearer x-axis labels in cases when the parameter is defined on transformed scale
- parameters plots: probability scale always ends at the same spot as is the last tick on the density scale
- adding warnings if any of the specified models has Rhat higher than 1.05 or the specified value
- grouping the same warnings messages together

## version 1.0.4
### Fixes:
- inability to run models without the silent = TRUE control

## version 1.0.3
### Features:
- x-axis rescaling for the weight function plot (by setting 'rescale_x = TRUE' in the 'plot.RoBMA' function)
- setting expected direction of the effect in for RoBMA function

### Fixes:
- marginal likelihood calculation for models with spike prior distribution on mean parameter which location was not set to 0
- some additional error messages 

### CRAM requested changes:
- changing information messages from 'cat' to 'message' from plot related functions
- saving and returning the 'par' settings to the user defined one in the base plot functions

## version 1.0.2
### Fixes:
- the summary and plot function now shows quantile based confidence intervals for individual models instead of the HPD provided before (this affects only 'summary'/'plot' with 'type = "individual"', all other confidence intervals were quantile based before)

## version 1.0.1
### Fixes:
- summary function returning median instead of mean

## version 1.0.0 (vs the osf version)
### Fixes:
- incorrectly weighted theta estimates
- models with non-zero point prior distribution incorrectly plotted using when "models" option in case that the mu parameter was transformed

### Additional features:
- analyzing OR
- distributions implemented using boost library (helps with convergence issues)
- ability to mute the non-suppressible "precision not achieved" warning messages by using "silent" = TRUE inside of the control argument
- vignettes

### Notable changes:
- the way how the seed is set before model fitting (the simulation study will not be reproducible with the new version of the package)
