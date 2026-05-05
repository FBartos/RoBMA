# Bayesian Meta-Analysis

Function for fitting normal-likelihood/effect-size random-effects,
meta-regression, multilevel, and location-scale meta-analytic models.
`brma.norm()` is an alias for `brma()`; raw-count GLMM models use
[`brma.glmm()`](https://fbartos.github.io/RoBMA/reference/brma.glmm.md).

## Usage

``` r
brma(
  yi,
  vi,
  sei,
  weights,
  ni,
  mods,
  scale,
  cluster,
  data,
  slab,
  subset,
  measure,
  prior_effect,
  prior_heterogeneity,
  prior_mods,
  prior_scale,
  prior_heterogeneity_allocation,
  standardize_continuous_predictors = TRUE,
  set_contrast_factor_predictors = "treatment",
  prior_unit_information_sd,
  rescale_priors = 1,
  prior_informed_field,
  prior_informed_subfield,
  sample = 5000,
  burnin = 2000,
  adapt = 500,
  chains = 3,
  thin = 1,
  parallel = FALSE,
  autofit = FALSE,
  autofit_control = set_autofit_control(),
  convergence_checks = set_convergence_checks(),
  seed = NULL,
  silent,
  ...
)
```

## Arguments

- yi:

  a vector of effect sizes, or a formula with the effect size on the
  left-hand side and location moderators on the right-hand side (for
  example `yi ~ x1 + x2`). If a formula is supplied, `mods` must not be
  specified.

- vi:

  a vector of sampling variances. Either `vi` or `sei` must be supplied
  for normal models.

- sei:

  a vector of standard errors. Either `vi` or `sei` must be supplied for
  normal models.

- weights:

  an optional vector of positive likelihood weights. For
  normal/effect-size models, each weight powers the estimate likelihood.
  For constructors with GLMM raw-count input, each weight powers the
  paired two-arm likelihood for one study.

- ni:

  an optional vector of sample sizes. Used for `measure = "GEN"` or when
  estimating `"UISD"`).

- mods:

  an optional matrix, data.frame, or formula specifying location
  moderators (meta-regressors). Formula input is evaluated in `data`.

- scale:

  an optional matrix, data.frame, or formula specifying scale predictors
  for location-scale models. Formula input is evaluated in `data`.

- cluster:

  an optional vector of cluster identifiers for multilevel
  meta-analysis.

- data:

  an optional data frame containing the variables.

- slab:

  an optional vector of study labels.

- subset:

  an optional logical or numeric vector specifying a subset of data to
  be used.

- measure:

  a character string specifying the effect size measure.
  Normal/effect-size constructors require an explicit value and support
  `"SMD"`, `"ZCOR"`, `"RR"`, `"OR"`, `"HR"`, `"RD"`, `"IRR"`, and
  `"GEN"`. Use `"GEN"` only for general effect sizes without a known
  unit information standard deviation. GLMM raw-count constructors
  support only `"OR"` and `"IRR"` and default to `"OR"`.

- prior_effect:

  prior distribution for the effect size (\\\mu\\) parameter (the
  intercept). If omitted, a default prior is constructed. In
  single-model functions, explicit `NULL` or `FALSE` sets a spike at
  zero.

- prior_heterogeneity:

  prior distribution for the heterogeneity (\\\tau\\) parameter. If
  omitted, a default prior is constructed. In single-model functions,
  explicit `NULL` or `FALSE` sets a spike at zero.

- prior_mods:

  prior distribution for the moderators (\\\beta\\) parameters. A single
  prior applies to all terms; a named list can specify term-specific
  priors. If omitted or `NULL`, default priors are used.

- prior_scale:

  prior distribution for the scale (\\\delta\\) parameters. A single
  prior applies to all terms; a named list can specify term-specific
  priors. If omitted or `NULL`, default priors are used.

- prior_heterogeneity_allocation:

  prior distribution for the fraction of heterogeneity allocated to the
  cluster-level component in multilevel models (\\\rho\\). If omitted or
  `NULL`, defaults to `Beta(1, 1)`.

- standardize_continuous_predictors:

  logical. Whether to standardize continuous predictors. Defaults to
  `TRUE`.

- set_contrast_factor_predictors:

  character. How to set contrast for factor predictors. Defaults are
  constructor-specific and shown in each function usage; single-model
  constructors use `"treatment"`, while model-averaging constructors use
  `"meandif"`.

- prior_unit_information_sd:

  numeric. The unit information standard deviation (\\\sigma\_{unit}\\).
  Cannot be used together with `prior_informed_field`.

- rescale_priors:

  numeric. A scaling factor for supported prior distributions. Point and
  none priors are unchanged; publication-bias priors are not rescaled
  except for the default PEESE prior's UISD adjustment. Defaults to 1.

- prior_informed_field:

  character. The field of the informed prior distributions. Omit to use
  the standard default prior specification; explicit `NULL` is invalid.

- prior_informed_subfield:

  character. The subfield of the informed prior distributions. Omit to
  use the field-specific default, such as `"Cochrane"` for
  `prior_informed_field = "medicine"`; explicit `NULL` is invalid.

- sample:

  numeric. Number of MCMC samples to save. Defaults to `5000`.

- burnin:

  numeric. Number of burn-in iterations. Defaults to `2000`.

- adapt:

  numeric. Number of adaptation iterations. Defaults to `500`.

- chains:

  numeric. Number of MCMC chains. Defaults to `3`.

- thin:

  numeric. Thinning interval. Defaults to `1`.

- parallel:

  logical. Whether to run MCMC chains in parallel. Defaults to `FALSE`.

- autofit:

  logical. Whether to automatically extend the MCMC chains if
  convergence is not met. Defaults to `FALSE`.

- autofit_control:

  list of autofit control settings. See
  [`set_autofit_control()`](https://fbartos.github.io/RoBMA/reference/RoBMA_control.md)
  for details.

- convergence_checks:

  list of convergence check settings. See
  [`set_convergence_checks()`](https://fbartos.github.io/RoBMA/reference/RoBMA_control.md)
  for details.

- seed:

  numeric. Random seed for reproducibility. Defaults to `NULL`.

- silent:

  logical. Whether to suppress output. Constructors with no explicit
  default use `RoBMA.get_option("silent")` when `silent` is omitted.
  Model-averaging wrappers default to `TRUE` unless explicitly changed.

- ...:

  additional advanced arguments. Fitting functions reject unused
  arguments; currently recognized internal arguments include
  `only_data`, `only_priors`, `is_JASP`, and `is_JASP_prefix`.

## Value

A fitted object of class `c("brma.norm", "brma")`. The object contains
checked `data`, checked `priors`, the JAGS `fit`, cached `summary`, and
cached `coefficients`. If the corresponding package options are enabled,
it can also contain cached LOO, WAIC, or marginal likelihood results.
The advanced internal `only_data = TRUE` and `only_priors = TRUE` paths
return partially constructed objects.

## Details

### Prior distributions

Prior distributions must be specified for all model parameters. This
typically includes the pooled effect \\\mu\\ and between-study
heterogeneity \\\tau\\. In the case of meta-regression, the pooled
effect \\\mu\\ corresponds to the intercept, and additional prior
distributions for the regression coefficients are required. In the case
of a location-scale model, the between-study heterogeneity corresponds
to the intercept of the scale regression, and additional prior
distributions for the scale regression coefficients are required.

There are several ways to specify the prior distributions:

1.  via a standardized effect size `measure` with known unit information
    standard deviation,

2.  by estimating unit information standard deviation using sample sizes
    `ni`,

3.  by manually setting `prior_unit_information_sd`,

4.  by specifying informed empirical prior distributions via
    `prior_informed_field` and `prior_informed_subfield`,

5.  or via fully custom specification using the `prior_effect`,
    `prior_heterogeneity`, `prior_mods`, `prior_scale`, and
    `prior_heterogeneity_allocation` arguments.

In all cases, the prior behavior can be further modified by the
`rescale_priors`, `standardize_continuous_predictors`, and
`set_contrast_factor_predictors` arguments.

#### (1) Specifying prior distributions via standardized effect size measures with known

unit information standard deviation This is the easiest way to specify
prior distributions. The width of prior distributions is based on a
fraction of the known unit information standard deviation (`UISD`)
(Röver et al. 2021) . The default prior distributions for the parameters
are set as follows:

|                           |                                    |
|---------------------------|------------------------------------|
| effect size:              | Normal(0, \\\frac{1}{2}\\ `UISD`)  |
| heterogeneity:            | Normal+(0, \\\frac{1}{4}\\ `UISD`) |
| effect moderation:        | Normal(0, \\\frac{1}{4}\\ `UISD`)  |
| heterogeneity moderation: | Normal(0, \\\frac{1}{2}\\)         |

The heterogeneity moderation parameters are multiplicative, as such they
are independent of `UISD`.

The default fraction of the `UISD` can be changed via
[`RoBMA.options()`](https://fbartos.github.io/RoBMA/reference/RoBMA_options.md)
using one of the following arguments: `"default_UISD.effect"`,
`"default_UISD.heterogeneity"`, `"default_UISD.mods"`,
`"default_UISD.scale"`.

The known `UISD` for standardized effect size measures (`measure`) are
set as follows:

|           |              |
|-----------|--------------|
| `"SMD"`:  | \\\sqrt{2}\\ |
| `"ZCOR"`: | \\1\\        |
| `"RR"`:   | \\\sqrt{4}\\ |
| `"OR"`:   | \\\sqrt{4}\\ |
| `"HR"`:   | \\\sqrt{4}\\ |
| `"IRR"`:  | \\\sqrt{4}\\ |

See Chapter 2.4 in Spiegelhalter et al. (2004) and Chapter 1 in Grieve
(2022) .

#### (2) Estimating unit information standard deviation using sample sizes

When effect sizes are on a non-standardized scale (`measure = "GEN"`) or
use a standardized effect size without known `UISD`, the `UISD` can be
estimated from sample sizes (`ni`) and standard errors (`sei`) following
Equation 6 in Röver et al. (2021) . The estimated `UISD` is then used to
scale the default prior distributions as described in section (1).

Note that the known `UISD` for standardized effect size measures
(section (1)) is used if available, even when `ni` is provided.

#### (3) Manually setting unit information standard deviation

Alternatively, the `UISD` can be specified directly via the
`prior_unit_information_sd` argument. This is useful when the
appropriate scale for prior distributions is known a priori or when
multiple analyses are to be performed on subsets of the same data
(re-estimating `UISD` on different subsets of the data can lead to
slightly different prior distributions for different subsets; see
[`estimate_unit_information_sd()`](https://fbartos.github.io/RoBMA/reference/estimate_unit_information_sd.md)).
The specified `prior_unit_information_sd` is then used to scale the
default prior distributions as described in section (1).

Note that the manually specified `prior_unit_information_sd` takes
precedence over the estimated `UISD` from `ni` (section (2)) and the
known `UISD` from `measure` (section (1)). It cannot be combined with
`prior_informed_field`.

#### (4) Specifying informed empirical prior distributions

Informed prior distributions can be specified via the
`prior_informed_field` and `prior_informed_subfield` arguments.
Currently, only `prior_informed_field = "medicine"` with subfields
defined in
[`BayesTools::prior_informed_medicine_names`](https://fbartos.github.io/BayesTools/reference/prior_informed_medicine_names.html)
is supported, which uses empirically derived prior distributions from
medical meta-analyses as described in Bartoš et al. (2021) and Bartoš et
al. (2023) .

When `prior_informed_field = "medicine"`, the default
`prior_informed_subfield` is `"Cochrane"` (i.e., using the whole CDSR
database as a reference). The informed prior distributions are available
for the following effect size measures:

|          |                              |
|----------|------------------------------|
| `"SMD"`: | standardized mean difference |
| `"OR"`:  | log odds ratio               |
| `"RR"`:  | log risk ratio               |
| `"RD"`:  | risk difference              |
| `"HR"`:  | log hazard ratio             |

Note that informed prior distributions are only available for the effect
size (\\\mu\\) and heterogeneity (\\\tau\\) parameters. For effect
moderation (`prior_mods`), the informed effect prior is scaled by a
factor of `RoBMA.get_option("default_informed_priors.mods")`. For
heterogeneity moderation (`prior_scale`), a normal prior with standard
deviation specified by
`RoBMA.get_option("default_informed_priors.scale")` is used.

#### (5) Fully custom prior distributions

Prior distributions can be fully customized by directly specifying the
`prior_effect`, `prior_heterogeneity`, `prior_mods`, `prior_scale`, and
`prior_heterogeneity_allocation` arguments. These should be prior
distribution objects created via
[`BayesTools::prior()`](https://fbartos.github.io/BayesTools/reference/prior.html)
or related functions (e.g.,
[`prior_factor()`](https://fbartos.github.io/RoBMA/reference/prior_factor.md)).

#### Rescaling prior distributions

The `rescale_priors` argument allows rescaling supported prior
distributions by a multiplicative factor. For example,
`rescale_priors = 2` doubles the standard deviations/scales of normal,
Cauchy, t, and inverse-gamma prior distributions, making them more
diffuse. Point and none priors are unchanged.

#### Handling of continuous and factor predictors

When `standardize_continuous_predictors = TRUE`, continuous moderator
and scale predictors are internally standardized before fitting.
Reported summaries are transformed back to the original predictor scale
by default; use `standardized_coefficients = TRUE` in
[`summary()`](https://rdrr.io/r/base/summary.html) or related methods to
inspect coefficients on the standardized scale.

Factor predictors use the contrast specified by
`set_contrast_factor_predictors`. The default `"treatment"` follows
standard treatment coding. Model-averaging functions default to
`"meandif"` so inclusion priors for factor levels are centered on
deviations from the grand mean.

## References

Bartoš F, Gronau QF, Timmers B, Otte WM, Ly A, Wagenmakers E (2021).
“Bayesian model-averaged meta-analysis in medicine.” *Statistics in
Medicine*, **40**(30), 6743–6761.
[doi:10.1002/sim.9170](https://doi.org/10.1002/sim.9170) .  
  
Bartoš F, Otte WM, Gronau QF, Timmers B, Ly A, Wagenmakers E (2023).
“Empirical prior distributions for Bayesian meta-analyses of binary and
time-to-event outcomes.”
[doi:10.48550/arXiv.2306.11468](https://doi.org/10.48550/arXiv.2306.11468)
. Preprint available at https://doi.org/10.48550/arXiv.2306.11468.  
  
Grieve AP (2022). *Hybrid frequentist/Bayesian power and Bayesian power
in planning clinical trials*. Chapman and Hall/CRC.  
  
Röver C, Bender R, Dias S, Schmid CH, Schmidli H, Sturtz S, Weber S,
Friede T (2021). “On weakly informative prior distributions for the
heterogeneity parameter in Bayesian random-effects meta-analysis.”
*Research Synthesis Methods*, **12**(4), 448–474.
[doi:10.1002/jrsm.1475](https://doi.org/10.1002/jrsm.1475) .  
  
Spiegelhalter DJ, Abrams KR, Myles JP (2004). *Bayesian approaches to
clinical trials and health-care evaluation*. John Wiley and Sons.
[doi:10.1002/0470092602](https://doi.org/10.1002/0470092602) .

## See also

[`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md),
[`BMA()`](https://fbartos.github.io/RoBMA/reference/BMA.md),
[`brma.glmm()`](https://fbartos.github.io/RoBMA/reference/brma.glmm.md),
[`summary.brma()`](https://fbartos.github.io/RoBMA/reference/summary.brma.md),
[`plot.brma()`](https://fbartos.github.io/RoBMA/reference/plot.brma.md),
[`predict.brma()`](https://fbartos.github.io/RoBMA/reference/predict.brma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE) &&
    requireNamespace("metafor", quietly = TRUE)) {
  data(dat.bcg, package = "metadat")
  dat <- metafor::escalc(
    measure = "RR",
    ai      = tpos,
    bi      = tneg,
    ci      = cpos,
    di      = cneg,
    data    = dat.bcg
  )

  fit <- brma(
    yi      = yi,
    vi      = vi,
    mods    = ~ ablat + year,
    data    = dat,
    measure = "RR",
    seed    = 1,
    silent  = TRUE
  )

  summary(fit)
  predict(fit, type = "terms")
}
} # }
```
