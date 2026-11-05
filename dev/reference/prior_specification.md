# Prior specification

Prior distributions can be supplied for all model parameters. When
omitted, the fitting functions construct defaults from the effect-size
measure, sample sizes, or informed-prior settings when available.

This typically includes the pooled effect \\\mu\\ and between-study
heterogeneity \\\tau\\. In the case of meta-regression, the pooled
effect \\\mu\\ corresponds to the intercept, and additional prior
distributions for the regression coefficients are required. In the case
of a location-scale model, the between-study heterogeneity corresponds
to the intercept of the scale regression, and additional prior
distributions for the scale regression coefficients are required.

Priors of parameters with a restricted support are truncated to it:
priors of the heterogeneity and of scale intercepts to \\\[0, \infty)\\,
priors of the heterogeneity allocation and of the base rate to \\\[0,
1\]\\. An unbounded truncation point, such as the default of
`prior("normal", ...)`, is moved to the support boundary with a warning;
a finite truncation point outside the support is an error.

## Arguments

- prior_effect:

  prior distribution for the effect size (\\\mu\\) parameter (the
  intercept). If omitted and consumed, a default prior is constructed.
  In single-model functions, explicit `NULL` or `FALSE` sets a spike at
  zero.

- prior_heterogeneity:

  prior distribution for the heterogeneity (\\\tau\\) parameter. If
  omitted and consumed, a default prior is constructed. In single-model
  functions, explicit `NULL` or `FALSE` sets a spike at zero. For
  [`brma.mv`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
  random-effect formulas, this can also be a partial or complete
  [`prior_random`](https://fbartos.github.io/RoBMA/dev/reference/prior_random.md)
  object; see
  [`random_effect_prior_specification`](https://fbartos.github.io/RoBMA/dev/reference/random_effect_prior_specification.md).

- prior_mods:

  prior distribution for the moderators (\\\beta\\) parameters. A single
  prior applies to all terms; a named list can specify term-specific
  priors. If omitted or `NULL`, default priors are used. This argument
  can be specified only when `mods` is specified. A formula without an
  intercept (`~ 0 + x`) fixes the location intercept at zero, so an
  `intercept` entry of the list is accepted only as a point prior at
  zero; any other `intercept` prior stops with an error of class
  `RoBMA_prior_intercept_conflict` (also of class `RoBMA_prior_error`).
  A named `intercept` entry supplies the complete intercept prior.
  Explicit competing global effect declarations are refused when this
  would discard them. Named `NULL` or `FALSE` intercept entries are
  invalid; specify the complete prior in one place.

- prior_scale:

  prior distribution for the scale (\\\delta\\) parameters. A single
  prior applies to all terms; a named list can specify term-specific
  priors. If omitted or `NULL`, default priors are used. This argument
  can be specified only when `scale` is specified. A named `intercept`
  entry supplies the complete positive SD prior used by
  `log(intercept)`. Global heterogeneity declarations remain meaningful
  when they supply another random root or its model-inclusion weights;
  declarations fully replaced without a consumer raise
  `RoBMA_prior_intercept_conflict`.

- prior_heterogeneity_allocation:

  prior distribution for the fraction of heterogeneity allocated to the
  cluster-level component in multilevel models (\\\rho\\, a variance
  fraction rather than a random-coefficient correlation). If omitted or
  `NULL`, defaults to `Beta(1, 1)`. This argument can be specified only
  when `cluster` is specified.

- prior_baserate:

  prior distribution for the estimate-specific midpoint base-rate
  probability in binomial GLMM models (`measure = "OR"`). If omitted or
  `NULL`, defaults to independent `Beta(1, 1)` priors. Continuous priors
  are restricted to `[0, 1]` (see `prior_specification`) and point
  priors must lie strictly inside `(0, 1)`. The estimate-unit
  log-likelihood of the model
  ([`log_lik()`](https://fbartos.github.io/RoBMA/dev/reference/log_lik.md),
  LOO, WAIC) integrates the base rate over any of these priors: point
  priors and untruncated beta priors adaptively (a uniform(0, 1) prior
  is a Beta(1, 1) prior), other families on a refined prior-CDF grid.
  The grid is accepted only when its log-likelihood changes by less than
  `1e-6` between refinements; any prior without an adaptive route (other
  than an untruncated beta prior; a truncated beta prior takes the grid
  too) can fail to reach this, especially with few events or very large
  arms, and then stops with an error of class
  `RoBMA_glmm_grid_nonconvergence`, as a lognormal base-rate prior does
  on cells with few events. A beta prior avoids it. The automatic
  computation of LOO/WAIC after fitting warns instead and returns the
  fit without them.

- prior_lograte:

  prior distribution for the estimate-specific midpoint log-rate in
  Poisson GLMM models (`measure = "IRR"`). If omitted or `NULL`, a
  normal prior centered at the pooled crude log-rate is used
  independently for each estimate. Its standard deviation defaults to
  `1` and can be changed via `RoBMA.options(default_lograte.sd = ...)`.
  Rescaling all person-times shifts the prior mean by the corresponding
  log conversion while leaving its standard deviation unchanged. The
  estimate-unit log-likelihood of the model
  ([`log_lik()`](https://fbartos.github.io/RoBMA/dev/reference/log_lik.md),
  LOO, WAIC) integrates the log-rate over the prior: point priors and
  untruncated normal priors adaptively, other families on a refined
  prior-CDF grid. The grid is accepted only when its log-likelihood
  changes by less than `1e-6` between refinements; any prior without an
  adaptive route (other than an untruncated normal prior; a truncated
  normal prior takes the grid too) can fail to reach this, especially
  with few events or very large arms, and then stops with an error of
  class `RoBMA_glmm_grid_nonconvergence`, as gamma, lognormal, and
  inverse-gamma log-rate priors do on cells with few events. A normal
  prior avoids it. The automatic computation of LOO/WAIC after fitting
  warns instead and returns the fit without them.

- prior_unit_information_sd:

  numeric. The unit information standard deviation (\\\sigma\_{unit}\\).
  Cannot be used together with `prior_informed_field`.

- rescale_priors:

  numeric. A scaling factor for supported prior distributions. Point and
  none priors are unchanged; supplied
  [`prior_random()`](https://fbartos.github.io/RoBMA/dev/reference/prior_random.md)
  scale architectures and publication-bias priors are not rescaled. The
  default PEESE prior separately uses its documented UISD adjustment,
  independently of `rescale_priors`. Defaults to 1.

- standardize_continuous_predictors:

  logical. Whether to standardize continuous predictors. Defaults to
  `TRUE`.

- set_contrast_factor_predictors:

  character. How to set contrast for factor predictors. Defaults are
  constructor-specific and shown in each function usage. Single-model
  constructors use `"treatment"`; model-averaging constructors use
  `"meandif"`. Accepted values are `"treatment"`, `"meandif"`,
  `"orthonormal"`, and `"independent"`. Independent coding creates one
  coefficient per factor level. Random-coefficient blocks reuse the
  resolved fixed-factor contrast unless `random_block(contrasts = ...)`
  overrides it. A warning is issued when default or inherited
  independent contrasts make a fixed or random design overspecified.
  With `"meandif"` or `"orthonormal"` contrasts, a formula term that
  interacts a factor with a predictor whose own term is missing, such as
  `g:x` in `~ g + g:x` or `~ g / x`, stops the fit: the term codes `g`
  by level indicators, with one coefficient per level, on which
  mean-difference and orthonormal priors are not defined. Specify an
  independent prior for that term, e.g.,
  `prior_mods = list("g:x" = prior_factor("normal", list(0, 1), contrast = "independent"))`
  (or `prior_scale` for a scale formula), include the missing term
  (`~ g * x`), or use `set_contrast_factor_predictors = "treatment"`.

- prior_informed_field:

  character. The field of the informed prior distributions. Omit to use
  the standard default prior specification; explicit `NULL` is invalid.

- prior_informed_subfield:

  character. The subfield of the informed prior distributions. Can only
  be specified together with `prior_informed_field`. Omit to use the
  field-specific default, such as `"Cochrane"` for
  `prior_informed_field = "medicine"`; explicit `NULL` is invalid.

## Details

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

### (1) Specifying prior distributions via standardized effect size measures with known

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
[`RoBMA.options()`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA_options.md)
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

### (2) Estimating unit information standard deviation using sample sizes

When effect sizes are on a non-standardized scale (`measure = "GEN"`) or
use a standardized effect size without known `UISD`, the `UISD` can be
estimated from sample sizes (`ni`) and standard errors (`sei`) following
Equation 6 in Röver et al. (2021) . The estimated `UISD` is then used to
scale the default prior distributions as described in section (1).

Note that the known `UISD` for standardized effect size measures
(section (1)) is used if available, even when `ni` is provided.

### (3) Manually setting unit information standard deviation

Alternatively, the `UISD` can be specified directly via the
`prior_unit_information_sd` argument. This is useful when the
appropriate scale for prior distributions is known a priori or when
multiple analyses are to be performed on subsets of the same data
(re-estimating `UISD` on different subsets of the data can lead to
slightly different prior distributions for different subsets; see
[`estimate_unit_information_sd()`](https://fbartos.github.io/RoBMA/dev/reference/estimate_unit_information_sd.md)).
The specified `prior_unit_information_sd` is then used to scale the
default prior distributions as described in section (1).

Note that the manually specified `prior_unit_information_sd` takes
precedence over the estimated `UISD` from `ni` (section (2)) and the
known `UISD` from `measure` (section (1)). It cannot be combined with
`prior_informed_field`.

### (4) Specifying informed empirical prior distributions

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

### (5) Fully custom prior distributions

Prior distributions can be fully customized by directly specifying the
`prior_effect`, `prior_heterogeneity`, `prior_mods`, `prior_scale`, and
`prior_heterogeneity_allocation` arguments. These should be prior
distribution objects created via
[`BayesTools::prior()`](https://fbartos.github.io/BayesTools/reference/prior.html)
or related functions (e.g.,
[`prior_factor()`](https://fbartos.github.io/RoBMA/dev/reference/prior_factor.md)).
For
[`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
random-effect formulas, see
[`random_effect_prior_specification`](https://fbartos.github.io/RoBMA/dev/reference/random_effect_prior_specification.md).

### Rescaling prior distributions

The `rescale_priors` argument allows rescaling supported prior
distributions by a multiplicative factor. For example,
`rescale_priors = 2` doubles the standard deviations/scales of normal,
Cauchy, t, and inverse-gamma prior distributions, making them more
diffuse. Point and none priors are unchanged. For publication-bias prior
distributions, see
[`publication_bias_prior_specification`](https://fbartos.github.io/RoBMA/dev/reference/publication_bias_prior_specification.md).

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

[`publication_bias_prior_specification`](https://fbartos.github.io/RoBMA/dev/reference/publication_bias_prior_specification.md),
[`random_effect_prior_specification`](https://fbartos.github.io/RoBMA/dev/reference/random_effect_prior_specification.md),
[`prior`](https://fbartos.github.io/BayesTools/reference/prior.html),
[`RoBMA.options`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA_options.md),
[`brma`](https://fbartos.github.io/RoBMA/dev/reference/brma.md)
