# Prints summary of `"RoBMA.reg"` ensemble implied by the specified priors and formula

`check_setup` prints summary of `"RoBMA.reg"` ensemble implied by the
specified prior distributions. It is useful for checking the ensemble
configuration prior to fitting all of the models.

`check_setup` prints summary of `"RoBMA.reg"` ensemble implied by the
specified prior distributions. It is useful for checking the ensemble
configuration prior to fitting all of the models.

## Usage

``` r
check_setup.reg(
  formula,
  data,
  test_predictors = TRUE,
  study_names = NULL,
  study_ids = NULL,
  transformation = if (any(colnames(data) != "y")) "fishers_z" else "none",
  prior_scale = if (any(colnames(data) != "y")) "cohens_d" else "none",
  standardize_predictors = TRUE,
  effect_direction = "positive",
  priors = NULL,
  model_type = NULL,
  priors_effect = prior(distribution = "normal", parameters = list(mean = 0, sd = 1)),
  priors_heterogeneity = prior(distribution = "invgamma", parameters = list(shape = 1,
    scale = 0.15)),
  priors_bias = list(prior_weightfunction(distribution = "two.sided", parameters =
    list(alpha = c(1, 1), steps = c(0.05)), prior_weights = 1/12),
    prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1,
    1), steps = c(0.05, 0.1)), prior_weights = 1/12), prior_weightfunction(distribution =
    "one.sided", parameters = list(alpha = c(1, 1), steps = c(0.05)), prior_weights =
    1/12), prior_weightfunction(distribution = "one.sided", parameters = list(alpha =
    c(1, 1, 1), steps = c(0.025, 0.05)), prior_weights = 1/12), 
    
    prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1,
    1), steps = c(0.05, 0.5)), prior_weights = 1/12), prior_weightfunction(distribution =
    "one.sided", parameters = list(alpha = c(1, 1, 1, 1), steps = c(0.025, 0.05, 0.5)),
    prior_weights = 1/12), prior_PET(distribution = "Cauchy", parameters = list(0, 1),
    truncation = list(0, Inf), prior_weights = 1/4), prior_PEESE(distribution = "Cauchy",
    parameters = list(0, 5), truncation = list(0, Inf), prior_weights = 1/4)),
  priors_effect_null = prior(distribution = "point", parameters = list(location = 0)),
  priors_heterogeneity_null = prior(distribution = "point", parameters = list(location =
    0)),
  priors_bias_null = prior_none(),
  priors_hierarchical = prior("beta", parameters = list(alpha = 1, beta = 1)),
  priors_hierarchical_null = NULL,
  prior_covariates = prior("normal", parameters = list(mean = 0, sd = 0.25)),
  prior_covariates_null = prior("spike", parameters = list(location = 0)),
  prior_factors = prior_factor("mnormal", parameters = list(mean = 0, sd = 0.25),
    contrast = "meandif"),
  prior_factors_null = prior("spike", parameters = list(location = 0)),
  models = FALSE,
  silent = FALSE,
  ...
)

check_setup.RoBMA.reg(
  formula,
  data,
  test_predictors = TRUE,
  study_names = NULL,
  study_ids = NULL,
  transformation = if (any(colnames(data) != "y")) "fishers_z" else "none",
  prior_scale = if (any(colnames(data) != "y")) "cohens_d" else "none",
  standardize_predictors = TRUE,
  effect_direction = "positive",
  priors = NULL,
  model_type = NULL,
  priors_effect = prior(distribution = "normal", parameters = list(mean = 0, sd = 1)),
  priors_heterogeneity = prior(distribution = "invgamma", parameters = list(shape = 1,
    scale = 0.15)),
  priors_bias = list(prior_weightfunction(distribution = "two.sided", parameters =
    list(alpha = c(1, 1), steps = c(0.05)), prior_weights = 1/12),
    prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1,
    1), steps = c(0.05, 0.1)), prior_weights = 1/12), prior_weightfunction(distribution =
    "one.sided", parameters = list(alpha = c(1, 1), steps = c(0.05)), prior_weights =
    1/12), prior_weightfunction(distribution = "one.sided", parameters = list(alpha =
    c(1, 1, 1), steps = c(0.025, 0.05)), prior_weights = 1/12), 
    
    prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1,
    1), steps = c(0.05, 0.5)), prior_weights = 1/12), prior_weightfunction(distribution =
    "one.sided", parameters = list(alpha = c(1, 1, 1, 1), steps = c(0.025, 0.05, 0.5)),
    prior_weights = 1/12), prior_PET(distribution = "Cauchy", parameters = list(0, 1),
    truncation = list(0, Inf), prior_weights = 1/4), prior_PEESE(distribution = "Cauchy",
    parameters = list(0, 5), truncation = list(0, Inf), prior_weights = 1/4)),
  priors_effect_null = prior(distribution = "point", parameters = list(location = 0)),
  priors_heterogeneity_null = prior(distribution = "point", parameters = list(location =
    0)),
  priors_bias_null = prior_none(),
  priors_hierarchical = prior("beta", parameters = list(alpha = 1, beta = 1)),
  priors_hierarchical_null = NULL,
  prior_covariates = prior("normal", parameters = list(mean = 0, sd = 0.25)),
  prior_covariates_null = prior("spike", parameters = list(location = 0)),
  prior_factors = prior_factor("mnormal", parameters = list(mean = 0, sd = 0.25),
    contrast = "meandif"),
  prior_factors_null = prior("spike", parameters = list(location = 0)),
  models = FALSE,
  silent = FALSE,
  ...
)

check_setup.reg(
  formula,
  data,
  test_predictors = TRUE,
  study_names = NULL,
  study_ids = NULL,
  transformation = if (any(colnames(data) != "y")) "fishers_z" else "none",
  prior_scale = if (any(colnames(data) != "y")) "cohens_d" else "none",
  standardize_predictors = TRUE,
  effect_direction = "positive",
  priors = NULL,
  model_type = NULL,
  priors_effect = prior(distribution = "normal", parameters = list(mean = 0, sd = 1)),
  priors_heterogeneity = prior(distribution = "invgamma", parameters = list(shape = 1,
    scale = 0.15)),
  priors_bias = list(prior_weightfunction(distribution = "two.sided", parameters =
    list(alpha = c(1, 1), steps = c(0.05)), prior_weights = 1/12),
    prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1,
    1), steps = c(0.05, 0.1)), prior_weights = 1/12), prior_weightfunction(distribution =
    "one.sided", parameters = list(alpha = c(1, 1), steps = c(0.05)), prior_weights =
    1/12), prior_weightfunction(distribution = "one.sided", parameters = list(alpha =
    c(1, 1, 1), steps = c(0.025, 0.05)), prior_weights = 1/12), 
    
    prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1,
    1), steps = c(0.05, 0.5)), prior_weights = 1/12), prior_weightfunction(distribution =
    "one.sided", parameters = list(alpha = c(1, 1, 1, 1), steps = c(0.025, 0.05, 0.5)),
    prior_weights = 1/12), prior_PET(distribution = "Cauchy", parameters = list(0, 1),
    truncation = list(0, Inf), prior_weights = 1/4), prior_PEESE(distribution = "Cauchy",
    parameters = list(0, 5), truncation = list(0, Inf), prior_weights = 1/4)),
  priors_effect_null = prior(distribution = "point", parameters = list(location = 0)),
  priors_heterogeneity_null = prior(distribution = "point", parameters = list(location =
    0)),
  priors_bias_null = prior_none(),
  priors_hierarchical = prior("beta", parameters = list(alpha = 1, beta = 1)),
  priors_hierarchical_null = NULL,
  prior_covariates = prior("normal", parameters = list(mean = 0, sd = 0.25)),
  prior_covariates_null = prior("spike", parameters = list(location = 0)),
  prior_factors = prior_factor("mnormal", parameters = list(mean = 0, sd = 0.25),
    contrast = "meandif"),
  prior_factors_null = prior("spike", parameters = list(location = 0)),
  models = FALSE,
  silent = FALSE,
  ...
)
```

## Arguments

- formula:

  a formula for the meta-regression model

- data:

  a data.frame containing the data for the meta-regression. Note that
  the column names have to correspond to the effect sizes (`d`, `logOR`,
  `OR`, `r`, `z`), a measure of sampling variability (`se`, `v`, `n`,
  `lCI`, `uCI`, `t`), and the predictors. See
  [`combine_data()`](reference/combine_data.md) for a complete list of
  reserved names and additional information about specifying input data.

- test_predictors:

  vector of predictor names to test for the presence of moderation
  (i.e., assigned both the null and alternative prior distributions).
  Defaults to `TRUE`, all predictors are tested using the default prior
  distributions (i.e., `prior_covariates`, `prior_covariates_null`,
  `prior_factors`, and `prior_factors_null`). To only estimate and
  adjust for the effect of predictors use `FALSE`. If `priors` is
  specified, any settings in `test_predictors` is overridden.

- study_names:

  an optional argument with the names of the studies

- study_ids:

  an optional argument specifying dependency between the studies (for
  using a multilevel model). Defaults to `NULL` for studies being
  independent.

- transformation:

  transformation to be applied to the supplied effect sizes before
  fitting the individual models. Defaults to `"fishers_z"`. We highly
  recommend using `"fishers_z"` transformation since it is the only
  variance stabilizing measure and does not bias PET and PEESE style
  models. The other options are `"cohens_d"`, correlation coefficient
  `"r"` and `"logOR"`. Supplying `"none"` will treat the effect sizes as
  unstandardized and refrain from any transformations.

- prior_scale:

  an effect size scale used to define priors. Defaults to `"cohens_d"`.
  Other options are `"fishers_z"`, correlation coefficient `"r"`, and
  `"logOR"`. The prior scale does not need to match the effect sizes
  measure - the samples from prior distributions are internally
  transformed to match the `transformation` of the data. The
  `prior_scale` corresponds to the effect size scale of default output,
  but can be changed within the summary function.

- standardize_predictors:

  whether continuous predictors should be standardized prior to
  estimating the model. Defaults to `TRUE`. Continuous predictor
  standardization is important for applying the default prior
  distributions for continuous predictors. Note that the resulting
  output corresponds to standardized meta-regression coefficients.

- effect_direction:

  the expected direction of the effect. Correctly specifying the
  expected direction of the effect is crucial for one-sided selection
  models, as they specify cut-offs using one-sided p-values. Defaults to
  `"positive"` (another option is `"negative"`).

- priors:

  named list of prior distributions for each predictor (with names
  corresponding to the predictors). It allows users to specify both the
  null and alternative hypothesis prior distributions for each predictor
  by assigning the corresponding element of the named list with another
  named list (with `"null"` and `"alt"`). If only one prior is specified
  for a given parameter, it is assumed to correspond to the alternative
  hypotheses and the default null hypothesis is specified (i.e.,
  `prior_covariates_null` or `prior_factors_null`). If a named list with
  only one named prior distribution is provided (either `"null"` or
  `"alt"`), only this prior distribution is used and no default
  distribution is filled in. Parameters without specified prior
  distributions are assumed to be only adjusted for using the default
  alternative hypothesis prior distributions (i.e., `prior_covariates`
  or `prior_factors`). If `priors` is specified, `test_predictors` is
  ignored.

- model_type:

  string specifying the RoBMA ensemble. Defaults to `NULL`. The other
  options are `"PSMA"`, `"PP"`, and `"2w"` which override settings
  passed to the `priors_effect`, `priors_heterogeneity`,
  `priors_effect`, `priors_effect_null`, `priors_heterogeneity_null`,
  `priors_bias_null`, and `priors_effect`. See details for more
  information about the different model types.

- priors_effect:

  list of prior distributions for the effect size (`mu`) parameter that
  will be treated as belonging to the alternative hypothesis. Defaults
  to a standard normal distribution
  `prior(distribution = "normal", parameters = list(mean = 0, sd = 1))`.

- priors_heterogeneity:

  list of prior distributions for the heterogeneity `tau` parameter that
  will be treated as belonging to the alternative hypothesis. Defaults
  to
  `prior(distribution = "invgamma", parameters = list(shape = 1, scale = .15))`
  that is based on heterogeneities estimates from psychology (van Erp et
  al. 2017) .

- priors_bias:

  list of prior distributions for the publication bias adjustment
  component that will be treated as belonging to the alternative
  hypothesis. Defaults to
  `list( prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1), steps = c(0.05)), prior_weights = 1/12), prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1, 1), steps = c(0.05, 0.10)), prior_weights = 1/12), prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1), steps = c(0.05)), prior_weights = 1/12), prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1), steps = c(0.025, 0.05)), prior_weights = 1/12), prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1), steps = c(0.05, 0.5)), prior_weights = 1/12), prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1, 1), steps = c(0.025, 0.05, 0.5)), prior_weights = 1/12), prior_PET(distribution = "Cauchy", parameters = list(0,1), truncation = list(0, Inf), prior_weights = 1/4), prior_PEESE(distribution = "Cauchy", parameters = list(0,5), truncation = list(0, Inf), prior_weights = 1/4) )`,
  corresponding to the RoBMA-PSMA model introduce by Bartoš et
  al. (2023) .

- priors_effect_null:

  list of prior distributions for the effect size (`mu`) parameter that
  will be treated as belonging to the null hypothesis. Defaults to a
  point null hypotheses at zero,
  `prior(distribution = "point", parameters = list(location = 0))`.

- priors_heterogeneity_null:

  list of prior distributions for the heterogeneity `tau` parameter that
  will be treated as belonging to the null hypothesis. Defaults to a
  point null hypotheses at zero (a fixed effect meta-analytic models),
  `prior(distribution = "point", parameters = list(location = 0))`.

- priors_bias_null:

  list of prior weight functions for the `omega` parameter that will be
  treated as belonging to the null hypothesis. Defaults no publication
  bias adjustment, [`prior_none()`](reference/prior_none.md).

- priors_hierarchical:

  list of prior distributions for the correlation of random effects
  (`rho`) parameter that will be treated as belonging to the alternative
  hypothesis. This setting allows users to fit a hierarchical
  (three-level) meta-analysis when `study_ids` are supplied. Note that
  this is an experimental feature and see News for more details.
  Defaults to a beta distribution
  `prior(distribution = "beta", parameters = list(alpha = 1, beta = 1))`.

- priors_hierarchical_null:

  list of prior distributions for the correlation of random effects
  (`rho`) parameter that will be treated as belonging to the null
  hypothesis. Defaults to `NULL`.

- prior_covariates:

  a prior distributions for the regression parameter of continuous
  covariates on the effect size under the alternative hypothesis (unless
  set explicitly in `priors`). Defaults to a relatively wide normal
  distribution
  `prior(distribution = "normal", parameters = list(mean = 0, sd = 0.25))`.

- prior_covariates_null:

  a prior distributions for the regression parameter of continuous
  covariates on the effect size under the null hypothesis (unless set
  explicitly in `priors`). Defaults to a no effect
  `prior("spike", parameters = list(location = 0))`.

- prior_factors:

  a prior distributions for the regression parameter of categorical
  covariates on the effect size under the alternative hypothesis (unless
  set explicitly in `priors`). Defaults to a relatively wide
  multivariate normal distribution specifying differences from the mean
  contrasts
  `prior_factor("mnormal", parameters = list(mean = 0, sd = 0.25), contrast = "meandif")`.

- prior_factors_null:

  a prior distributions for the regression parameter of categorical
  covariates on the effect size under the null hypothesis (unless set
  explicitly in `priors`). Defaults to a no effect
  `prior("spike", parameters = list(location = 0))`.

- models:

  should the models' details be printed.

- silent:

  do not print the results.

- ...:

  additional arguments.

## Value

`check_setup.reg` invisibly returns list of summary tables.

`check_setup.reg` invisibly returns list of summary tables.

## See also

[`check_setup()`](reference/check_setup.md)
[`RoBMA.reg()`](reference/RoBMA.reg.md)

[`check_setup()`](reference/check_setup.md)
[`RoBMA.reg()`](reference/RoBMA.reg.md)

## Examples

``` r
# check regression setup with example data
data(Andrews2021)
check_setup.reg(~ measure + age, data = Andrews2021)
#> Warning: You are about to estimate 144 models based on the model formula and prior specification.
#> Robust Bayesian meta-regression (set-up)
#> Components summary:
#>                Models Prior prob.
#> Effect         72/144       0.500
#> Heterogeneity  72/144       0.500
#> Bias          128/144       0.500
#> 
#> Meta-regression components summary:
#>         Models Prior prob.
#> measure 72/144       0.500
#> age     72/144       0.500

# check setup with custom priors
check_setup.reg(~ measure + age, data = Andrews2021,
                priors_effect = prior("normal", list(mean = 0, sd = 0.5),
                prior_weight = 2))
#> Warning: You are about to estimate 144 models based on the model formula and prior specification.
#> Robust Bayesian meta-regression (set-up)
#> Components summary:
#>                Models Prior prob.
#> Effect         72/144       0.667
#> Heterogeneity  72/144       0.500
#> Bias          128/144       0.500
#> 
#> Meta-regression components summary:
#>         Models Prior prob.
#> measure 72/144       0.500
#> age     72/144       0.500
```
