# Estimate a Robust Bayesian Meta-Analysis

`RoBMA` is used to estimate a robust Bayesian meta-analysis. The
interface allows a complete customization of the ensemble with different
prior (or list of prior) distributions for each component.

## Usage

``` r
RoBMA(
  d = NULL,
  r = NULL,
  logOR = NULL,
  OR = NULL,
  z = NULL,
  y = NULL,
  se = NULL,
  v = NULL,
  n = NULL,
  lCI = NULL,
  uCI = NULL,
  t = NULL,
  study_names = NULL,
  study_ids = NULL,
  data = NULL,
  weight = NULL,
  transformation = if (is.null(y)) "fishers_z" else "none",
  prior_scale = if (is.null(y)) "cohens_d" else "none",
  effect_direction = "positive",
  model_type = NULL,
  rescale_priors = 1,
  priors_effect = set_default_priors("effect", rescale = rescale_priors),
  priors_heterogeneity = set_default_priors("heterogeneity", rescale = rescale_priors),
  priors_bias = set_default_priors("bias", rescale = rescale_priors),
  priors_effect_null = set_default_priors("effect", null = TRUE),
  priors_heterogeneity_null = set_default_priors("heterogeneity", null = TRUE),
  priors_bias_null = set_default_priors("bias", null = TRUE),
  priors_hierarchical = set_default_priors("hierarchical"),
  priors_hierarchical_null = set_default_priors("hierarchical", null = TRUE),
  algorithm = "bridge",
  chains = 3,
  sample = 5000,
  burnin = 2000,
  adapt = 500,
  thin = 1,
  parallel = FALSE,
  autofit = TRUE,
  autofit_control = set_autofit_control(),
  convergence_checks = set_convergence_checks(),
  save = "all",
  seed = NULL,
  silent = TRUE,
  ...
)
```

## Arguments

- d:

  a vector of effect sizes measured as Cohen's d / Hedges' g
  (standardized mean differences)

- r:

  a vector of effect sizes measured as correlations

- logOR:

  a vector of effect sizes measured as log odds ratios

- OR:

  a vector of effect sizes measured as odds ratios

- z:

  a vector of effect sizes measured as Fisher's z

- y:

  a vector of unspecified effect sizes (note that effect size
  transformations are unavailable with this type of input)

- se:

  a vector of standard errors of the effect sizes

- v:

  a vector of variances of the effect sizes

- n:

  a vector of overall sample sizes

- lCI:

  a vector of lower bounds of confidence intervals

- uCI:

  a vector of upper bounds of confidence intervals

- t:

  a vector of t/z-statistics

- study_names:

  an optional argument with the names of the studies

- study_ids:

  an optional argument specifying dependency between the studies (for
  using a multilevel model). Defaults to `NULL` for studies being
  independent.

- data:

  a data object created by the `combine_data` function. This is an
  alternative input entry to specifying the `d`, `r`, `y`, etc...
  directly. I.e., RoBMA function does not allow passing a data.frame and
  referencing to the columns.

- weight:

  specifies likelihood weights of the individual estimates. Notes that
  this is an untested experimental feature.

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

- effect_direction:

  the expected direction of the effect. Correctly specifying the
  expected direction of the effect is crucial for one-sided selection
  models, as they specify cut-offs using one-sided p-values. Defaults to
  `"positive"` (another option is `"negative"`).

- model_type:

  string specifying the RoBMA ensemble. Defaults to `NULL`. The other
  options are `"PSMA"`, `"PP"`, and `"2w"` which override settings
  passed to the `priors_effect`, `priors_heterogeneity`,
  `priors_effect`, `priors_effect_null`, `priors_heterogeneity_null`,
  `priors_bias_null`, and `priors_effect`. See details for more
  information about the different model types.

- rescale_priors:

  a re-scaling factor for the prior distributions. The re-scaling factor
  allows to adjust the width of all default priors simultaneously.
  Defaults to `1`.

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

- algorithm:

  a string specifying the algorithm used for the model averaging.
  Defaults to `"bridge"` which results in estimating individual models
  using JAGS and computing the marginal likelihood using bridge
  sampling. An alternative is `"ss"` which uses spike and slab like
  parameterization to approximate the Bayesian model averaging with a
  single model. Note that significantly more `sample`, `burnin`, and
  `adapt` iterations are needed for the `"ss"` algorithm.

- chains:

  a number of chains of the MCMC algorithm.

- sample:

  a number of sampling iterations of the MCMC algorithm. Defaults to
  `5000`.

- burnin:

  a number of burnin iterations of the MCMC algorithm. Defaults to
  `2000`.

- adapt:

  a number of adaptation iterations of the MCMC algorithm. Defaults to
  `500`.

- thin:

  a thinning of the chains of the MCMC algorithm. Defaults to `1`.

- parallel:

  whether the individual models should be fitted in parallel. Defaults
  to `FALSE`. The implementation is not completely stable and might
  cause a connection error.

- autofit:

  whether the model should be fitted until the convergence criteria
  (specified in `autofit_control`) are satisfied. Defaults to `TRUE`.

- autofit_control:

  allows to pass autofit control settings with the
  [`set_autofit_control()`](reference/RoBMA_control.md) function. See
  [`?set_autofit_control`](reference/RoBMA_control.md) for options and
  default settings.

- convergence_checks:

  automatic convergence checks to assess the fitted models, passed with
  [`set_convergence_checks()`](reference/RoBMA_control.md) function. See
  [`?set_convergence_checks`](reference/RoBMA_control.md) for options
  and default settings.

- save:

  whether all models posterior distributions should be kept after
  obtaining a model-averaged result. Defaults to `"all"` which does not
  remove anything. Set to `"min"` to significantly reduce the size of
  final object, however, some model diagnostics and further manipulation
  with the object will not be possible.

- seed:

  a seed to be set before model fitting, marginal likelihood
  computation, and posterior mixing for reproducibility of results.
  Defaults to `NULL` - no seed is set.

- silent:

  whether all print messages regarding the fitting process should be
  suppressed. Defaults to `TRUE`. Note that `parallel = TRUE` also
  suppresses all messages.

- ...:

  additional arguments.

## Value

`RoBMA` returns an object of class 'RoBMA'.

## Details

The default settings of the RoBMA 2.0 package corresponds to the
RoBMA-PSMA ensemble proposed by Bartoš et al. (2023) . The previous
versions of the package (i.e., RoBMA \< 2.0) used specifications
proposed by Maier et al. (2023) (this specification can be easily
obtained by setting `model_type = "2w"`. The RoBMA-PP specification from
Bartoš et al. (2023) can be obtained by setting `model_type = "PP"`. The
complete list of default prior distributions is described at
[`set_default_priors()`](reference/set_default_priors.md). Note that
inclusion of the PET and PEESE style publication bias adjustments models
might pick up on small-study effects. To remove true heterogeneity due
to study design, sub-populations, treatments etc. potentially causing
small-study effects, use meta-regression via the
[`RoBMA.reg()`](reference/RoBMA.reg.md) function, or remove the PET and
PEESE style models from the publication bias adjustment component of the
ensemble.

The
[[`vignette("CustomEnsembles", package = "RoBMA")`](articles/CustomEnsembles.md)](doc/CustomEnsembles.md)
and
[[`vignette("ReproducingBMA", package = "RoBMA")`](articles/ReproducingBMA.md)](doc/ReproducingBMA.md)
vignettes describe how to use `RoBMA()` to fit custom meta-analytic
ensembles (see [`prior()`](reference/prior.md),
[`prior_weightfunction()`](reference/prior_weightfunction.md),
[`prior_PET()`](reference/prior_PET.md), and
[`prior_PEESE()`](reference/prior_PEESE.md) for more information about
prior distributions).

The RoBMA function first generates models from a combination of the
provided priors for each of the model parameters. Then, the individual
models are fitted using
[autorun.jags](https://rdrr.io/pkg/runjags/man/autorun.jags.html)
function. A marginal likelihood is computed using
[bridge_sampler](https://rdrr.io/pkg/bridgesampling/man/bridge_sampler.html)
function. The individual models are then combined into an ensemble using
the posterior model probabilities using
[BayesTools](https://fbartos.github.io/BayesTools/reference/BayesTools.html)
package.

Generic [`summary.RoBMA()`](reference/summary.RoBMA.md),
[`print.RoBMA()`](reference/print.RoBMA.md), and
[`plot.RoBMA()`](reference/plot.RoBMA.md) functions are provided to
facilitate manipulation with the ensemble. A visual check of the
individual model diagnostics can be obtained using the
[`diagnostics()`](reference/diagnostics.md) function. The fitted model
can be further updated or modified by
[`update.RoBMA()`](reference/update.RoBMA.md) function.

## References

Bartoš F, Maier M, Wagenmakers E, Doucouliagos H, Stanley TD (2023).
“Robust Bayesian meta-analysis: Model-averaging across complementary
publication bias adjustment methods.” *Research Synthesis Methods*,
**14**(1), 99–116.
[doi:10.1002/jrsm.1594](https://doi.org/10.1002/jrsm.1594) .  
  
Maier M, Bartoš F, Wagenmakers E (2023). “Robust Bayesian Meta-Analysis:
Addressing publication bias with model-averaging.” *Psychological
Methods*, **28**(1), 107–122.
[doi:10.1037/met0000405](https://doi.org/10.1037/met0000405) .  
  
van Erp S, Verhagen J, Grasman RP, Wagenmakers E (2017). “Estimates of
between-study heterogeneity for 705 meta-analyses reported in
Psychological Bulletin from 1990–2013.” *Journal of Open Psychology
Data*, **5**(1), 1–5.
[doi:10.5334/jopd.33](https://doi.org/10.5334/jopd.33) .

## See also

[`summary.RoBMA()`](reference/summary.RoBMA.md),
[`update.RoBMA()`](reference/update.RoBMA.md),
[`check_setup()`](reference/check_setup.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# using the example data from Bem 2011 and fitting the default (RoBMA-PSMA) model
fit <- RoBMA(d = Bem2011$d, se = Bem2011$se, study_names = Bem2011$study)

# in order to speed up the process, we can turn the parallelization on
fit <- RoBMA(d = Bem2011$d, se = Bem2011$se, study_names = Bem2011$study, parallel = TRUE)

# we can get a quick overview of the model coefficients just by printing the model
fit

# a more detailed overview using the summary function (see '?summary.RoBMA' for all options)
summary(fit)

# the model-averaged effect size estimate can be visualized using the plot function
# (see ?plot.RoBMA for all options)
plot(fit, parameter = "mu")

# forest plot can be obtained with the forest function (see ?forest for all options)
forest(fit)

# plot of the individual model estimates can be obtained with the plot_models function
#  (see ?plot_models for all options)
plot_models(fit)

# diagnostics for the individual parameters in individual models can be obtained using diagnostics
# function (see 'diagnostics' for all options)
diagnostics(fit, parameter = "mu", type = "chains")

# the RoBMA-PP can be fitted with addition of the 'model_type' argument
fit_PP <- RoBMA(d = Bem2011$d, se = Bem2011$se, study_names = Bem2011$study, model_type = "PP")

# as well as the original version of RoBMA (with two weightfunctions)
fit_original <- RoBMA(d = Bem2011$d, se = Bem2011$se, study_names = Bem2011$study,
                      model_type = "2w")

# or different prior distribution for the effect size (e.g., a half-normal distribution)
# (see 'vignette("CustomEnsembles")' for a detailed guide on specifying a custom model ensemble)
fit <- RoBMA(d = Bem2011$d, se = Bem2011$se, study_names = Bem2011$study,
             priors_effect = prior("normal", parameters = list(0, 1),
                                   truncation = list(0, Inf)))
} # }
```
