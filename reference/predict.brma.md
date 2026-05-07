# Predict From brma Object

`predict.brma` predicts values

## Usage

``` r
# S3 method for class 'brma'
predict(
  object,
  newdata = NULL,
  type = "terms",
  as_measure = TRUE,
  output_measure = NULL,
  transform = NULL,
  probs = c(0.025, 0.975),
  bias_adjusted = FALSE,
  quiet = FALSE,
  conditional = FALSE,
  ...
)
```

## Arguments

- object:

  a fitted brma object

- newdata:

  specification for prediction data. Defaults to `NULL` which
  corresponds to prediction for the observed data. Alternatives are:

  - `TRUE` returns a single aggregated prediction (either the scalar
    parameter for models without moderators/scale, or the average across
    the model matrix for regression models). Not available for
    `type = "response"` since observation-level sampling variance is
    required.

  - A data.frame (for meta-regression) or a named list with effect size
    measure and variability metrics (for meta-analysis) for new studies.
    The input must contain the variables used by the requested
    prediction type. Outcome columns are optional for `type = "terms"`
    unless PET/PEESE bias terms are included via
    `bias_adjusted = FALSE`.

- type:

  type of prediction to be performed. Options are:

  - `"terms"` (alias: `"marginal"`): Fixed-effect parameters only (mu).
    Produces the mean effect size estimate at the given predictor
    levels, not accounting for random effects.

  - `"cluster"`: Fixed effects plus cluster-level random effects (mu +
    gamma). Only available for multilevel (3-level) models.

  - `"estimate"` (aliases: `"effect"`, `"blup"`): True effects (mu +
    gamma + theta). For existing normal data, returns posterior draws of
    conditional BLUP means, not simulated latent-effect draws. For new
    data, samples estimate-level random effects.

  - `"response"` (alias: `"outcome"`): Predicted observed values (yi).
    Incorporates both heterogeneity and sampling variability.

  - `"terms.scale"`: Scale parameter (tau), incorporating scale
    regression if present.

- as_measure:

  logical; whether to convert GLMM response predictions from simulated
  counts to continuity-corrected effect-size estimators (logOR for
  binomial, logIRR for Poisson). Defaults to `TRUE`. Only relevant for
  GLMM models with `type = "response"`. When `FALSE`, returns raw
  frequency data (counts). Use `type = "estimate"` for latent
  logOR/logIRR predictions.

- output_measure:

  effect-size measure for location/effect predictions. Defaults to the
  fitted measure. Supported conversions are among `"SMD"`, `"COR"`,
  `"ZCOR"`, and `"OR"`; `"RR"`, `"HR"`, `"IRR"`, `"RD"`, and `"GEN"` can
  only be returned on their fitted measure. Use `transform = "EXP"` for
  ratio-scale output from log-scale measures.

- transform:

  optional display transformation. Currently `"EXP"` exponentiates
  log-scale measures `"OR"`, `"RR"`, `"HR"`, and `"IRR"`.

- probs:

  quantiles of the posterior samples to be displayed when the returned
  `brma_samples` object is printed. Defaults to `c(.025, .975)`.

- bias_adjusted:

  whether predictions should adjust for publication bias. Defaults to
  `FALSE`. When `TRUE`:

  - PET/PEESE terms are NOT added to the mean parameter (mu), returning
    the bias-corrected effect estimate.

  - For `type = "response"` with selection models, samples from the
    ordinary normal predictive distribution instead of the
    selected-normal distribution, simulating what would be observed
    without publication bias.

  When `FALSE`:

  - PET/PEESE terms ARE added to mu, returning predictions that include
    the expected bias (i.e., what we expect to observe given publication
    bias).

  - For `type = "response"` with selection models, samples from the
    selected-normal distribution reflecting the selective publishing
    process.

- quiet:

  logical; whether to suppress informational messages about prediction
  scale and bias adjustment.

- conditional:

  whether to return conditional posterior predictions for RoBMA
  product-space objects. For location predictions, samples are
  conditioned on the effect component; for `type = "terms.scale"`,
  samples are conditioned on the heterogeneity component. Conditional
  samples are flattened to one chain after subsetting posterior rows.

- ...:

  additional arguments

## Value

A `brma_samples` object containing posterior samples. When printed,
displays a summary table via
[`BayesTools::ensemble_estimates_table`](https://fbartos.github.io/BayesTools/reference/BayesTools_ensemble_tables.html).
The underlying samples matrix can be accessed directly (the object
inherits from matrix) or via
[`summary()`](https://rdrr.io/r/base/summary.html) to obtain the summary
table. The samples can also be converted to posterior draws formats
using
[`as_draws()`](https://fbartos.github.io/RoBMA/reference/as_draws.brma.md)
and related functions.

## Details

**Type hierarchy:**

- `"terms"`: mu (fixed effects only)

- `"cluster"`: mu + gamma (adds cluster-level random effect)

- `"estimate"`: mu + gamma + theta (adds estimate-level random effect)

- `"response"`: mu + gamma + theta + epsilon (adds sampling error)

For existing normal observations, `type = "estimate"` reports the
conditional BLUP mean \\E(\theta_i \mid y_i, \mu_i, \tau_i)\\ for each
posterior draw. It is therefore an empirical-Bayes summary, not a draw
from the full latent-effect posterior \\\theta_i \mid y_i\\. For new
data, estimate-level random effects are sampled from their model
distribution.

For RoBMA product-space objects, conditional posterior predictions
subset posterior rows according to model indicators. This removes the
original chain structure, so returned `brma_samples` objects are
intentionally stored as one flattened chain.

Note that in contrast to
[predict](https://wviechtb.github.io/metafor/reference/predict.rma.html),
the `type = "response"` produces predictions for the new effect size
estimates. To obtain results corresponding to metafor's predict
function, use `type = "terms"` for the mean effect size and
`type = "estimate"` for true effects (prediction interval).

**Likelihood weights:** If the model was fitted with `weights`, the
weights affect the posterior fit, log-likelihood, LOO/WAIC, and
existing-data conditional diagnostics such as BLUP shrinkage and
leverage. They do not change the observation-level sampling error used
by `type = "response"` for normal models: response predictions simulate
raw future effect-size estimates with the supplied `sei`, not
`sei / sqrt(weight)`.

## See also

[`pooled_effect()`](https://fbartos.github.io/RoBMA/reference/pooled_effect.md),
[`pooled_heterogeneity()`](https://fbartos.github.io/RoBMA/reference/pooled_heterogeneity.md),
[`blup()`](https://fbartos.github.io/RoBMA/reference/blup.md)

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

  predict(fit, type = "terms")
  predict(fit, newdata = TRUE, type = "terms")
  predict(fit, type = "estimate")
}
} # }
```
