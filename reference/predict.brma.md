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
  V_new = NULL,
  conditioning_depth = "marginal",
  ...
)
```

## Arguments

- object:

  a fitted brma object

- newdata:

  prediction data. Defaults to `NULL`, which predicts the observed
  design rows. `TRUE` returns a single aggregated prediction at the
  average fitted design (the scalar parameter for models without
  moderators or scale predictors, otherwise the arithmetic mean across
  the fitted rows of the location and of the observation-specific
  heterogeneity): the pooled location (`type = "terms"`), the pooled
  heterogeneity (`"terms.scale"`), a new true effect (`"estimate"`), or
  a new cluster effect added to the pooled location (`"cluster"`,
  multilevel models only). It is the quantity reported by
  [`pooled_effect`](https://fbartos.github.io/RoBMA/reference/pooled_effect.md)
  and
  [`pooled_heterogeneity`](https://fbartos.github.io/RoBMA/reference/pooled_heterogeneity.md)
  without the prediction interval columns, requires marginal depth, and
  is not available for `type = "response"`, which needs
  observation-level sampling variance. Marginal new-effect predictions
  may use an explicit data frame or named list with one row per target
  effect and every variable used by the requested fixed, scale, and
  random-effect designs. Missing ordinary grouping-only identifiers are
  synthesized as row-specific new levels; known-`R` grouping identifiers
  must be supplied. Outcome columns are optional for `type = "terms"`
  unless PET/PEESE bias terms are included via `bias_adjusted = FALSE`.
  Missing values are rejected; prediction never drops rows silently.
  Binomial GLMM response predictions require arm totals `n1i` and `n2i`;
  supply the totals alone, optionally with event counts `ai` and `ci`,
  or supply the four cell counts `ai`, `bi`, `ci`, and `di`. Supplied
  cells and totals must agree. Poisson GLMM response predictions require
  exposures `t1i` and `t2i`; event counts `x1i` and `x2i` are optional.
  For known-`V`
  [`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
  response prediction, explicit `newdata` also requires sampling
  covariance through `V_new`. Product-only selection models require no
  publication-group column. Unadjusted normal selected latent
  predictions (`type = "estimate"`) require `sei` or `vi` in explicit
  `newdata`. If an active selection branch uses `weight_rule = "best"`,
  its publication groups must also be resolvable in the prediction data.

- type:

  type of prediction to be performed. Options are:

  - `"terms"` (alias: `"marginal"`): fixed location parameters only,
    \\X\beta\\.

  - `"estimate"` (alias: `"effect"`): latent true effects. The
    random-effect distribution is selected by `conditioning_depth`.

  - `"response"` (alias: `"outcome"`): Predicted observed values (yi).
    Adds the outcome sampling distribution to the corresponding
    latent-effect target.

  - `"terms.scale"`: Scale parameter (tau), incorporating scale
    regression if present. For
    [`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
    random-formula models, returns a named list of component-specific
    `brma_samples` matrices, whether or not scale formulas are present.

  The released `"cluster"` and `"blup"` spellings are retained as
  compatibility shortcuts for conditional location samples: they fix
  `conditioning_depth = "cluster"` and `"estimate"`, respectively. Use
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) or
  [`blup()`](https://fbartos.github.io/RoBMA/reference/blup.md) when
  conditional means are the intended quantity.

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

  For new outcomes, selection predictions draw retained contexts from
  their original Gaussian law, then sample each full publication event
  conditional on that context. Selected latent predictions additionally
  reconstruct the integrated true effects conditional on the selected
  response. At estimate depth, response prediction is a new replication
  event conditional on the fitted latent true effects; LOO-PIT instead
  uses deletion within the original publication event. Retained
  estimate-level effects use their fitted posterior values at estimate
  depth and fresh Gaussian draws for a new estimate within a fitted
  cluster. The fitted selection specification and prediction
  conditioning depth are separate choices.

- quiet:

  logical; whether to suppress the informational message when
  conditional posterior samples are flattened to one chain.

- conditional:

  whether to return conditional posterior predictions for RoBMA
  product-space objects. For location predictions, samples are
  conditioned on the effect component; for `type = "terms.scale"`,
  samples are conditioned on the heterogeneity component. Conditional
  samples are flattened to one chain after subsetting posterior rows.

- V_new:

  optional sampling covariance for explicit `newdata` response or
  selected latent predictions from known-`V`
  [`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
  models. Accepts a vector of diagonal variances, a square covariance
  matrix, a list of block covariance matrices, or a
  [`known_v_factor`](https://fbartos.github.io/RoBMA/reference/known_v_factor.md)
  declaration. Its total dimension must match `nrow(newdata)`. It is
  required for explicit-newdata response predictions and unadjusted
  selected latent predictions. Cross-covariance with observed rows is
  not supported. Numerical symmetry is handled as for `V`: accepted
  off-diagonal pairs are averaged internally; the supplied object and
  diagonal are unchanged.

- conditioning_depth:

  random-effect conditioning depth for `type = "estimate"` and
  `"response"`. `"marginal"` (default) predicts new latent effects,
  `"cluster"` conditions on the fitted cluster effect and predicts a new
  estimate within that cluster, and `"estimate"` predicts the fitted
  latent effect. Cluster depth is available only for specialized
  multilevel models fitted with `cluster`. Non-marginal prediction is
  currently restricted to `newdata = NULL`, where fitted identities are
  unambiguous.

- ...:

  unused additional arguments. Supplied arguments trigger a warning.

## Value

A `brma_samples` object containing posterior samples. For
[`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
random-formula `type = "terms.scale"` predictions, returns a named list
of component-specific `brma_samples` objects; with `newdata = TRUE` this
component result is the `component = "all"` result of
[`pooled_heterogeneity`](https://fbartos.github.io/RoBMA/reference/pooled_heterogeneity.md)
(a named list for several components, one `brma_samples` object for a
single component), computed with the arithmetic mean of the
observation-specific tau of a scale formula; its components are titled
“Pooled Heterogeneity (`<block>`)” and “Total Pooled Heterogeneity”.
When printed, each `brma_samples` object displays a summary table via
[`BayesTools::ensemble_estimates_table`](https://fbartos.github.io/BayesTools/reference/BayesTools_ensemble_tables.html).
The underlying samples matrix can be accessed directly because the
object inherits from matrix. For a single `brma_samples` object, use
[`summary()`](https://rdrr.io/r/base/summary.html), or
`as.data.frame(summary(x))`, to obtain the summary table with component
and parameter identifiers; `as.data.frame(x)` returns the posterior
draws. A list of results (a `brma_samples_list`) coerces directly:
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) of the
list gives the summary table, and `as.data.frame(format = "list")`
retains separate tables. The samples can also be converted to posterior
draws formats using
[`as_draws()`](https://fbartos.github.io/RoBMA/reference/as_draws.brma.md)
and related functions.

## Details

Prediction has two independent axes. `type` selects the quantity: fixed
location (`"terms"`), latent true effect (`"estimate"`), or observed
outcome (`"response"`). `conditioning_depth` selects how much of the
fitted random-effect hierarchy is retained. In the specialized
multilevel model fitted with `cluster`, \\y\_{ij} = X\_{ij}\beta + u_j +
v\_{ij} + \epsilon\_{ij}\\, the targets are:

- `"marginal"`: condition on neither \\u_j\\ nor \\v\_{ij}\\;

- `"cluster"`: condition on fitted \\u_j\\, but draw a new \\v\_{ij}\\;

- `"estimate"`: draw the fitted latent effect from its posterior.

`type = "response"` adds \\\epsilon\\ to the selected latent-effect
target. For normal fitted effects this includes the conditional latent
uncertainty, not only the BLUP mean.
[`blup()`](https://fbartos.github.io/RoBMA/reference/blup.md) and
`fitted(..., conditioning_depth = "estimate")` remain the
conditional-mean interfaces.

For admitted weighted Gaussian selection models, reconstruction of
fitted integrated true effects uses the normalized conditional kernel
\\N(u;0,\tau^2)N(y;c+u,s^2)^w\\. Its conditioning variance is \\s^2/w\\;
the fitted global posterior and retained latent contributions remain the
supplied posterior. This post-fit augmentation does not assert that
likelihood powering commutes with integrating a latent effect.
Estimate-depth estimate and response targets record this policy and the
row-aligned likelihood weights in `RoBMA_target$latent_augmentation`.
Response replication adds the original sampling variance \\s^2\\.

Ordinary independent Gaussian selection predictions use exact
diagonal-cell chunk pricing after conservative fixed-geometry and output
reserves. Other model structures retain the dense covariance charge
under the same budget. The `RoBMA.known_v_covariance_max_bytes` option
controls covariance chunk memory. Selected predictions are reproducible
for a fixed seed and chunk budget, but changing this option can change
the returned seeded draw sequence because chunking changes the order of
random-number consumption.

`newdata` supplies design values and grouping identities; it does not
select an estimand. Consequently, marginal prediction has the same law
for `newdata = NULL` and an equivalent explicit design. Matching fitted
grouping labels retains the requested dependence among new marginal
draws but does not reuse fitted random effects. Cluster- and
estimate-depth predictions require fitted identities and therefore
currently reject explicit `newdata`.

For binomial and Poisson GLMM response predictions, estimate depth uses
the fitted posterior nuisance base rate \\\pi_i\\ or rate \\\phi_i\\.
Marginal and cluster depths describe a new estimate and therefore draw a
new nuisance rate independently from its model prior. Those predictions
contain a prior-predictive nuisance-rate component, which can dominate
their width when the nuisance prior is diffuse; fitted
\\\pi_i\\/\\\phi_i\\ values are never silently reused for a marginal or
cluster target.

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

For known-`V`
[`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
models, explicit `newdata` response and unadjusted selected latent
prediction require `V_new`. Random-formula models preserve the joint
random-effect structure selected by `conditioning_depth` and then add
`V` or `V_new` sampling noise. Cross-covariance with observed rows is
not supported. Explicit `newdata` estimate/response predictions are not
available for
[`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
objects with marginalized known-`R` blocks, because their row
multipliers are fitted-row metadata.
[`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
prediction results carry a `brma_mv_prediction_target` attribute
recording the formula, random-effect, and response-generation target
used.

For known-`V` models, the conditional random-effect predictions require
applying the inverse of the marginal covariance. A covariance that is
too ill-conditioned for that (a reciprocal condition number of about
`2e-10` or less in any posterior draw used) is unavailable: the call
stops with an error of class `RoBMA_known_v_precision_unavailable`,
whose parent class is `RoBMA_numerical_unavailable`, instead of
returning a rounded number. The condition carries the draw positions
(`rows`) and their `reciprocal_condition`.

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
  predict(fit, newdata = data.frame(ablat = 50, year = 1970),
          type = "estimate")
  predict(fit, type = "estimate")
}
} # }
```
