# Random Effects for brma Objects

Extracts random effect deviations from a fitted brma object. These are
posterior-sample offsets from the fixed-effect predictions, analogous to
random-effect deviations returned by
[`metafor::ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html).

## Usage

``` r
# S3 method for class 'brma'
ranef(
  object,
  bias_adjusted = FALSE,
  probs = c(0.025, 0.975),
  component = "all",
  simplify = TRUE,
  expand = TRUE,
  ...
)
```

## Arguments

- object:

  a fitted brma object

- bias_adjusted:

  whether to adjust for publication bias. Defaults to `FALSE`. See
  [`blup.brma`](https://fbartos.github.io/RoBMA/dev/reference/blup.brma.md)
  for details.

- probs:

  quantiles of the posterior distribution to be displayed. Defaults to
  `c(.025, .975)` for 95% credible intervals.

- component:

  random-effect component to return. Defaults to `"all"`. For decomposed
  models, use a component name such as `"cluster"`, `"estimate"`, or a
  [`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
  random-effect block name. Use `"total"` for the summed random-effect
  deviation.

- simplify:

  whether `component = "all"` should return a single `brma_samples`
  object instead of a one-element list. Defaults to `TRUE`, matching
  standard 2-level
  [`brma()`](https://fbartos.github.io/RoBMA/dev/reference/brma.md)
  behavior.

- expand:

  whether to repeat random-effect contributions for every fitted
  observation. Defaults to `TRUE` (one column per estimate; the cluster
  block of a multilevel model labels its columns
  `u_cluster[<cluster>|<estimate>]`), except for
  [`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
  random-formula models, where it defaults to `FALSE`. `FALSE` returns
  one column per unique grouping-level contribution in the fitted
  grouping order used by
  [`metafor::ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html)
  (`u_cluster[<cluster>]` for a multilevel model). Nested groups follow
  their original outer-to-inner order. Indicator-coded random
  coefficients retain one column per observed grouping-level and
  coefficient combination. Set to `TRUE` for observation-aligned output.

- ...:

  additional arguments forwarded to
  [`predict.brma`](https://fbartos.github.io/RoBMA/dev/reference/predict.brma.md)
  for supported options such as `conditional`. `newdata`, `type`,
  `quiet`, `output_measure`, and `transform` are controlled by this
  method.

## Value

A `brma_samples` object for a single selected/simplified component, or a
named list of `brma_samples` objects for decomposed `component = "all"`
output. Use [`summary()`](https://rdrr.io/r/base/summary.html) on a
single component
([`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) of it
returns its draws), or
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) on the
list, to obtain the displayed summary table in long form, or
`as.data.frame(format = "list")` to retain separate component tables.

## Details

Random effects are computed as the difference between the true effects
(BLUPs) and the fixed-effect predictions: \$\$u_i = \hat{\theta}\_i -
\hat{\mu}\_i\$\$

For standard (2-level) models, returns a single `brma_samples` object
with the estimate-level random effects.

For multilevel (3-level) models, returns a list with two `brma_samples`
matrices:

- `cluster`:

  Cluster-level random effects (\\\gamma_j \cdot \tau\_{between}\\),
  representing between-cluster deviations from the fixed effects.

- `estimate`:

  Estimate-level random effects (\\\theta_i - \mu_i - \gamma_j \cdot
  \tau\_{between}\\), representing within-cluster deviations from the
  cluster means.

For
[`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
random-formula models, decomposes the Gaussian conditional means by
random-effect block. The result is invariant to whether each block was
sampled or marginalized during fitting. If there is only one block and
`simplify = TRUE`, returns a single `brma_samples` object. Multiple
blocks are returned in one flat list under their canonical fitted block
names; no additional location layer or component-name prefix is added.
Unique-level output is available when a block's contribution is constant
within each grouping level or observed grouping-level and indicator-
coefficient combination. Other random slopes generally require
`expand = TRUE`. When the random-effect standard deviation of a block
varies across the rows of one returned cell (a grouping level, or a
grouping level and indicator coefficient) through a `scale` formula
(also for the cluster component of a clustered model with a `scale`
formula, whether normal, GLMM, or with a weight function), the
unique-level output is refused with an error of class
`RoBMA_ranef_expand_required` (also of class
`RoBMA_unsupported_request`); use `expand = TRUE`. Only the components
that are returned are judged: `component = "estimate"` of a clustered
model and a block name whose own scale does not vary are returned,
whereas `"all"` and `"total"` include every block. Summed
`component = "total"` output across multiple blocks also requires
`expand = TRUE` because different blocks need not share grouping levels.

## See also

[`blup.brma()`](https://fbartos.github.io/RoBMA/dev/reference/blup.brma.md),
[`predict.brma()`](https://fbartos.github.io/RoBMA/dev/reference/predict.brma.md),
[`pooled_effect()`](https://fbartos.github.io/RoBMA/dev/reference/pooled_effect.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")
  fit <- brma(
    yi      = yi,
    vi      = vi,
    data    = dat.lehmann2018,
    measure = "SMD",
    seed    = 1,
    silent  = TRUE
  )

  ranef(fit)
}
} # }
```
