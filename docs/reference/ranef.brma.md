# Random Effects for brma Objects

Extracts random effect deviations from a fitted brma object. These are
posterior-sample offsets from the fixed-effect predictions, analogous to
random-effect deviations returned by
[`metafor::ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html).

## Usage

``` r
# S3 method for class 'brma'
ranef(object, bias_adjusted = FALSE, probs = c(0.025, 0.975), ...)
```

## Arguments

- object:

  a fitted brma object

- bias_adjusted:

  whether to adjust for publication bias. Defaults to `FALSE`. See
  [`blup.brma`](https://fbartos.github.io/RoBMA/reference/blup.brma.md)
  for details.

- probs:

  quantiles of the posterior distribution to be displayed. Defaults to
  `c(.025, .975)` for 95% credible intervals.

- ...:

  additional arguments forwarded to
  [`predict.brma`](https://fbartos.github.io/RoBMA/reference/predict.brma.md)
  for supported options such as `conditional`. `newdata`, `type`,
  `quiet`, `output_measure`, and `transform` are controlled by this
  method.

## Value

For 2-level models, a `brma_samples` object. For 3-level models, a named
list of `brma_samples` objects (one per variance component).

## Details

Random effects are computed as the difference between the true effects
(BLUPs) and the fixed-effect predictions: \$\$u_i = \hat{\theta}\_i -
\hat{\mu}\_i\$\$

For standard (2-level) models, returns a single `brma_samples` object
with the estimate-level random effects.

For multilevel (3-level) models, returns a list with two
observation-aligned `brma_samples` matrices, one column per estimate
row:

- `cluster`:

  Cluster-level random effects (\\\gamma_j \cdot \tau\_{between}\\),
  representing between-cluster deviations from the fixed effects.

- `estimate`:

  Estimate-level random effects (\\\theta_i - \mu_i - \gamma_j \cdot
  \tau\_{between}\\), representing within-cluster deviations from the
  cluster means.

## See also

[`blup.brma()`](https://fbartos.github.io/RoBMA/reference/blup.brma.md),
[`predict.brma()`](https://fbartos.github.io/RoBMA/reference/predict.brma.md),
[`pooled_effect()`](https://fbartos.github.io/RoBMA/reference/pooled_effect.md)

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
