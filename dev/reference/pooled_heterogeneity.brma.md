# Pooled Heterogeneity for brma Objects

Computes heterogeneity (tau) of a fitted brma object at the average
fitted scale design (the mean of the observation-specific tau) and, for
random-formula models, at the average expanded random-effect design.

## Usage

``` r
# S3 method for class 'brma'
pooled_heterogeneity(
  object,
  probs = c(0.025, 0.975),
  conditional = FALSE,
  component = "all",
  ...
)
```

## Arguments

- object:

  a fitted brma object

- probs:

  quantiles of the posterior distribution to be displayed. Defaults to
  `c(.025, .975)` for 95% credible intervals.

- conditional:

  whether to return the pooled heterogeneity conditional on the
  heterogeneity component for RoBMA product-space objects. Defaults to
  `FALSE`.

- component:

  heterogeneity component to return for
  [`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
  models. Defaults to `"all"`. Use `"total"` for the variance-additive
  total heterogeneity.

- ...:

  additional arguments passed to
  [`predict.brma`](https://fbartos.github.io/RoBMA/dev/reference/predict.brma.md)
  (`as_measure`, `V_new`, and `conditioning_depth`); wrapper arguments
  such as `newdata`, `type`, and `quiet` are fixed (supplying one is R's
  error for an argument matched twice, as for
  [`predict()`](https://rdrr.io/r/stats/predict.html)), and other
  unknown arguments are unused and trigger a warning.

## Value

A `brma_samples` object containing posterior samples. When printed,
displays a summary table. For decomposed
[`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
models, a named list of `brma_samples` objects is returned. Use
[`summary()`](https://rdrr.io/r/base/summary.html) on an individual
object ([`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
of an individual object returns its draws) or
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) on the
complete result to obtain a summary data frame with component and
parameter identifiers. Use `as.data.frame(format = "list")` to retain
separate component tables. The samples can be converted to posterior
draws formats using
[`as_draws()`](https://fbartos.github.io/RoBMA/dev/reference/as_draws.brma.md).

## Details

For location-scale models, the function returns the arithmetic mean over
the fitted rows of the observation-specific tau, as in 4.0.0; it equals
`rowMeans(predict(object, type = "terms.scale"))`. It is not the root
mean square heterogeneity across the observed rows. Use
[`summary_heterogeneity()`](https://fbartos.github.io/RoBMA/dev/reference/summary_heterogeneity.md)
for the latter.

For models without scale regression, this returns the single tau
parameter.

For multilevel (3-level) models, the returned tau is the total
heterogeneity: `tau = sqrt(tau_within^2 + tau_between^2)`.

For
[`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
models, one ordinary selected component is named `tau`; a genuine
variance-additive total is named `tau_total`; and a mean-variance
allocation anchor is named `tau_common`. These quantities are available
when a random formula is declared. For random-formula models,
`component = "all"` returns one `brma_samples` object when there is a
single heterogeneity component and a named list when there are multiple
components. `component = "total"` computes the variance-additive total
at the average random-effect design. Random slopes are evaluated as
\\\bar q^T G \bar q\\, including within-block correlations. Default and
named allocation components retain their declared scale anchors;
explicit totals and pooled-effect prediction intervals use the full
covariance of the average design, including allocation-owned blocks.
With a known group covariance \\R\\, the reported `sd` is the fitted
multiplier of that kernel. It need not equal every row's marginal
standard deviation when \\\mathrm{diag}(R)\\ is not one. Row-specific
marginal standard deviations are available from
`predict(type = "terms.scale")`. A scale formula of a random-formula
model enters the average random-effect design as the arithmetic mean of
the observation-specific tau, as for the univariate models.
`predict(object, newdata = TRUE, type = "terms.scale")` returns the same
component result: a named list for several components and one
`brma_samples` object for a single component.

For
[`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
models and models with factor-specific or random-slope terms, the
average expanded random-effect design can be synthetic and need not
correspond to a single realizable outcome or study. The result is the
heterogeneity of that linear composite, not typical heterogeneity across
the observed designs. Interpret it only when the composite has a
meaningful scientific interpretation; use
[`summary_heterogeneity()`](https://fbartos.github.io/RoBMA/dev/reference/summary_heterogeneity.md)
for observed-design heterogeneity or `predict(..., newdata = ...)` for a
specified realizable design.

## See also

[`predict.brma()`](https://fbartos.github.io/RoBMA/dev/reference/predict.brma.md),
[`pooled_effect()`](https://fbartos.github.io/RoBMA/dev/reference/pooled_effect.md),
[`blup()`](https://fbartos.github.io/RoBMA/dev/reference/blup.md)

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

  pooled_heterogeneity(fit)
}
} # }
```
