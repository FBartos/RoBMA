# Transform brma Object to Zplot

Transforms an estimated brma model into a zplot object that can be
summarized and plotted to assess replicability.

## Usage

``` r
# S3 method for class 'brma'
as_zplot(
  object,
  significance_level = stats::qnorm(0.975),
  max_samples = 1000,
  ...
)
```

## Arguments

- object:

  a normal-outcome `brma` object.

- significance_level:

  z-value threshold for significance. Defaults to `qnorm(0.975)`
  (two-sided alpha = 0.05).

- max_samples:

  maximum number of posterior samples for estimation. Defaults to 1000.

- ...:

  additional arguments (currently unused).

## Value

The input object with added class `"zplot_brma"` and a new `zplot` list
component containing:

- estimates:

  a list with posterior samples for `EDR` and missing-study `weights`

- data:

  a list with `significance_level`, observed `z`-statistics,
  `N_significant`, and `N_observed`

## Details

Zplot analysis estimates the Expected Discovery Rate (EDR), the
posterior mean probability that an exact replication would be
statistically significant at the supplied threshold. This provides
insight into the replicability of findings in a literature.

The implementation uses extrapolation mode by default, which removes
selection-model weights when evaluating the model-implied z-value
density. PET/PEESE regression terms remain part of the fitted location
model. Zplot diagnostics are available only for normal outcome models.
GLMM objects are rejected because their raw likelihood is on a count
scale while zplot diagnostics require observed effect-size z-statistics
with standard errors.

The resulting object retains all original brma properties while adding
zplot results, enabling both standard meta-analytic summaries and zplot
diagnostics on the same object.

## See also

[`summary.zplot_brma()`](https://fbartos.github.io/RoBMA/reference/summary.zplot_brma.md),
[`plot.zplot_brma()`](https://fbartos.github.io/RoBMA/reference/plot.zplot_brma.md),
[`hist.zplot_brma()`](https://fbartos.github.io/RoBMA/reference/hist.zplot_brma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")
  fit <- bPET(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")

  zfit <- as_zplot(fit)
  summary(zfit)
  plot(zfit)
}
} # }
```
