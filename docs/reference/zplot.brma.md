# Plot Zplot Diagnostics Directly

Convenience wrapper for creating and plotting zplot diagnostics from a
fitted `brma` object.

## Usage

``` r
# S3 method for class 'brma'
zplot(
  object,
  significance_level = stats::qnorm(0.975),
  summary_max_samples = 10000,
  ...
)

# S3 method for class 'zplot_brma'
zplot(object, ...)
```

## Arguments

- object:

  a normal-outcome `brma` object, or a `zplot_brma` object.

- significance_level:

  z-value threshold for significance. Defaults to `qnorm(0.975)`
  (two-sided alpha = 0.05).

- summary_max_samples:

  maximum number of posterior samples used for the EDR and missing-study
  summaries stored in the generated zplot object. This is separate from
  the plot-density `max_samples` argument accepted by
  [`plot.zplot_brma()`](https://fbartos.github.io/RoBMA/reference/plot.zplot_brma.md)
  through `...`. Defaults to 10000. Use `Inf` to use all posterior
  samples.

- ...:

  arguments passed to
  [`plot.zplot_brma()`](https://fbartos.github.io/RoBMA/reference/plot.zplot_brma.md).

## Value

`NULL` invisibly for base graphics, or a ggplot2 object.

## Details

When `object` already inherits from `zplot_brma`, `zplot()` dispatches
directly to
[`plot.zplot_brma()`](https://fbartos.github.io/RoBMA/reference/plot.zplot_brma.md)
without recomputing stored summaries.

## See also

[`as_zplot.brma()`](https://fbartos.github.io/RoBMA/reference/as_zplot.brma.md),
[`plot.zplot_brma()`](https://fbartos.github.io/RoBMA/reference/plot.zplot_brma.md),
[`summary.zplot_brma()`](https://fbartos.github.io/RoBMA/reference/summary.zplot_brma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")
  fit <- bPET(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")

  zplot(fit)
}
} # }
```
