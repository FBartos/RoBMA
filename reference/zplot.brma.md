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
  conditioning_depth = "marginal",
  integration_control = set_selection_likelihood_control(),
  ...
)

# S3 method for class 'zplot_brma'
zplot(object, ...)
```

## Arguments

- object:

  a normal-outcome `brma` object, or a `zplot_brma` object.

- significance_level:

  positive z-value threshold for significance. Defaults to
  `qnorm(0.975)` (two-sided alpha = 0.05). The threshold must be finite
  and its two-sided alpha,
  `2 * pnorm(significance_level, lower.tail = FALSE)`, strictly between
  0 and 1. Zero (alpha = 1), a threshold so small that its alpha rounds
  to 1 (below about 1e-16), and one so large that its alpha underflows
  to 0 (above about 38.4) are refused because the Soric false discovery
  rate is undefined there.

- summary_max_samples:

  retained for compatibility and validated as a posterior-sample budget,
  but unused by the default conversion path: direct plotting does not
  compute EDR or missing-study summaries. Use
  `as_zplot(max_samples = ...)` to compute those summaries. The
  plot-density `max_samples` argument is accepted through `...`; it
  defaults to 1000 for `brma.mv` models and 10000 for other models. Use
  `Inf` to use all posterior samples.

- conditioning_depth:

  predictive conditioning depth of the displayed densities. Defaults to
  `"marginal"`; see
  [`as_zplot()`](https://fbartos.github.io/RoBMA/reference/as_zplot.brma.md)
  for details.

- integration_control:

  numerical integration settings created by
  [`set_selection_likelihood_control()`](https://fbartos.github.io/RoBMA/reference/set_selection_likelihood_control.md).
  These control integration within each posterior draw, independently of
  `max_samples`. Factor blocks and supported full-covariance
  product-selection blocks use deterministic Gauss-Hermite quadrature
  where possible, with the existing adaptive randomized quasi-Monte
  Carlo fallback. Density-error checks use the peak within each
  posterior draw and also check the selection normalizer. Bounded CDF
  interpolation, Gaussian-mixture approximation, and omitted
  contributions are charged to the same numerical error budget.

- ...:

  arguments passed to
  [`plot.zplot_brma()`](https://fbartos.github.io/RoBMA/reference/plot.zplot_brma.md).

## Value

`NULL` invisibly for base graphics, or a ggplot2 object.

## Details

For models using the default
[`as_zplot.brma()`](https://fbartos.github.io/RoBMA/reference/as_zplot.brma.md)
method, direct plotting prepares the observed z-statistics and
predictive target without computing EDR or missing-study summary
integrals. It uses the same density computation and numerical checks as
`plot(as_zplot(object))`. A failure confined to a summary integral
therefore does not prevent a plot; density-integration failures still
stop plotting. Use
[`as_zplot()`](https://fbartos.github.io/RoBMA/reference/as_zplot.brma.md)
followed by [`summary()`](https://rdrr.io/r/base/summary.html) to obtain
the summaries.

Custom
[`as_zplot()`](https://fbartos.github.io/RoBMA/reference/as_zplot.brma.md)
methods for subclasses retain their conversion and plotting dispatch,
including any summaries those methods compute.

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
