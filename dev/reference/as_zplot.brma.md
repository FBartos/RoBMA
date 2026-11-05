# Transform brma Object to Zplot

Transforms an estimated brma model into a zplot object that can be
summarized and plotted to assess replicability.

## Usage

``` r
# S3 method for class 'brma'
as_zplot(
  object,
  significance_level = stats::qnorm(0.975),
  max_samples = 10000,
  conditioning_depth = "marginal",
  integration_control = set_selection_likelihood_control(),
  ...
)
```

## Arguments

- object:

  a normal-outcome `brma` object.

- significance_level:

  positive z-value threshold for significance. Defaults to
  `qnorm(0.975)` (two-sided alpha = 0.05). The threshold must be finite
  and its two-sided alpha,
  `2 * pnorm(significance_level, lower.tail = FALSE)`, strictly between
  0 and 1. Zero (alpha = 1), a threshold so small that its alpha rounds
  to 1 (below about 1e-16), and one so large that its alpha underflows
  to 0 (above about 38.4) are refused because the Soric false discovery
  rate is undefined there.

- max_samples:

  maximum number of posterior samples for the EDR and missing-study
  summaries. Defaults to 10000. Use `Inf` to use all posterior samples.
  The density budget is supplied separately to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) or
  [`lines()`](https://rdrr.io/r/graphics/lines.html) and defaults to
  1000 for `brma.mv` models and 10000 for other models.

- conditioning_depth:

  predictive conditioning depth. Options are `"marginal"` (default),
  which integrates all newly realized random effects; `"cluster"`, which
  retains fitted cluster effects and integrates within-cluster
  heterogeneity; and `"estimate"`, which retains each fitted latent true
  effect, integrates its posterior conditional uncertainty, and adds new
  sampling error. Cluster conditioning is available only for multilevel
  models. Selection models currently support only `"marginal"`
  conditioning.

- integration_control:

  numerical integration settings created by
  [`set_selection_likelihood_control()`](https://fbartos.github.io/RoBMA/dev/reference/set_selection_likelihood_control.md).
  These control integration within each posterior draw, independently of
  `max_samples`. Factor blocks and supported full-covariance
  product-selection blocks use deterministic Gauss-Hermite quadrature
  where possible, with the existing adaptive randomized quasi-Monte
  Carlo fallback. Density-error checks use the peak within each
  posterior draw and also check the selection normalizer. Bounded CDF
  interpolation, Gaussian-mixture approximation, and omitted
  contributions are charged to the same numerical error budget.

- ...:

  additional arguments (currently unused).

## Value

The input object with added class `"zplot_brma"` and a new `zplot` list
component containing:

- estimates:

  a list with posterior samples for `EDR` and extrapolation `weights`,
  when a missing-count convention applies

- data:

  a list with `significance_level`, `conditioning_depth`, stored
  `integration_control`, observed `z`-statistics, `N_significant`, and
  `N_observed`

- target:

  the fitted reporting specification, reference distribution, source and
  publication identities, and missing-count applicability

## Details

Zplot is the pushforward into z-space of the design-conditional
posterior predictive distribution. It retains the empirical mix of
observed standard errors and moderators, evaluates one row-specific
predictive distribution for each estimate, transforms each through \\Z =
Y / SE\\, and averages the resulting densities with equal weight per
estimate.

The default marginal target integrates every latent effect that would be
newly realized in a replicated literature. For a specialized multilevel
model fitted with `cluster`, the predictive variance before sampling
error is the sum of between- and within-cluster heterogeneity. Cluster
conditioning instead retains the fitted cluster effects and includes
only within-cluster heterogeneity. Estimate conditioning uses the
conditional true-effect means returned by `predict(type = "blup")` and
additionally integrates the posterior conditional variance of those
latent effects before adding new sampling error. It therefore represents
a posterior predictive replication of the same uncertain latent effects
rather than a BLUP-mean plug-in distribution. Thus the marginal zplot
collapses the row-specific slices of the marginal funnel target onto the
common significance scale. Unlike LOO diagnostics, it is a
full-posterior predictive projection and does not leave observations or
clusters out when estimating model parameters.

The Expected Discovery Rate (EDR) is the posterior average probability
that a new outcome from the bias-adjusted reference distribution is
statistically significant at the supplied threshold. It is therefore not
an exact-repeat estimand unless `conditioning_depth = "estimate"` is
requested.

The EDR uses the bias-adjusted reference distribution: selection and
PET/PEESE offsets are removed. For multivariate selection models, models
retaining estimate-level or sampling-error sources or integrating
applicable contextual sources, and best-rule selection models, this
reference is the normalized pre-selection Gaussian marginal, integrating
all newly realized contexts. Its integral is one and EDR is its
significance probability. Relative selection weights do not identify an
absolute publication probability or missing count, so `Missing N` is
unavailable for these targets.

The fitted selected curve preserves the original publication event and
the fitted model's normalization order. It integrates the model's
integrated sources inside each full-event normalizer and averages
conditionally normalized projections over newly realized retained
contexts. A best rule still selects the full original event when
displaying one estimate. Supported full-covariance product blocks reuse
numerator and normalizer calculations while refining their quadrature
separately. Computation by posterior row and covariance block preserves
the full publication event.

Numerical diagnostics distinguish raw normalizer-mass ratios from bounds
on omitted contributions. A raw mass ratio does not certify that the
compressed displayed curve integrates to one. An omitted-mass bound
covers only the omitted components, not the total integrated error of
the compressed curve. Retained component weights are not renormalized to
hide omissions.

Independent univariate product-weight models and the released
conditional cluster display retain their reference-bin-normalized
inverse-weight extrapolation convention and corresponding `Missing N`
summary. Product selection does not use publication groups; its
covariance dependencies and source settings determine the reference
distribution and diagnostic. Best-p-value selection preserves the
resolved publication groups. A literal count interpretation additionally
treats the reference bin as fully observed and requires a
reporting/stopping design relating attempted to observed estimates.
These assumptions are not supplied by the fitted relative-weight
likelihood, which conditions on the reported vector and does not
identify an absolute publication probability or finite missing pool.
Valid above-one relative weights can give a negative diagnostic; such a
value is not a literal missing count. Restricting a plot with `from` or
`to` does not restrict the summary integrals.

Zplot diagnostics are available only for normal outcome models. GLMM
objects are rejected because their raw likelihood is on a count scale
while zplot diagnostics require observed effect-size z-statistics with
standard errors. For correlated known-`V`
[`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
models, zplot is a scalar display using \\z_i = y_i / se_i\\. The
histogram does not display dependence, but the fitted selected marginals
retain the full covariance. Numerical integration can be substantially
more expensive for correlated selection than for ordinary Gaussian
models. Density plotting uses its own `max_samples` argument,
independently of the summaries stored here.

The resulting object retains all original brma properties while adding
zplot results, enabling both standard meta-analytic summaries and zplot
diagnostics on the same object.

## See also

[`summary.zplot_brma()`](https://fbartos.github.io/RoBMA/dev/reference/summary.zplot_brma.md),
[`plot.zplot_brma()`](https://fbartos.github.io/RoBMA/dev/reference/plot.zplot_brma.md),
[`hist.zplot_brma()`](https://fbartos.github.io/RoBMA/dev/reference/hist.zplot_brma.md),
[`predict.brma()`](https://fbartos.github.io/RoBMA/dev/reference/predict.brma.md),
[`funnel.brma()`](https://fbartos.github.io/RoBMA/dev/reference/funnel.md),
[`regplot.brma()`](https://fbartos.github.io/RoBMA/dev/reference/regplot.md)

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
