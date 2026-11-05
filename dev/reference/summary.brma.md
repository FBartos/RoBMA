# Summarize brma Object

`summary.brma` creates summary tables for a brma object. For RoBMA
objects, inclusion summaries are printed before parameter estimates.
Random-effect inclusion rows are included in `Component Inclusion`,
immediately before `Publication Bias` when present. Their labels prefix
the corresponding SD quantity with `Heterogeneity:`, retaining component
names where needed. Data-frame exports use the same layout.

## Usage

``` r
# S3 method for class 'brma'
summary(
  object,
  probs = c(0.025, 0.5, 0.975),
  include_mcmc_diagnostics = TRUE,
  standardized_coefficients = FALSE,
  conditional = FALSE,
  logBF = FALSE,
  BF01 = FALSE,
  ...
)

# S3 method for class 'summary.brma'
print(x, ...)

# S3 method for class 'brma'
print(x, ...)
```

## Arguments

- object:

  a fitted brma object

- probs:

  quantiles of the posterior samples to be displayed. Defaults to
  `c(.025, .50, .975)`

- include_mcmc_diagnostics:

  whether to include MCMC diagnostics in the output. Defaults to `TRUE`.

- standardized_coefficients:

  whether to show standardized formula coefficients. Defaults to
  `FALSE`. When set to `TRUE`, meta-regression coefficients and
  random-component summaries are returned on the standardized scale on
  which prior distributions are specified by default (i.e.,
  `standardize_continuous_predictors = TRUE`). A random-effect row whose
  original-scale value the fitted random block cannot represent (for
  example the intercept of a block with factor-specific slopes of a
  standardized predictor) is `NA` in the default original-scale table,
  with a footnote; its value on the fitted scale is shown with
  `standardized_coefficients = TRUE`.

- conditional:

  whether to include conditional estimates for RoBMA product-space
  objects. Defaults to `FALSE`.

- logBF:

  whether to show inclusion Bayes factors on the log scale. Defaults to
  `FALSE`.

- BF01:

  whether to show inverse inclusion Bayes factors. Defaults to `FALSE`.

- ...:

  unused additional arguments. Supplied arguments trigger a warning.

- x:

  a `summary.brma` or fitted `brma` object.

## Value

A list of class `summary.brma` with model name, optional RoBMA inclusion
tables, common estimates, moderator estimates, scale estimates,
publication-bias estimates, and optional conditional estimates.
`known_v_backend` records the known-sampling-covariance backend.
`inclusion_random` contains random-effect inclusion tables;
`selection_model`, `selection_sampling`, and
`selection_sensitivity_diagnostics` retain selection-target metadata and
optional sensitivity diagnostics. The printed form displays the
non-empty tables. Objects created with `only_data = TRUE` return the
resolved data, and objects created with `only_priors = TRUE` return a
`summary.brma` object that prints only the model title, as at 4.0.0, and
holds the resolved prior list in its `priors` element, because neither
carries a posterior. Intercept-only multivariate location models label
their sole coefficient `mu`, consistently with ordinary meta-analysis
models. Once location moderators are present, the coefficient is labeled
`intercept`; an intercept fixed at zero is omitted. Multivariate models
without `mods` or `scale` display location and random estimates together
in one `Estimates` table. With predictors, random estimates join
`Common Estimates`, location coefficients appear in `Meta-Regression`
(or `Location` when `scale` is present), and scale coefficients appear
in `Scale`. Without location moderators, `mu` remains in the common
table. Conditional estimates follow the same layout. Data-frame
conversion follows these displayed sections; the returned list retains
its separate raw tables. Scale tables label the baseline SD as
`exp(intercept)`; its estimates are already exponentiated. Other scale
coefficients remain on the log-SD scale. Multivariate scale rows
identify their targeted random-effect heterogeneity (`tau`), total
heterogeneity (`tau_total`), or mean-variance allocation scale
(`tau_common`). The baseline corresponds to zero non-intercept design
columns under the fitted contrasts and the requested predictor
standardization. The random table reports the quantities aligned with
prior specification; use
[`summary_heterogeneity()`](https://fbartos.github.io/RoBMA/dev/reference/summary_heterogeneity.md)
for aggregate variances and the complete family of deterministic
allocation transforms. This reporting rule is the same for
[`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
and
[`BMA.mv()`](https://fbartos.github.io/RoBMA/dev/reference/BMA.mv.md):
component `tau` values derived from a variance allocation appear only in
[`summary_heterogeneity()`](https://fbartos.github.io/RoBMA/dev/reference/summary_heterogeneity.md).
Optional selection-conditioning comparisons are inspected explicitly
with
[`selection_sensitivity_diagnostics()`](https://fbartos.github.io/RoBMA/dev/reference/selection_sensitivity_diagnostics.md).
Printing a model or its summary does not emit sensitivity notifications.
Fitted objects store their default summary, so
[`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) without further
arguments return it without recomputation;
[`update()`](https://fbartos.github.io/RoBMA/dev/reference/update.brma.md)
and
[`add_selection_sensitivity_diagnostics()`](https://fbartos.github.io/RoBMA/dev/reference/selection_sensitivity_diagnostics.md)
store it again.

## See also

[`brma()`](https://fbartos.github.io/RoBMA/dev/reference/brma.md),
[`brma.glmm()`](https://fbartos.github.io/RoBMA/dev/reference/brma.glmm.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")

  fit <- bPET(
    yi      = yi,
    vi      = vi,
    data    = dat.lehmann2018,
    measure = "SMD",
    seed    = 1,
    silent  = TRUE
  )

  summary(fit)
}
} # }

```
