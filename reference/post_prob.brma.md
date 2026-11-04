# Posterior Model Probabilities for brma Objects

Compute posterior model probabilities from marginal likelihoods of brma
models.

## Usage

``` r
# S3 method for class 'brma'
post_prob(x, ..., prior_prob = NULL, model_names = NULL)
```

## Arguments

- x:

  a brma model object.

- ...:

  additional brma model objects.

- prior_prob:

  numeric vector with prior model probabilities or weights. Values must
  be finite, nonnegative, have the same length as the retained models,
  and have a positive total. If omitted, a uniform prior is used. Finite
  comparisons combine the original weights with marginal likelihoods in
  log space and normalize the resulting probabilities. An exactly zero
  prior weight excludes that model.

- model_names:

  character vector with model names. If `NULL` (the default), names will
  be derived from deparsing the call.

## Value

A named numeric vector with posterior model probabilities for single
estimates, or a repetition-by-model numeric matrix when any model has
repeated estimates. Each successfully evaluated comparison sums to one.

## Details

The marginal likelihoods must first be computed using
[`add_marglik`](https://fbartos.github.io/RoBMA/reference/add_marglik.brma.md).
`x` and at least one additional `brma` model must be supplied.
Non-`brma` objects in `...` are ignored with a warning. All retained
models must be fitted to the same outcome target/data, including outcome
type and, when present, likelihood weights. Cluster identifiers define
model structure and are not part of the outcome-data identity. A known
sampling covariance without correlated blocks (e.g., `V = diag(sei^2)`
in [`brma.mv`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)) is
the outcome data of its standard errors, whereas a correlated `V` is a
different target.

If any model contains repeated bridge estimates, comparisons retain
those repetitions regardless of model order. Models with fewer
repetitions are recycled to the largest repetition count with a warning,
following bridgesampling. Standardized marginal-likelihood results must
be finite. Legacy raw `"bridge_list"` results retain non-finite
repetitions and use upstream comparison behavior through bridgesampling,
with normalized prior probabilities and its missing-value propagation,
warnings, and errors.

## See also

[`add_marglik`](https://fbartos.github.io/RoBMA/reference/add_marglik.brma.md),
[`bridge_sampler.brma`](https://fbartos.github.io/RoBMA/reference/bridge_sampler.brma.md),
[`bf.brma`](https://fbartos.github.io/RoBMA/reference/bf.brma.md),
[`logml.brma`](https://fbartos.github.io/RoBMA/reference/logml.brma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")
  fit1 <- brma(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
  fit2 <- brma(
    yi           = yi,
    vi           = vi,
    data         = dat.lehmann2018,
    measure      = "SMD",
    prior_effect = FALSE
  )

  fit1 <- add_marglik(fit1)
  fit2 <- add_marglik(fit2)

  post_prob(fit1, fit2)
}
} # }
```
