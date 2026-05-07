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
  and have a positive total. If omitted, a uniform prior is used.
  Supplied values are normalized internally.

- model_names:

  character vector with model names. If `NULL` (the default), names will
  be derived from deparsing the call.

## Value

A named numeric vector with posterior model probabilities (i.e., which
sum to one).

## Details

The marginal likelihoods must first be computed using
[`add_marglik`](https://fbartos.github.io/RoBMA/reference/add_marglik.brma.md).
`x` and at least one additional `brma` model must be supplied.
Non-`brma` objects in `...` are ignored with a warning. All retained
models must be fitted to the same outcome target/data, including outcome
type and, when present, weights and cluster identifiers.

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
