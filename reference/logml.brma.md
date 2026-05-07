# Log Marginal Likelihood for brma Objects

Extract the log marginal likelihood from a brma model. The marginal
likelihood must first be computed using
[`add_marglik`](https://fbartos.github.io/RoBMA/reference/add_marglik.brma.md).

## Usage

``` r
# S3 method for class 'brma'
logml(x, ...)
```

## Arguments

- x:

  a brma model object.

- ...:

  additional arguments (currently not used).

## Value

A scalar numeric value representing the log marginal likelihood.

## Details

This function extracts the log marginal likelihood from the bridge
sampling object that was previously computed and stored using
[`add_marglik`](https://fbartos.github.io/RoBMA/reference/add_marglik.brma.md).
Product-space model-averaging objects (`BMA.norm`, `BMA.glmm`, and
`RoBMA`) do not expose bridge-sampling marginal likelihoods.

## See also

[`add_marglik`](https://fbartos.github.io/RoBMA/reference/add_marglik.brma.md),
[`bridge_sampler.brma`](https://fbartos.github.io/RoBMA/reference/bridge_sampler.brma.md),
[`bf.brma`](https://fbartos.github.io/RoBMA/reference/bf.brma.md),
[`post_prob.brma`](https://fbartos.github.io/RoBMA/reference/post_prob.brma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")
  fit <- brma(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")

  fit <- add_marglik(fit)

  logml(fit)
}
} # }
```
