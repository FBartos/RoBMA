# Bridge Sampling for brma Objects

Extract the marginal likelihood bridge sampling object from a brma
model. The marginal likelihood must first be computed using
[`add_marglik`](https://fbartos.github.io/RoBMA/reference/add_marglik.brma.md).

## Usage

``` r
# S3 method for class 'brma'
bridge_sampler(samples, ...)
```

## Arguments

- samples:

  a brma model object.

- ...:

  additional arguments (currently not used).

## Value

An object of class `"bridge"` as returned by
[`bridge_sampler`](https://rdrr.io/pkg/bridgesampling/man/bridge_sampler.html).

## Details

This function extracts the bridge sampling object that was previously
computed and stored using
[`add_marglik`](https://fbartos.github.io/RoBMA/reference/add_marglik.brma.md).
If the marginal likelihood has not been computed, an error is thrown.
Product-space model-averaging objects (`BMA.norm`, `BMA.glmm`, and
`RoBMA`) do not expose bridge-sampling marginal likelihoods.

The returned object can be used for Bayesian model comparison via
[`bf`](https://rdrr.io/pkg/bridgesampling/man/bf.html) and
[`post_prob`](https://rdrr.io/pkg/bridgesampling/man/post_prob.html).

## See also

[`add_marglik`](https://fbartos.github.io/RoBMA/reference/add_marglik.brma.md),
[`bridge_sampler`](https://rdrr.io/pkg/bridgesampling/man/bridge_sampler.html),
[`logml.brma`](https://fbartos.github.io/RoBMA/reference/logml.brma.md),
[`bf.brma`](https://fbartos.github.io/RoBMA/reference/bf.brma.md),
[`post_prob.brma`](https://fbartos.github.io/RoBMA/reference/post_prob.brma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")
  fit <- brma(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")

  fit <- add_marglik(fit)

  bridge <- bridge_sampler(fit)
  print(bridge)
}
} # }
```
