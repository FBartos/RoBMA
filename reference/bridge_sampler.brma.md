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

An object of class `"bridge"`, or `"bridge_list"` for repeated bridge
estimates, as returned by
[`bridge_sampler`](https://rdrr.io/pkg/bridgesampling/man/bridge_sampler.html).

## Details

This function extracts the raw bridgesampling object underlying the
marginal-likelihood result stored by
[`add_marglik`](https://fbartos.github.io/RoBMA/reference/add_marglik.brma.md).
If the marginal likelihood has not been computed, an error is thrown.
Product-space model-averaging objects (`BMA.norm`, `BMA.glmm`, `BMA.mv`,
`RoBMA`, and `RoBMA.mv`) do not expose bridge-sampling marginal
likelihoods.

Fully fixed models can have a zero-dimensional marginal likelihood
evaluated exactly. Such fits have no bridge-sampling object, so this
method raises a `RoBMA_exact_marglik_no_bridge` error; use
[`logml`](https://rdrr.io/pkg/bridgesampling/man/logml.html) to retrieve
their exact log marginal likelihood. For Bayesian model comparison, pass
fitted `brma` objects directly to
[`bf`](https://rdrr.io/pkg/bridgesampling/man/bf.html) or
[`post_prob`](https://rdrr.io/pkg/bridgesampling/man/post_prob.html).

When a stored marginal-likelihood result has no raw bridge object,
extraction raises `RoBMA_bridge_unavailable`. Its `reason` is
`"exact_zero_dimensional"` for an exact result, or `"missing_upstream"`
when a stored result lacks its raw bridge object. For the latter,
recompute the result with
[`add_marglik`](https://fbartos.github.io/RoBMA/reference/add_marglik.brma.md).

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
