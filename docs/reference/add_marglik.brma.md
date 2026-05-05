# Add Marginal Likelihood to brma Objects

Compute the marginal likelihood of a brma model using bridge sampling
and store the result in the object.

## Usage

``` r
# S3 method for class 'brma'
add_marglik(object, ...)
```

## Arguments

- object:

  a brma model object.

- ...:

  additional arguments (currently not used).

## Value

The brma object with the marginal likelihood result stored in
`object[["marglik"]]`.

## Details

The marginal likelihood is computed using the `bridgesampling` package
via the
[`BayesTools::JAGS_bridgesampling`](https://fbartos.github.io/BayesTools/reference/JAGS_bridgesampling.html)
wrapper. The result is stored in the object and can be extracted using
[`bridge_sampler.brma`](https://fbartos.github.io/RoBMA/reference/bridge_sampler.brma.md).

Product-space model-averaging objects (`BMA.norm`, `BMA.glmm`, and
`RoBMA`) do not expose a bridge-sampling marginal likelihood; use
predictive comparison methods such as
[`loo.brma`](https://fbartos.github.io/RoBMA/reference/loo.brma.md)
instead.

## See also

[`bridge_sampler.brma`](https://fbartos.github.io/RoBMA/reference/bridge_sampler.brma.md),
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

  logml(fit)
}
} # }
```
