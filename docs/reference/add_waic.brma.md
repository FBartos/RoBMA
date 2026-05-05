# Add WAIC to brma Objects

Compute the Widely Applicable Information Criterion (WAIC) for brma
model objects and store the result in the object.

## Usage

``` r
# S3 method for class 'brma'
add_waic(object, unit = "estimate", ...)
```

## Arguments

- object:

  a brma model object.

- unit:

  output/deletion unit. See
  [`add_loo`](https://fbartos.github.io/RoBMA/reference/add_loo.brma.md);
  the same accepted values and multilevel constraint apply.

- ...:

  additional arguments passed to
  [`waic`](https://mc-stan.org/loo/reference/waic.html).

## Value

The brma object with the WAIC result stored in
`object[["waic"]][[unit]]`.

## Details

WAIC is an alternative to LOO-CV for estimating out-of-sample predictive
accuracy. Like LOO, it evaluates expected predictive performance for new
observations.

In most cases, LOO-PSIS (via
[`add_loo`](https://fbartos.github.io/RoBMA/reference/add_loo.brma.md))
is preferred over WAIC because it provides better estimates and includes
diagnostics (Pareto k values) that indicate when the approximation may
be unreliable.

## See also

[`waic.brma`](https://fbartos.github.io/RoBMA/reference/waic.brma.md),
[`add_loo`](https://fbartos.github.io/RoBMA/reference/add_loo.brma.md),
[`waic`](https://mc-stan.org/loo/reference/waic.html)
