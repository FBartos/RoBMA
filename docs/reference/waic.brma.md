# WAIC for brma Objects

Extract the WAIC object from a brma model object. The WAIC must first be
computed using
[`add_waic`](https://fbartos.github.io/RoBMA/reference/add_waic.brma.md).

## Usage

``` r
# S3 method for class 'brma'
waic(x, unit = "estimate", ...)
```

## Arguments

- x:

  a brma model object.

- unit:

  output/deletion unit. See
  [`add_loo`](https://fbartos.github.io/RoBMA/reference/add_loo.brma.md).

- ...:

  additional arguments (currently unused).

## Value

An object of class `"waic"` as returned by
[`waic`](https://mc-stan.org/loo/reference/waic.html).

## Details

This function extracts the WAIC object that was previously computed and
stored using `object <- add_waic(object, unit = unit)`. If WAIC has not
been computed for the requested unit, an error is thrown.

This is the RoBMA S3 generic and `brma` method. The method is also
registered for [`waic`](https://mc-stan.org/loo/reference/waic.html), so
`loo::waic(fit)` extracts the cached WAIC object for `brma` fits. Use
[`waic`](https://mc-stan.org/loo/reference/waic.html) directly for raw
log-likelihood arrays or matrices.

In most cases, LOO-PSIS (via
[`loo.brma`](https://fbartos.github.io/RoBMA/reference/loo.brma.md)) is
preferred over WAIC because it provides better estimates and includes
diagnostics (Pareto k values) that indicate when the approximation may
be unreliable.

## See also

[`add_waic`](https://fbartos.github.io/RoBMA/reference/add_waic.brma.md),
[`loo.brma`](https://fbartos.github.io/RoBMA/reference/loo.brma.md),
[`waic`](https://mc-stan.org/loo/reference/waic.html)
