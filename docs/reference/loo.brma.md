# LOO-PSIS for brma Objects

Extract the LOO-PSIS object from a brma model object. The LOO must first
be computed using
[`add_loo`](https://fbartos.github.io/RoBMA/reference/add_loo.brma.md).

## Usage

``` r
# S3 method for class 'brma'
loo(x, unit = "estimate", ...)
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

An object of class `c("psis_loo", "loo")` as returned by
[`loo`](https://mc-stan.org/loo/reference/loo.html).

## Details

This function extracts the LOO object that was previously computed and
stored using `object <- add_loo(object, unit = unit)`. If LOO has not
been computed for the requested unit, an error is thrown.

This is the RoBMA S3 generic and `brma` method. Use
[`loo`](https://mc-stan.org/loo/reference/loo.html) directly for raw
log-likelihood arrays or matrices.

## See also

[`add_loo`](https://fbartos.github.io/RoBMA/reference/add_loo.brma.md),
[`loo`](https://mc-stan.org/loo/reference/loo.html),
[`loo_compare`](https://mc-stan.org/loo/reference/loo_compare.html),
[`pareto_k_ids`](https://mc-stan.org/loo/reference/pareto-k-diagnostic.html)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")
  fit <- bPET(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
  fit <- add_loo(fit)

  loo_fit <- loo(fit)
  print(loo_fit)
}
} # }
```
