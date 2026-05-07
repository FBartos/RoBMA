# Compare brma Models Using LOO

Compare multiple brma models using LOO-PSIS cross-validation. This is a
convenience wrapper around
[`loo_compare`](https://mc-stan.org/loo/reference/loo_compare.html).

## Usage

``` r
# S3 method for class 'brma'
loo_compare(x, ..., unit = "estimate")
```

## Arguments

- x:

  a brma model object (the first model to compare).

- ...:

  additional brma model objects or `loo` objects to compare.

- unit:

  output/deletion unit used when extracting LOO from brma objects.

## Value

A matrix of class `"compare.loo"` as returned by
[`loo_compare`](https://mc-stan.org/loo/reference/loo_compare.html).

## Details

This function compares models based on their expected out-of-sample
predictive performance (ELPD).

**Important for model comparison:** When comparing models via
[`loo_compare`](https://mc-stan.org/loo/reference/loo_compare.html), the
selection is based on expected out-of-sample predictive performance.
This evaluates how well models predict *new* observations, not how well
they fit the observed data. RoBMA rejects comparisons with different
outcome targets/data, `unit`, or implied `conditioning_depth`.

## See also

[`loo.brma`](https://fbartos.github.io/RoBMA/reference/loo.brma.md),
[`loo_compare`](https://mc-stan.org/loo/reference/loo_compare.html)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")

  fit_bias <- RoBMA(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
  fit_nobias <- BMA(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")

  fit_bias <- add_loo(fit_bias)
  fit_nobias <- add_loo(fit_nobias)

  loo_compare(fit_bias, fit_nobias)
  loo_compare(loo(fit_bias), loo(fit_nobias))
}
} # }
```
