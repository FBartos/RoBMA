# Unsupported Information Criteria for brma Objects

`brma` objects expose posterior pointwise log-likelihood draws
([`logLik.brma`](https://fbartos.github.io/RoBMA/dev/reference/logLik.brma.md),
[`log_lik`](https://fbartos.github.io/RoBMA/dev/reference/log_lik.md)),
not a maximized scalar log-likelihood. AIC and BIC are therefore
undefined and these methods always stop with an explicit error.

## Usage

``` r
# S3 method for class 'brma'
AIC(object, ..., k = 2)

# S3 method for class 'brma'
BIC(object, ...)
```

## Arguments

- object:

  a brma model object.

- ...:

  ignored.

- k:

  ignored; retained for the
  [`stats::AIC()`](https://rdrr.io/r/stats/AIC.html) method signature.

## Value

These methods do not return a value.

## See also

[`log_lik`](https://fbartos.github.io/RoBMA/dev/reference/log_lik.md),
[`logLik.brma`](https://fbartos.github.io/RoBMA/dev/reference/logLik.brma.md),
[`loo.brma`](https://fbartos.github.io/RoBMA/dev/reference/loo.brma.md),
[`waic.brma`](https://fbartos.github.io/RoBMA/dev/reference/waic.brma.md)
