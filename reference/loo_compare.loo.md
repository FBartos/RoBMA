# Compare loo Objects Using LOO

Method for comparing RoBMA-targeted `loo` or `waic` objects directly.

## Usage

``` r
# S3 method for class 'loo'
loo_compare(x, ..., unit = "estimate")
```

## Arguments

- x:

  a RoBMA-targeted `loo` or `waic` object (the first model to compare).

- ...:

  additional RoBMA-targeted `loo`/`waic` or `brma` objects to compare.

- unit:

  output/deletion unit used when extracting LOO from brma objects.

## Value

The comparison returned by
[`loo_compare`](https://mc-stan.org/loo/reference/loo_compare.html) for
the stored results: a data frame of class `"compare.loo"` with one row
per model, ordered by expected predictive performance, which prints as
loo prints it. The models are named `model1`, `model2`, ... unless the
call names them.

## See also

[`loo.brma`](https://fbartos.github.io/RoBMA/reference/loo.brma.md),
[`loo_compare`](https://mc-stan.org/loo/reference/loo_compare.html)
