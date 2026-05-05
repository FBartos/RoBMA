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

A matrix of class `"compare.loo"` as returned by
[`loo_compare`](https://mc-stan.org/loo/reference/loo_compare.html).

## See also

[`loo.brma`](https://fbartos.github.io/RoBMA/reference/loo.brma.md),
[`loo_compare`](https://mc-stan.org/loo/reference/loo_compare.html)
