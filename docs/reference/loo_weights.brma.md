# Extract Normalized PSIS Weights from brma Object

Extract the normalized Pareto smoothed importance sampling (PSIS)
weights from a brma model object.

## Usage

``` r
# S3 method for class 'brma'
loo_weights(object, unit = "estimate", ...)
```

## Arguments

- object:

  a brma model object.

- unit:

  output/deletion unit. See
  [`add_loo`](https://fbartos.github.io/RoBMA/reference/add_loo.brma.md).

- ...:

  currently unused.

## Value

An `S x K` matrix for estimate-unit LOO, or `S x G` matrix for
cluster-unit LOO, with posterior samples in rows and LOO targets in
columns. Columns are normalized to sum to one.

## Details

LOO must first be computed with
`object <- add_loo(object, unit = unit)`. This method extracts the
stored PSIS object and does not compute LOO.
