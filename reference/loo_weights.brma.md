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
stored PSIS object and does not compute LOO. It does not compute weights
across models. For predictive model weights, use
`loo_model_weights(model_1, model_2)`; for a model comparison table, use
`loo_compare(model_1, model_2)`.

When the score integrates local random effects, these weights reweight
global-parameter draws. They must not be applied to full-data sampled
local effects as though those were draws from the deletion posterior;
local-effect calculations require the corresponding retained-data
reconstruction.

## See also

[`loo.brma`](https://fbartos.github.io/RoBMA/reference/loo.brma.md),
[`loo_compare.brma`](https://fbartos.github.io/RoBMA/reference/loo_compare.brma.md),
[`loo_model_weights.brma`](https://fbartos.github.io/RoBMA/reference/loo_model_weights.brma.md)
