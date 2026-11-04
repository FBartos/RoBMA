# Ordered Factor Prior

Create priors on ordered factor effects by separating a scalar total
effect from its allocation across ordered-level increments.

## Usage

``` r
prior_ordered(
  total,
  allocation = NULL,
  contrast = "cumulative",
  id = NULL,
  prior_weights = 1
)
```

## Arguments

- total:

  scalar prior for the sum of ordered increments. With `"cumulative"`,
  this is the first-to-last level difference; with
  `"cumulative_levels"`, it is the last level's effect relative to the
  zero baseline preceding the first level.

- allocation:

  allocation prior or fixed split. `NULL` creates a late-bound flat
  Dirichlet allocation with dimension determined by the formula term.
  Numeric vectors specify fixed splits. Dirichlet priors specify random
  allocations. For terms with multiple ordered factors, use a named
  list.

- contrast:

  ordered contrast. `"cumulative"` uses `nlevels(f) - 1` increments and
  sets the first level effect to zero. `"cumulative_levels"` uses
  `nlevels(f)` increments and gives the first level a non-zero
  allocation share.

- id:

  optional allocation-sharing id. Equal ids share compatible allocation
  weights, not the total-effect prior. For terms with multiple ordered
  factors, use a named character vector.

- prior_weights:

  numeric prior model weight.

## Value

An object inheriting from `prior`.

## Details

This is RoBMA's re-export of
[`BayesTools::prior_ordered()`](https://fbartos.github.io/BayesTools/reference/prior_ordered.html).
Spike-and-slab or mixture behavior belongs on `total`, e.g.,
`prior_ordered(prior_spike_and_slab(...))`. Unsupported ordered
combinations raise typed `BayesTools_ordered_unavailable` conditions;
constructor support does not imply support for every post-fit operation.

## See also

[`prior_factor`](https://fbartos.github.io/RoBMA/reference/prior_factor.md)
