# Print summary.brma_samples Object

Prints a summary table of posterior samples using
[`BayesTools::ensemble_estimates_table`](https://fbartos.github.io/BayesTools/reference/BayesTools_ensemble_tables.html).
Returns the summary table invisibly.

## Usage

``` r
# S3 method for class 'summary.brma_samples'
print(x, probs = NULL, ...)
```

## Arguments

- x:

  a `summary.brma_samples` object

- probs:

  quantiles for credible intervals. If `NULL`, uses the default stored
  in the object (typically `c(.025, .975)`)

- ...:

  additional arguments passed to
  [`BayesTools::ensemble_estimates_table`](https://fbartos.github.io/BayesTools/reference/BayesTools_ensemble_tables.html)

## Value

Returns the summary table invisibly.
