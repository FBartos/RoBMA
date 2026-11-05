# Summarize brma_samples Object

Creates and returns a summary table of posterior samples using
[`BayesTools::ensemble_estimates_table`](https://fbartos.github.io/BayesTools/reference/BayesTools_ensemble_tables.html).

## Usage

``` r
# S3 method for class 'brma_samples'
summary(object, probs = NULL, ...)
```

## Arguments

- object:

  a `brma_samples` object

- probs:

  quantiles for credible intervals. If `NULL`, uses the default stored
  in the object (typically `c(.025, .975)`)

- ...:

  additional arguments passed to
  [`BayesTools::ensemble_estimates_table`](https://fbartos.github.io/BayesTools/reference/BayesTools_ensemble_tables.html)

## Value

A `BayesTools_table` object containing the summary statistics.
