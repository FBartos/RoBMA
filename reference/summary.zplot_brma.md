# Summarize Zplot Results

Creates summary tables for zplot estimates including EDR, Soric FDR, and
estimated missing studies.

## Usage

``` r
# S3 method for class 'zplot_brma'
summary(object, probs = c(0.025, 0.975), ...)
```

## Arguments

- object:

  a zplot_brma object.

- probs:

  quantiles of the posterior distribution to display. Defaults to
  `c(.025, .975)`.

- ...:

  additional arguments (currently unused).

## Value

An object of class `"summary.zplot_brma"` containing the estimates
table.

## Details

The summary includes:

- EDR:

  Expected Discovery Rate - average power of significant studies

- Soric FDR:

  Expected proportion of false discoveries among significant results,
  computed from EDR following Soric (1989)

- Missing N:

  Estimated number of studies suppressed by publication bias, computed
  from the selection model weights

The footer reports the Observed Discovery Rate (ODR) with 95\\
comparison with the model-estimated EDR.

## See also

[`as_zplot.brma()`](https://fbartos.github.io/RoBMA/reference/as_zplot.brma.md),
[`plot.zplot_brma()`](https://fbartos.github.io/RoBMA/reference/plot.zplot_brma.md)
