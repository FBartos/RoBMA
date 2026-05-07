# Summarize brma Object

`summary.brma` creates summary tables for a brma object. For RoBMA
objects, inclusion summaries are printed before parameter estimates.

## Usage

``` r
# S3 method for class 'brma'
summary(
  object,
  probs = c(0.025, 0.5, 0.975),
  include_mcmc_diagnostics = TRUE,
  standardized_coefficients = FALSE,
  conditional = FALSE,
  logBF = FALSE,
  BF01 = FALSE,
  ...
)

# S3 method for class 'summary.brma'
print(x, ...)

# S3 method for class 'brma'
print(x, ...)
```

## Arguments

- object:

  a fitted brma object

- probs:

  quantiles of the posterior samples to be displayed. Defaults to
  `c(.025, .50, .975)`

- include_mcmc_diagnostics:

  whether to include MCMC diagnostics in the output. Defaults to `TRUE`.

- standardized_coefficients:

  whether to show standardized meta-regression coefficients. Defaults to
  `FALSE`. When set to `TRUE`, standardized meta-regression coefficients
  are returned for the intercept and continuous predictors. These
  coefficients correspond to the standardized scale on which prior
  distributions are specified by default (i.e.,
  `standardize_continuous_predictors = TRUE`).

- conditional:

  whether to include conditional estimates for RoBMA product-space
  objects. Defaults to `FALSE`.

- logBF:

  whether to show inclusion Bayes factors on the log scale. Defaults to
  `FALSE`.

- BF01:

  whether to show inverse inclusion Bayes factors. Defaults to `FALSE`.

- ...:

  additional arguments

- x:

  a `summary.brma` or fitted `brma` object.

## Value

A list of class `summary.brma` with model name, optional RoBMA inclusion
tables, common estimates, moderator estimates, scale estimates,
publication-bias estimates, and optional conditional estimates. The
printed form displays the non-empty tables.

## See also

[`brma()`](https://fbartos.github.io/RoBMA/reference/brma.md),
[`brma.glmm()`](https://fbartos.github.io/RoBMA/reference/brma.glmm.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")

  fit <- bPET(
    yi      = yi,
    vi      = vi,
    data    = dat.lehmann2018,
    measure = "SMD",
    seed    = 1,
    silent  = TRUE
  )

  summary(fit)
}
} # }

```
