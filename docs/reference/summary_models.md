# Summarize Model-Averaged Component Weights

Creates marginal or individual model-weight summaries for RoBMA-class
product-space objects.

## Usage

``` r
summary_models(object, ...)

# S3 method for class 'RoBMA'
summary_models(object, type = "marginal", include_mcmc_diagnostics = TRUE, ...)

# S3 method for class 'summary_models.RoBMA'
print(x, ...)
```

## Arguments

- object:

  a fitted RoBMA-class product-space object, including `RoBMA`,
  `BMA`/`BMA.norm`, and `BMA.glmm`.

- ...:

  additional arguments

- type:

  whether to summarize marginal component prior distributions
  (`"marginal"`) or individual model combinations (`"individual"`).

- include_mcmc_diagnostics:

  whether to include Bayes factor MCMC diagnostics in the output.
  Defaults to `TRUE`.

- x:

  a `summary_models.RoBMA` object.

## Value

A list of class `summary_models.RoBMA` with elements `name` and `type`.
For `type = "marginal"`, element `marginal` contains component tables
with columns such as `prior_prob`, `post_prob`, and `inclusion_BF`. For
`type = "individual"`, element `individual` contains individual model
combinations and posterior probabilities.

## Details

Only mixture-prior components are summarized; non-mixture components are
omitted.

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")
  fit <- RoBMA(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")

  summary_models(fit)
  summary_models(fit, type = "individual")
}
} # }
```
