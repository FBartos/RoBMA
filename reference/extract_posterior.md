# Extract Posterior Samples from a RoBMA Model

Extracts posterior samples for a specified parameter from a RoBMA model
object.

## Usage

``` r
extract_posterior(
  x,
  parameter = "mu",
  conditional = FALSE,
  output_scale = NULL,
  ...
)
```

## Arguments

- x:

  a fitted RoBMA object

- parameter:

  a parameter for which posterior samples should be extracted. Defaults
  to `"mu"` (for the effect size). The additional options are `"tau"`
  (for the heterogeneity), `"weightfunction"` (for the estimated
  weightfunction), or `"PET"` and `"PEESE"` (for the PET-PEESE
  coefficients).

- conditional:

  whether conditional estimates should be extracted. Defaults to `FALSE`
  which extracts the model-averaged estimates. Note that both
  `"weightfunction"` and `"PET-PEESE"` are always ignoring the other
  type of publication bias adjustment.

- output_scale:

  transform the effect sizes and the meta-analytic effect size estimate
  to a different scale. Defaults to `NULL` which returns the same scale
  as the model was estimated on.

- ...:

  additional arguments passed to the method.

## Value

A matrix containing the posterior samples for the specified parameter.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming 'fit' is a fitted RoBMA model:
posterior_mu <- extract_posterior(fit, parameter = "mu")
} # }
```
