# Summarizes heterogeneity of a RoBMA model

Computes the prediction interval, the absolute heterogeneity (tau,
tau^2), and relative measures of heterogeneity (I^2, H^2) for a fitted
RoBMA object.

## Usage

``` r
summary_heterogeneity(
  object,
  type = "ensemble",
  conditional = FALSE,
  output_scale = NULL,
  probs = c(0.025, 0.975),
  short_name = FALSE,
  remove_spike_0 = FALSE
)
```

## Arguments

- object:

  a fitted RoBMA object

- type:

  whether to show the overall RoBMA results (`"ensemble"`) or a detailed
  summary of the individual models (`"individual"`). Can be abbreviated
  to first letters.

- conditional:

  show the conditional estimates (assuming that the alternative is
  true). Defaults to `FALSE`. Only available for `type == "ensemble"`.

- output_scale:

  transform the meta-analytic estimates to a different scale. Defaults
  to `NULL` which returns the same scale as the model was estimated on.

- probs:

  quantiles of the posterior samples to be displayed. Defaults to
  `c(.025, .975)`

- short_name:

  whether priors names should be shortened to the first (couple) of
  letters. Defaults to `FALSE`.

- remove_spike_0:

  whether spike prior distributions with location at zero should be
  omitted from the summary. Defaults to `FALSE`.

## Value

`summary_heterogeneity` returns a list of tables of class
'BayesTools_table'.

## Details

The `conditional` argument allows for computing the conditional
prediction interval based on models assuming the presence of the effect
and the conditional heterogeneity estimates tau, tau^2, I^2, and H^2
assuming the presence of the heterogeneity.

Relative heterogeneity measures (I^2 and H^2) are not available for
BiBMA models.
