# Compute estimated true effect sizes

`true_effects` computes the estimated true effect size for a fitted
RoBMA object. These estimates correspond to the frequentist "Best Linear
Unbiased Predictions (BLUPs)"
([blup](https://wviechtb.github.io/metafor/reference/blup.html)). Only
available for normal-normal models estimated using the spike-and-slab
algorithm (i.e., `algorithm = "ss"`).

## Usage

``` r
true_effects(
  object,
  conditional = FALSE,
  output_scale = NULL,
  probs = c(0.025, 0.975),
  as_samples = FALSE
)
```

## Arguments

- object:

  a fitted RoBMA object

- conditional:

  show the conditional estimates (assuming that the alternative is
  true). Defaults to `FALSE`. Only available for `type == "ensemble"`.

- output_scale:

  transform the meta-analytic estimates to a different scale. Defaults
  to `NULL` which returns the same scale as the model was estimated on.

- probs:

  quantiles of the posterior samples to be displayed. Defaults to
  `c(.025, .975)`

- as_samples:

  whether posterior samples instead of a summary table should be
  returned. Defaults to `FALSE`.

## Value

`pooled_effect` returns a list of tables of class 'BayesTools_table'.

## Details

The true effect size estimates are computed under the normal likelihood.
As such, they might not approximate the true effect size estimates under
selection models well.

The conditional estimate is calculated conditional on the presence of
the effect (in meta-analysis) or the intercept (in meta-regression).
