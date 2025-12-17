# Compute pooled effect size

`pooled_effect` computes the pooled effect size for a fitted RoBMA.reg
and BiBMA.reg object. Only available for models estimated using the
spike-and-slab algorithm (i.e., `algorithm = "ss"`).

## Usage

``` r
pooled_effect(
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

The meta-regression specification results in the intercept corresponding
to the adjusted effect estimate (i.e., adjusting for the effect of
moderators). In case of moderators inbalance, the adjusted effect
estimate might not be representative of the sample of studies. The
pooled effect size function averages the effect size estimate across the
moderators proportionately to the moderators levels observed in the data
set. Note that there is no Bayes factor test for the presence of the
pooled effect (the summary function provides the adjusted effect and the
test for the presence of the adjusted effect).

The conditional estimate is calculated conditional on the presence of
the adjusted effect (i.e., the intercept).

## See also

[`adjusted_effect()`](reference/adjusted_effect.md)
