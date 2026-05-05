# Compute adjusted effect size

`adjusted_effect` computes the adjusted effect size for a fitted
RoBMA.reg and BiBMA.reg object.

## Usage

``` r
adjusted_effect(
  object,
  conditional = FALSE,
  output_scale = NULL,
  probs = c(0.025, 0.975),
  ...
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

- ...:

  additional arguments

## Value

`pooled_effect` returns a list of tables of class 'BayesTools_table'.

## Details

Non-default meta-regression specification (i.e., using treatment
contrasts for predictors) might results in the intercept corresponding
to the effect estimate in the baseline group. (i.e., adjusting for the
effect of moderators). The adjusted effect size function averages the
effect size estimate across the moderators levels. Note that there is no
Bayes factor test for the presence of the adjusted effect (the summary
function provides the effect estimate in the baseline group and the test
for the presence of the effect in the baseline group if a treatment
contrasts are specified).

The conditional estimate is calculated conditional on the presence of
the baseline group effect (i.e., the intercept).

## See also

[`pooled_effect()`](https://fbartos.github.io/RoBMA/reference/pooled_effect.md)
