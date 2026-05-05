# Summarize fitted RoBMA object

`summary.RoBMA` creates summary tables for a RoBMA object.

## Usage

``` r
# S3 method for class 'RoBMA'
summary(
  object,
  type = "ensemble",
  conditional = FALSE,
  output_scale = NULL,
  probs = c(0.025, 0.975),
  logBF = FALSE,
  BF01 = FALSE,
  short_name = FALSE,
  remove_spike_0 = FALSE,
  ...
)
```

## Arguments

- object:

  a fitted RoBMA object

- type:

  whether to show the overall RoBMA results (`"ensemble"`), an overview
  of the individual models (`"models"`), an overview of the individual
  models MCMC diagnostics (`"diagnostics"`), or a detailed summary of
  the individual models (`"individual"`). Can be abbreviated to first
  letters.

- conditional:

  show the conditional estimates (assuming that the alternative is
  true). Defaults to `FALSE`. Only available for `type == "ensemble"`.

- output_scale:

  transform the meta-analytic estimates to a different scale. Defaults
  to `NULL` which returns the same scale as the model was estimated on.

- probs:

  quantiles of the posterior samples to be displayed. Defaults to
  `c(.025, .975)`

- logBF:

  show log of Bayes factors. Defaults to `FALSE`.

- BF01:

  show Bayes factors in support of the null hypotheses. Defaults to
  `FALSE`.

- short_name:

  whether priors names should be shortened to the first (couple) of
  letters. Defaults to `FALSE`.

- remove_spike_0:

  whether spike prior distributions with location at zero should be
  omitted from the summary. Defaults to `FALSE`.

- ...:

  additional arguments

## Value

`summary.RoBMA` returns a list of tables of class 'BayesTools_table'.

## Note

See
[`diagnostics()`](https://fbartos.github.io/RoBMA/reference/diagnostics.md)
for visual convergence checks of the individual models.

## See also

[`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md),
[`diagnostics()`](https://fbartos.github.io/RoBMA/reference/diagnostics.md),
[`check_RoBMA()`](https://fbartos.github.io/RoBMA/reference/check_RoBMA.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# using the example data from Anderson et al. 2010 and fitting the default model
# (note that the model can take a while to fit)
fit <- RoBMA(r = Anderson2010$r, n = Anderson2010$n, study_names = Anderson2010$labels)

# summary can provide many details about the model
summary(fit)

# estimates from the conditional models can be obtained with
summary(fit, conditional = TRUE)

# overview of the models and their prior and posterior probability, marginal likelihood,
# and inclusion Bayes factor can be obtained with
summary(fit, type = "models")

# diagnostics overview, containing the maximum R-hat, minimum ESS, maximum MCMC error, and
# maximum MCMC error / sd across parameters for each individual model can be obtained with
summary(fit, type = "diagnostics")

# summary of individual models and their parameters can be further obtained by
summary(fit, type = "individual")
} # }
```
