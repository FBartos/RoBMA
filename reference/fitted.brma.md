# Fitted Values for brma Objects

Extract in-sample fitted values from a fitted `brma` object.

## Usage

``` r
# S3 method for class 'brma'
fitted(
  object,
  unit = "estimate",
  conditioning_depth = "marginal",
  component = "location",
  bias_adjusted = FALSE,
  output_measure = NULL,
  transform = NULL,
  conditional = FALSE,
  ...
)
```

## Arguments

- object:

  a fitted brma object.

- unit:

  output unit. Only `"estimate"` is implemented currently.

- conditioning_depth:

  conditioning depth for location fitted values. `"marginal"` uses fixed
  effects only, `"cluster"` conditions on cluster-level random effects,
  and `"estimate"` conditions on the full estimate-level fitted value.

- component:

  fitted component to return. `"location"` returns location fitted
  values, `"scale"` returns fitted heterogeneity \\\tau_i\\, and `"all"`
  returns both as a named list.

- bias_adjusted:

  whether location fitted values should adjust for publication bias.
  Defaults to `FALSE`.

- output_measure:

  effect-size measure for location/effect predictions. Defaults to the
  fitted measure. Supported conversions are among `"SMD"`, `"COR"`,
  `"ZCOR"`, and `"OR"`; `"RR"`, `"HR"`, `"IRR"`, `"RD"`, and `"GEN"` can
  only be returned on their fitted measure. Use `transform = "EXP"` for
  ratio-scale output from log-scale measures.

- transform:

  optional display transformation. Currently `"EXP"` exponentiates
  log-scale measures `"OR"`, `"RR"`, `"HR"`, and `"IRR"`.

- conditional:

  whether to return fitted values from conditional posterior predictions
  for RoBMA product-space objects.

- ...:

  additional arguments. Currently only `quiet` is honored.

## Value

A numeric vector of fitted values, or a named list with `location` and
`scale` components when `component = "all"`.

## Details

This method is a compact adapter around
[`predict.brma`](https://fbartos.github.io/RoBMA/reference/predict.brma.md).
It summarizes posterior prediction draws by their column means and
returns a base numeric vector, matching the usual
[`fitted`](https://rdrr.io/r/stats/fitted.values.html) contract. Use
[`predict()`](https://rdrr.io/r/stats/predict.html) directly when
posterior draws or intervals are needed.

The default `conditioning_depth = "marginal"` corresponds to
`predict(object, type = "terms")` and matches the usual fitted-value
convention for meta-regression. For normal models,
`conditioning_depth = "estimate"` corresponds to BLUP means for the
observed estimates.

For `component = "all"`, `conditioning_depth`, `output_measure`, and
`transform` apply only to the `location` component. The `scale`
component always returns fitted \\\tau_i\\ values.

## See also

[`predict.brma()`](https://fbartos.github.io/RoBMA/reference/predict.brma.md),
[`residuals.brma()`](https://fbartos.github.io/RoBMA/reference/residuals.brma.md),
[`blup.brma()`](https://fbartos.github.io/RoBMA/reference/blup.brma.md)
