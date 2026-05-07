# Plot Estimated Marginal Means

Plots estimated marginal means stored in a `marginal_means.brma` object
using
[`BayesTools::plot_marginal()`](https://fbartos.github.io/BayesTools/reference/plot_marginal.html).

## Usage

``` r
# S3 method for class 'marginal_means.brma'
plot(
  x,
  parameter,
  type = NULL,
  prior = FALSE,
  plot_type = "base",
  dots_prior = NULL,
  output_measure = NULL,
  transform = NULL,
  ...
)
```

## Arguments

- x:

  a `marginal_means.brma` object.

- parameter:

  moderator term to plot. Use the original term name, for example
  `"measure"`, `"intercept"` for the intercept when available, `"mu"` as
  an intercept alias, or the internal parameter name, for example
  `"mu_measure"`.

- type:

  for RoBMA product-space objects, whether to plot model-averaged
  (`"averaged"`) or conditional (`"conditional"`) marginal means.
  Defaults to `"averaged"` and is available only for RoBMA marginal
  means.

- prior:

  whether the marginal prior distribution should be added to the plot.
  Defaults to `FALSE`.

- plot_type:

  whether to use base R graphics (`"base"`) or ggplot2 (`"ggplot"`).
  Defaults to `"base"`.

- dots_prior:

  list of additional graphical arguments passed to the prior plotting
  function.

- output_measure:

  effect-size measure for location/effect predictions. Defaults to the
  fitted measure. Supported conversions are among `"SMD"`, `"COR"`,
  `"ZCOR"`, and `"OR"`; `"RR"`, `"HR"`, `"IRR"`, `"RD"`, and `"GEN"` can
  only be returned on their fitted measure. Use `transform = "EXP"` for
  ratio-scale output from log-scale measures.

- transform:

  optional display transformation. Currently `"EXP"` exponentiates
  log-scale measures `"OR"`, `"RR"`, `"HR"`, and `"IRR"`.

- ...:

  additional graphical arguments passed to
  [`BayesTools::plot_marginal()`](https://fbartos.github.io/BayesTools/reference/plot_marginal.html).

## Value

`NULL` invisibly if `plot_type = "base"` or a ggplot object if
`plot_type = "ggplot"`.
