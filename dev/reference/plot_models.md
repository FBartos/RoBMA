# Plot posterior model probabilities

Display the posterior probabilities of the individual models in a fitted
model-averaged meta-analysis. Each column is one full model and its
width is its posterior probability. Rows describe the model's fitted
components.

## Usage

``` r
plot_models(x, ...)

# S3 method for class 'RoBMA'
plot_models(
  x,
  parameter = NULL,
  component = "auto",
  factor_view = "levels",
  mixture_view = "inclusion",
  group_labels = FALSE,
  labels = NULL,
  decreasing = TRUE,
  col = NULL,
  col_bias = NULL,
  legend = FALSE,
  legend_position = "bottom",
  plot_type = "base",
  ...
)
```

## Arguments

- x:

  a fitted model-averaged object, such as a fit from
  [`RoBMA()`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.md) or
  [`BMA()`](https://fbartos.github.io/RoBMA/dev/reference/BMA.md).

- ...:

  additional arguments for methods.

- parameter:

  an ordered character vector of model component, term, or factor-level
  selectors. Use `NULL` to show the model components, including
  random-effect rows, followed by the location and scale terms.
  Random-effect inclusion rows use the labels printed by
  [`summary.brma()`](https://fbartos.github.io/RoBMA/dev/reference/summary.brma.md)
  and appear before publication bias. Use qualified selectors such as
  `"mu_year"`, `"log_tau_year"`, or `"mu_delivery[online]"` to
  disambiguate terms or levels.

- component:

  component to select from: `"auto"`, `"mods"` (also `"location"`),
  `"scale"`, `"bias"`, or `"random"`.

- factor_view:

  `"levels"` shows a row per fitted factor level and `"terms"` shows one
  neutral inclusion row per factor term. Supply a single setting or
  named overrides for individual terms. Unspecified terms use
  `"levels"`.

- mixture_view:

  `"inclusion"` displays each member's declared null or alternative
  status. `"model"` shows one categorical row for the mixture.
  `"individual"` shows one row for each declared member, including null
  members. Supply a single setting or named overrides for individual
  terms. Unspecified terms use `"inclusion"`.

- group_labels:

  whether to display group headings beside the rows.

- labels:

  named character vector of row label overrides. Use stable selectors
  such as `"mu"`, `"mu_year"`, `"mu_delivery[online]"`, or `"bias[PET]"`
  (when the mixture member is named `PET`). Overrides also update the
  mixture labels in the legend.

- decreasing:

  whether to draw models in decreasing posterior probability order.

- col:

  named colors for `positive`, `negative`, `included`, and `excluded`.
  Defaults are blue, red, black, and white, respectively.

- col_bias:

  named colors for `bias` and `no_bias` or for publication-bias mixture
  members, using their names from
  [`summary_models()`](https://fbartos.github.io/RoBMA/dev/reference/summary_models.md)
  or qualified member selectors of the form `"bias[PET]"`.
  Member-specific colors override bias and no-bias colors in categorical
  and individual views. Use qualified selectors for members named `bias`
  or `no_bias`. This argument requires a displayed publication-bias
  component.

- legend:

  whether to show the color legend.

- legend_position:

  location of the legend, `"bottom"` or `"top"`.

- plot_type:

  `"base"` for base graphics or `"ggplot"` for ggplot2.

## Value

Invisibly `NULL` for base graphics or a ggplot object for ggplot2. The
ggplot object has an integer `n_rows` attribute giving the number of
displayed rows after selection and expansion, excluding the group
headings and legend.

## Details

The column widths come from the same joint indicator draws as
`summary_models(x, type = "individual")`. Choosing rows changes the
displayed information while keeping the fitted models and their widths.
Ungated coefficients can appear alongside the model components. Select
the main components by `"mu"`, `"tau"`, `"bias"`, and `"rho"` for an
available heterogeneity allocation mixture. Use qualified selectors for
moderators with those names.

Inclusion uses the null or alternative assignment of the active prior
member, or the fitted inclusion gate. In the inclusion view, directional
colors indicate the sign of the posterior median within that full model
on the fitted predictor scale. If the predictors were centered or
standardized during fitting, intercept and moderator directions refer to
that parameterization. Heterogeneity, heterogeneity allocation,
publication bias, random-effect gates, whole factors, and coefficients
with a median of zero use neutral inclusion colors. Factor level rows
share the factor's inclusion gate. Treatment contrasts show coefficients
relative to the reference level and omit the reference row.
Mean-difference and orthonormal contrasts show deviations from the
equally weighted mean of the factor levels.

## See also

[`summary_models()`](https://fbartos.github.io/RoBMA/dev/reference/summary_models.md),
[`plot.brma()`](https://fbartos.github.io/RoBMA/dev/reference/plot.brma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- RoBMA(yi = d, sei = se, data = Bem2011, measure = "SMD", seed = 1)
p <- plot_models(fit, parameter = c("mu", "bias"),
                 mixture_view = c(bias = "individual"), plot_type = "ggplot")
attr(p, "n_rows")
p
} # }
```
