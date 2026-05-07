# Plots Weight Function of brma Object

`plot_weightfunction.brma` visualizes the posterior (and optionally
prior) publication-bias weight function of a brma object.

## Usage

``` r
plot_weightfunction(x, ...)

# S3 method for class 'brma'
plot_weightfunction(
  x,
  rescale_p_values = TRUE,
  prior = FALSE,
  plot_type = "base",
  show_data = TRUE,
  dots_data = NULL,
  dots_prior = NULL,
  ...
)
```

## Arguments

- x:

  a fitted `brma` object with a weightfunction/selection component, such
  as
  [`bselmodel()`](https://fbartos.github.io/RoBMA/reference/bselmodel.md)
  or a [`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md)
  object with weightfunction priors.

- ...:

  list of additional graphical arguments to be passed to the plotting
  function. Supported arguments are `lwd`, `lty`, `col`, `col.fill`,
  `xlab`, `ylab`, `main`, `xlim`, `ylim` to adjust the line thickness,
  line type, line color, fill color, x-label, y-label, title, x-axis
  range, and y-axis range respectively.

- rescale_p_values:

  whether to rescale p-values to the interval shown by the
  weightfunction plot. Defaults to `TRUE`.

- prior:

  whether prior distribution should be added to figure. Defaults to
  `FALSE`.

- plot_type:

  whether to use a base plot `"base"` or ggplot2 `"ggplot"` for
  plotting. Defaults to `"base"`.

- show_data:

  whether observed one-sided p-values should be shown as rug marks on
  the weightfunction axis. Defaults to `TRUE`.

- dots_data:

  list of additional graphical arguments for observed p-value rug marks.
  Supported arguments include `col`/`color`, `alpha`,
  `lwd`/`linewidth`/`size`, `side`/`rug_side`, and
  `height`/`rug_height`/`ticksize`.

- dots_prior:

  list of additional graphical arguments to be passed to the plotting
  function of the prior distribution. Supported arguments are `lwd`,
  `lty`, `col`, and `col.fill`, to adjust the line thickness, line type,
  line color, and fill color of the prior distribution respectively.

## Value

`plot_weightfunction.brma` returns either `NULL` if `plot_type = "base"`
or a `ggplot2` object if `plot_type = "ggplot"`. The method errors for
fitted objects without a weightfunction component.

## See also

[`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")
  fit <- bselmodel(
    yi      = yi,
    vi      = vi,
    data    = dat.lehmann2018,
    measure = "SMD",
    seed    = 1,
    silent  = TRUE
  )

  plot_weightfunction(fit)
  plot_weightfunction(fit, prior = TRUE)
}
} # }

```
