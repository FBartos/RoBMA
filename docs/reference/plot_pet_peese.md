# Plot PET-PEESE Fit of brma Object

`plot_pet_peese` visualizes posterior (and prior) PET-PEESE fit of a
brma object.

## Usage

``` r
plot_pet_peese(x, ...)

# S3 method for class 'brma'
plot_pet_peese(
  x,
  show_data = TRUE,
  prior = FALSE,
  plot_type = "base",
  dots_prior = NULL,
  ...
)
```

## Arguments

- x:

  a fitted brma object with a PET or PEESE component.

- ...:

  list of additional graphical arguments to be passed to the plotting
  function. Supported arguments are `lwd`, `lty`, `col`, `col.fill`,
  `xlab`, `ylab`, `main`, `xlim`, `ylim`, `pch`, `bg`, `cex`, and `size`
  to adjust the line thickness, line type, line color, fill color,
  x-label, y-label, title, x-axis range, y-axis range, and observed data
  point style respectively.

- show_data:

  whether the observed effect sizes should be added to figure. Defaults
  to `TRUE`.

- prior:

  whether prior distribution should be added to figure. Defaults to
  `FALSE`.

- plot_type:

  whether to use a base plot `"base"` or ggplot2 `"ggplot"` for
  plotting. Defaults to `"base"`.

- dots_prior:

  list of additional graphical arguments to be passed to the plotting
  function of the prior distribution. Supported arguments are `lwd`,
  `lty`, `col`, and `col.fill`, to adjust the line thickness, line type,
  line color, and fill color of the prior distribution respectively.

## Value

`plot_pet_peese` returns either `NULL` invisibly if `plot_type = "base"`
or a ggplot2 object if `plot_type = "ggplot"`.

## Details

The plot shows observed `yi` values against `sei`. PET regression uses
\\\mu + PET \cdot se_i\\; PEESE regression uses \\\mu + PEESE \cdot
se_i^2\\, with the fitted effect direction.

## See also

[`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")

  fit <- bPET(
    yi      = yi,
    vi      = vi,
    data    = dat.lehmann2018,
    measure = "SMD",
    seed    = 1,
    silent  = TRUE
  )

  plot_pet_peese(fit)
}
} # }

```
