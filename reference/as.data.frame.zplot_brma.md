# Convert Z-Plot Summaries to a Data Frame

Converts the table displayed for a `zplot_brma` or `summary.zplot_brma`
object to a component-aware long data frame. Printed quantile labels are
returned as syntactic `CI_` column names.

## Usage

``` r
# S3 method for class 'zplot_brma'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  stringsAsFactors = FALSE,
  ...
)

# S3 method for class 'summary.zplot_brma'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  stringsAsFactors = FALSE,
  ...
)
```

## Arguments

- x:

  a `zplot_brma` or `summary.zplot_brma` object.

- row.names:

  `NULL` or a character vector giving the row names.

- optional:

  logical; passed to the final data-frame coercion.

- stringsAsFactors:

  accepted for compatibility with
  [`data.frame()`](https://rdrr.io/r/base/data.frame.html).

- ...:

  additional arguments passed to
  [`summary()`](https://rdrr.io/r/base/summary.html) for a `zplot_brma`
  object and otherwise unused.

## Value

A plain `data.frame` with leading `component` and `parameter` columns.
