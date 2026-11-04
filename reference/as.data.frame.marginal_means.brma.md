# Convert Marginal Means to a Data Frame

Converts the summary displayed for a `marginal_means.brma` object to a
plain long data frame with leading `component` and `parameter` columns.
Printed quantile labels are returned as syntactic `CI_` column names.

## Usage

``` r
# S3 method for class 'marginal_means.brma'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  stringsAsFactors = FALSE,
  ...
)

# S3 method for class 'summary.marginal_means.brma'
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

  a `marginal_means.brma` or `summary.marginal_means.brma` object.

- row.names:

  `NULL` or a character vector giving the row names.

- optional:

  logical; passed to the final data-frame coercion.

- stringsAsFactors:

  accepted for compatibility with
  [`data.frame()`](https://rdrr.io/r/base/data.frame.html).

- ...:

  additional arguments passed to
  [`summary()`](https://rdrr.io/r/base/summary.html).

## Value

A plain long `data.frame` containing the displayed marginal-means
summary and leading `component` and `parameter` columns.
