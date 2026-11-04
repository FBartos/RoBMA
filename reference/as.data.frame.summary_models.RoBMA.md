# Convert Model-Weight Summaries to a Data Frame

Converts every table displayed by a `summary_models.RoBMA` object to one
component-aware long data frame.

## Usage

``` r
# S3 method for class 'summary_models.RoBMA'
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

  a `summary_models.RoBMA` object.

- row.names:

  `NULL` or a character vector giving the row names.

- optional:

  logical; passed to the final data-frame coercion.

- stringsAsFactors:

  accepted for compatibility with
  [`data.frame()`](https://rdrr.io/r/base/data.frame.html).

- ...:

  unused additional arguments.

## Value

A plain `data.frame` with leading `component` and `parameter` columns.
