# Convert Influence Results to a Data Frame

Converts the influence and DFBETAS tables displayed by an `infl.brma`
object to one component-aware long data frame.

## Usage

``` r
# S3 method for class 'infl.brma'
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

  an `infl.brma` object.

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
