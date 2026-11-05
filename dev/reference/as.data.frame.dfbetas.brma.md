# Convert DFBETAS Results to a Data Frame

Converts a `dfbetas.brma` object to a component-aware long data frame
while retaining its displayed coefficient columns.

## Usage

``` r
# S3 method for class 'dfbetas.brma'
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

  a `dfbetas.brma` object.

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
