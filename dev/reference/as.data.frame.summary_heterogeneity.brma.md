# Convert Heterogeneity Summaries to Data Frames

Converts heterogeneity summaries to plain long data frames with leading
`component` and `parameter` columns. A multi-component result can
instead be returned as a named list with `format = "list"`.

## Usage

``` r
# S3 method for class 'summary_heterogeneity.brma'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  stringsAsFactors = FALSE,
  ...
)

# S3 method for class 'summary_heterogeneity.brma_list'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  format = c("long", "list"),
  stringsAsFactors = FALSE,
  ...
)
```

## Arguments

- x:

  a `summary_heterogeneity.brma` or `summary_heterogeneity.brma_list`
  object.

- row.names:

  `NULL` or a character vector giving the row names for the result.
  Custom row names are unsupported when `format = "list"`.

- optional:

  logical; passed to the final data-frame coercion.

- stringsAsFactors:

  accepted for compatibility with
  [`data.frame()`](https://rdrr.io/r/base/data.frame.html).

- ...:

  unused additional arguments.

- format:

  for multi-component results, whether to return one `"long"` data frame
  or a named `"list"` of data frames.

## Value

A plain `data.frame`, or a named list of data frames when
`format = "list"`.
