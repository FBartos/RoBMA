# Convert brma_samples Results to Data Frames

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) and
[`data.frame()`](https://rdrr.io/r/base/data.frame.html) of a
`brma_samples` result return its posterior draws (one row per draw, one
column per parameter), as for the underlying matrix. The posterior
summary table displayed by
[`print()`](https://rdrr.io/r/base/print.html) is
`as.data.frame(summary(x))`: a long data frame with component and
parameter identifiers (multi-component paths use `/` as a separator),
which is also the result for a multi-component `brma_samples_list`.
Printed quantile labels such as `0.025` and `PI 0.025` are returned as
the syntactic column names `CI_0.025` and `PI_0.025`, respectively.

## Usage

``` r
# S3 method for class 'brma_samples'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'summary.brma_samples'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  stringsAsFactors = FALSE,
  ...
)

# S3 method for class 'brma_samples_list'
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

  a `brma_samples`, `summary.brma_samples`, or multi-component
  `brma_samples_list` result

- row.names:

  `NULL` or a character vector giving the row names for the resulting
  data frame. Custom row names are unsupported when `format = "list"`.

- optional:

  logical; passed to the final data-frame coercion method.

- ...:

  for a `brma_samples` result, additional arguments passed to the matrix
  coercion; for a `brma_samples_list`, additional arguments passed to
  [`summary.brma_samples()`](https://fbartos.github.io/RoBMA/dev/reference/summary.brma_samples.md)

- stringsAsFactors:

  logical; passed to the output data frame.

- format:

  for multi-component results, whether to return a single `"long"` data
  frame or a named, possibly nested `"list"` with a plain data frame at
  each leaf. Defaults to `"long"`.

## Value

For an individual `brma_samples` object, a `data.frame` of the posterior
draws. For a `summary.brma_samples` object, a plain long `data.frame`
containing the displayed posterior summary statistics and leading
`component` and `parameter` columns. For a `brma_samples_list`, a long
data frame by default or, with `format = "list"`, a named, possibly
nested list of such data frames.
