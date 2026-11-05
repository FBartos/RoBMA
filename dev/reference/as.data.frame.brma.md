# Convert brma Summaries to a Data Frame

Converts the non-empty tables displayed by a fitted `brma` object or its
`summary.brma` result to one plain, long-form data frame. The leading
`component` column identifies the summary section and `parameter`
retains the displayed row label. Printed quantile labels are returned as
syntactic `CI_` column names.

## Usage

``` r
# S3 method for class 'brma'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  stringsAsFactors = FALSE,
  ...
)

# S3 method for class 'summary.brma'
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

  a fitted `brma` or `summary.brma` object.

- row.names:

  `NULL` or a character vector giving the row names.

- optional:

  logical; passed to the final data-frame coercion.

- stringsAsFactors:

  accepted for compatibility with
  [`data.frame()`](https://rdrr.io/r/base/data.frame.html).

- ...:

  for a fitted object, additional arguments passed to
  [`summary.brma()`](https://fbartos.github.io/RoBMA/dev/reference/summary.brma.md);
  otherwise unused.

## Value

A plain `data.frame` containing all displayed summary rows. An object
created with `only_data = TRUE` or `only_priors = TRUE` has no fitted
model to export and stops with an error of class `RoBMA_unfitted_object`
(also of class `RoBMA_unsupported_request`);
[`summary()`](https://rdrr.io/r/base/summary.html) shows the data of an
`only_data` object, and the summary of an `only_priors` object shows the
priors in its `priors` element.
