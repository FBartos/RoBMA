# Check LOO Diagnostics for brma Objects

Check Pareto k diagnostics for a brma model object and warn if any
values are high.

## Usage

``` r
# S3 method for class 'brma'
check_loo(object, unit = "estimate", ...)
```

## Arguments

- object:

  a brma model object.

- unit:

  output/deletion unit. See
  [`add_loo`](https://fbartos.github.io/RoBMA/reference/add_loo.brma.md).

- ...:

  currently unused.

## Value

NULL (throws warning if diagnostics are unreliable).

## Details

LOO must first be computed with
`object <- add_loo(object, unit = unit)`. The method warns when any
Pareto \\k\\ diagnostic is greater than 0.7.
