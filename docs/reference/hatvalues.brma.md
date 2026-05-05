# Hat Values for brma Objects

Computes hat values (leverages) from a fitted brma object. Returns
posterior mean leverages, one for each observation.

## Usage

``` r
# S3 method for class 'brma'
hatvalues(model, ...)
```

## Arguments

- model:

  a fitted brma object

- ...:

  additional arguments (currently ignored)

## Value

A numeric vector of length `K`, one posterior mean leverage per
observation.

## Details

Hat values (leverages) measure the influence of each observation on the
fitted values. In a Bayesian meta-analysis, the random effects variance
\\\tau^2\\ is uncertain, so the hat matrix depends on the posterior
samples of \\\tau^2\\.

This function computes the diagonal elements of the hat matrix: \$\$h_i
= (X (X^T W X)^{-1} X^T W)\_{ii}\$\$ where \\W\\ is the weight matrix
inverse to the marginal variance matrix.

The hat matrix is computed for each posterior draw and then averaged
over draws, matching the vector output shape used by `metafor`.

This method is available only for normal outcome models without
weightfunction selection.

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE) &&
    requireNamespace("metafor", quietly = TRUE)) {
  data(dat.bcg, package = "metadat")
  dat <- metafor::escalc(
    measure = "RR",
    ai      = tpos,
    bi      = tneg,
    ci      = cpos,
    di      = cneg,
    data    = dat.bcg
  )

  fit <- brma(yi = yi, vi = vi, data = dat, measure = "RR")
  hatvalues(fit)
}
} # }
```
