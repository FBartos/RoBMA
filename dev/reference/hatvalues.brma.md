# Hat Values for brma Objects

Computes hat values (leverages) from a fitted brma object. Returns
posterior mean leverages, one for each observation.

## Usage

``` r
# S3 method for class 'brma'
hatvalues(model, max_samples = Inf, ...)
```

## Arguments

- model:

  a fitted brma object

- max_samples:

  maximum posterior draws used for known-`V`
  [`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
  marginal covariance computations. Defaults to `Inf`. Finite values
  deterministically thin draws across posterior row order.

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
= (X (X^T W X)^{-1} X^T W)\_{ii}\$\$ where \\W\\ is the fitted precision
matrix. Without additional likelihood weights it is the inverse marginal
covariance; likelihood weights modify estimation precision without
redefining the original outcome covariance used to standardize
residuals.

The hat matrix is computed for each posterior draw and then averaged
over draws, matching the vector output shape used by `metafor`. The
members of a
[`RoBMA()`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.md)
ensemble can differ in their PET and PEESE terms; each draw then uses
the design of the member it belongs to (the location columns plus the
`PET` column of PET members or the `PEESE` column of PEESE members), and
the reported leverage is the posterior mean over all draws. For
known-`V`
[`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
random-formula models, leverages use the marginal GLS covariance \\V +
ZGZ'\\.

This method is available only for normal outcome models without
weightfunction selection.

For known-`V` models, the leverages require applying the inverse of the
marginal covariance. A covariance that is too ill-conditioned for that
(a reciprocal condition number of about `2e-10` or less in any posterior
draw used) is unavailable: the call stops with an error of class
`RoBMA_known_v_precision_unavailable`, whose parent class is
`RoBMA_numerical_unavailable`, instead of returning a rounded number.
The condition carries the draw positions (`rows`) and their
`reciprocal_condition`.

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
