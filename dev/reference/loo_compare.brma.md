# Compare brma Models Using LOO

Compare multiple brma models using LOO-PSIS cross-validation. This is a
convenience wrapper around
[`loo_compare`](https://mc-stan.org/loo/reference/loo_compare.html).

## Usage

``` r
# S3 method for class 'brma'
loo_compare(x, ..., unit = "estimate")
```

## Arguments

- x:

  a brma model object (the first model to compare).

- ...:

  additional brma model objects or `loo` objects to compare.

- unit:

  output/deletion unit used when extracting LOO from brma objects.

## Value

The comparison returned by
[`loo_compare`](https://mc-stan.org/loo/reference/loo_compare.html) for
the stored results: a data frame of class `"compare.loo"` with one row
per model, ordered by expected predictive performance, which prints as
loo prints it. The models are named `model1`, `model2`, ... unless the
call names them.

## Details

This function compares models based on their expected out-of-sample
predictive performance (ELPD).

**Important for model comparison:** When comparing models via
[`loo_compare`](https://mc-stan.org/loo/reference/loo_compare.html), the
selection is based on expected predictive performance for the documented
deletion target, not in-sample fit. For dependent estimates, deleting
one estimate may retain other observations from its fitted dependency
block; this is not new-group prediction. RoBMA rejects comparisons with
different outcome targets/data, `unit`, or retained context. A known
sampling covariance without correlated blocks (e.g., `V = diag(sei^2)`
in
[`brma.mv`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md))
is the outcome data of its standard errors, whereas a correlated `V` is
a different target. The same compatibility key is applied separately to
LOO and WAIC; the two criteria cannot be mixed in one comparison table.
Sampled versus marginalized local-effect provenance does not by itself
block comparison.

## See also

[`loo.brma`](https://fbartos.github.io/RoBMA/dev/reference/loo.brma.md),
[`loo_compare`](https://mc-stan.org/loo/reference/loo_compare.html)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")

  fit_bias <- RoBMA(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
  fit_nobias <- BMA(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")

  fit_bias <- add_loo(fit_bias)
  fit_nobias <- add_loo(fit_nobias)

  loo_compare(fit_bias, fit_nobias)
  loo_compare(loo(fit_bias), loo(fit_nobias))
}
} # }
```
