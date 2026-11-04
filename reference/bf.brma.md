# Bayes Factor for brma Objects

Compute the Bayes factor comparing two brma models.

## Usage

``` r
# S3 method for class 'brma'
bf(x1, x2, log = FALSE, ...)

# S3 method for class 'brma'
bayes_factor(x1, x2, log = FALSE, ...)
```

## Arguments

- x1:

  a brma model object (numerator).

- x2:

  a brma model object (denominator).

- log:

  logical; if `TRUE`, the log Bayes factor is returned. Default is
  `FALSE`.

- ...:

  additional arguments (currently not used).

## Value

For single estimates, a list of class `"bf_default"` with components:

- `bf`: (scalar) value of the Bayes factor in favor of `x1` over `x2`.

- `log`: Boolean indicating whether `bf` corresponds to the log Bayes
  factor.

For repeated estimates, a list of class `c("bf.brma", "bf_bridge_list")`
with a vector `bf`, scalar `bf_median_based`, and `log`. It uses the
upstream repeated-Bayes-factor print method. Both
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) and
[`data.frame()`](https://rdrr.io/r/base/data.frame.html) return a long
data frame with `component`, `parameter`, and full-precision numeric
`value` columns, containing every repetition and the displayed
median-based, range, and interquartile-range summaries.

## Details

Computes the Bayes factor in favor of the model `x1` over the model
`x2`. The marginal likelihoods must first be computed using
[`add_marglik`](https://fbartos.github.io/RoBMA/reference/add_marglik.brma.md).
Both models must be fitted to the same outcome target/data, including
outcome type and, when present, likelihood weights. Cluster identifiers
define model structure and are not part of the outcome-data identity. A
known sampling covariance without correlated blocks (e.g.,
`V = diag(sei^2)` in
[`brma.mv`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)) is
the outcome data of its standard errors, whereas a correlated `V` is a
different target.

When either model contains repeated bridge estimates, the shorter
sequence is recycled with a warning, following bridgesampling.
Repetitions are compared in order. The median-based estimate is computed
from the separate medians of the log marginal likelihoods before
recycling, rather than the median of the repetition-wise Bayes factors.
Standardized marginal-likelihood results must be finite. Legacy raw
`"bridge_list"` results retain failed repetitions in the returned
estimates, with upstream missing-value behavior when printed.

## See also

[`add_marglik`](https://fbartos.github.io/RoBMA/reference/add_marglik.brma.md),
[`bridge_sampler.brma`](https://fbartos.github.io/RoBMA/reference/bridge_sampler.brma.md),
[`post_prob.brma`](https://fbartos.github.io/RoBMA/reference/post_prob.brma.md),
[`logml.brma`](https://fbartos.github.io/RoBMA/reference/logml.brma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")
  fit1 <- brma(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
  fit2 <- brma(
    yi           = yi,
    vi           = vi,
    data         = dat.lehmann2018,
    measure      = "SMD",
    prior_effect = FALSE
  )

  fit1 <- add_marglik(fit1)
  fit2 <- add_marglik(fit2)

  bf(fit1, fit2)
}
} # }
```
