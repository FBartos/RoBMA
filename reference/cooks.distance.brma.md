# Cook's Distance for brma Objects

Computes Cook's distance for a fitted brma object. Cook's distance
measures the aggregate influence of each observation on the model
coefficients.

## Usage

``` r
# S3 method for class 'brma'
cooks.distance(model, ...)
```

## Arguments

- model:

  a fitted normal-outcome `brma` object without a weightfunction
  component.

- ...:

  additional arguments (currently ignored).

## Value

A numeric vector of Cook's distance values, one for each observation.

## Details

Cook's distance is computed as a PSIS leave-one-out deletion diagnostic.
For each observation \\i\\, normalized PSIS weights estimate the fitted
values under the leave-one-out posterior. The distance is the squared
posterior Mahalanobis distance between the full-data and leave-one-out
fitted-value vectors: \$\$D_i = \Delta_i' V\_\mu^+ \Delta_i\$\$

where \\\Delta_i = \hat{\mu} - \hat{\mu}\_{(-i)}\\, \\V\_\mu^+\\ is the
generalized inverse of the full-posterior fitted-value covariance. This
uses the unscaled meta-analytic convention of metafor and Viechtbauer
and Cheung (2010) , without division by the fixed-effect model rank.
Values can differ from metafor because RoBMA uses posterior moments and
PSIS deletion rather than classical refits. For
[`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
known-`V` models, Cook's distance uses estimate-unit PSIS weights. With
correlated known-`V`, deletion is conditional estimate deletion and the
reported fitted-value target is the fixed-location mean \\\mu =
X\beta\\; sampled or marginalized random effects are not included in the
reported fitted value. Selection
[`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
influence is currently unavailable. The PSIS fixed-location moments do
not require the inverse marginal covariance used by
[`hatvalues()`](https://rdrr.io/r/stats/influence.measures.html).

## References

(Viechtbauer and Cheung 2010)

## See also

[`influence.brma`](https://fbartos.github.io/RoBMA/reference/influence.brma.md),
[`dffits.brma`](https://fbartos.github.io/RoBMA/reference/dffits.brma.md),
[`hatvalues.brma`](https://fbartos.github.io/RoBMA/reference/hatvalues.brma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")
  fit <- bPET(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
  fit <- add_loo(fit)

  cooks.distance(fit)
}
} # }
```
