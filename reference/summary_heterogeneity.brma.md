# Summary of Heterogeneity for brma Objects

Computes the absolute heterogeneity (tau, tau^2) and relative measures
of heterogeneity (I^2, H^2) for a fitted brma object.

## Usage

``` r
# S3 method for class 'brma'
summary_heterogeneity(object, probs = c(0.025, 0.975), ...)
```

## Arguments

- object:

  a fitted brma object

- probs:

  quantiles of the posterior distribution to be displayed. Defaults to
  `c(.025, .975)` for 95% credible intervals.

- ...:

  additional arguments (currently ignored)

## Value

A list of class `summary_heterogeneity.brma` containing:

- `estimates`: A `BayesTools_table` with heterogeneity statistics

## Details

For standard (2-level) random-effects models, the function reports:

- `tau`: Between-study standard deviation

- `tau2`: Between-study variance

- `I2`: Percentage of total variance due to heterogeneity

- `H2`: Ratio of total to sampling variance

For multilevel (3-level) models with nested effects, the function
additionally partitions heterogeneity into estimate-level and
cluster-level components:

- `rho`: Proportion of heterogeneity variance allocated to clusters

- `tau [within]`: Estimate-level standard deviation

- `tau [between]`: Cluster-level standard deviation

- `tau2 [within]`: Estimate-level variance

- `tau2 [between]`: Cluster-level variance

- `I2 [within]`: Percentage of variance due to estimate-level
  heterogeneity

- `I2 [between]`: Percentage of variance due to cluster-level
  heterogeneity

For location-scale models, `tau2` aggregates the observation-specific
heterogeneity variances \\\tau_i^2\\; the corresponding `tau` summary is
the square root of this aggregate variance. The relative \\I^2\\ and
\\H^2\\ measures average the observation-specific indices.

The I^2 and H^2 statistics are computed following the metafor package
implementation, using the "typical" sampling variance formula from
Higgins and Thompson (2002) . For multilevel models, the partitioned I^2
follows the approach described in the metafor documentation.

## References

Higgins JP, Thompson SG (2002). “Quantifying heterogeneity in a
meta-analysis.” *Statistics in Medicine*, **21**(11), 1539–1558.
[doi:10.1002/sim.1186](https://doi.org/10.1002/sim.1186) .

## See also

[`pooled_heterogeneity()`](https://fbartos.github.io/RoBMA/reference/pooled_heterogeneity.md),
[`summary.brma()`](https://fbartos.github.io/RoBMA/reference/summary.brma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")
  fit <- brma(
    yi      = yi,
    vi      = vi,
    data    = dat.lehmann2018,
    measure = "SMD",
    seed    = 1,
    silent  = TRUE
  )

  summary_heterogeneity(fit)
}
} # }
```
