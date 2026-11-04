# Summary of Heterogeneity for brma Objects

Computes the absolute heterogeneity (tau, tau^2) and relative measures
of heterogeneity (I^2, H^2) for a fitted brma object.

## Usage

``` r
# S3 method for class 'brma'
summary_heterogeneity(object, probs = c(0.025, 0.975), component = "all", ...)

# S3 method for class 'summary_heterogeneity.brma_list'
print(x, ...)
```

## Arguments

- object:

  a fitted brma object

- probs:

  quantiles of the posterior distribution to be displayed. Defaults to
  `c(.025, .975)` for 95% credible intervals.

- component:

  heterogeneity component to return for
  [`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
  models. Defaults to `"all"`. Use `"total"` for the variance-additive
  total heterogeneity. Random-formula allocation nodes can be selected
  by their displayed name, such as `"study/esid"`.

- ...:

  additional arguments (currently ignored)

- x:

  a `summary_heterogeneity.brma_list` object

## Value

A list of class `summary_heterogeneity.brma` containing:

- `estimates`: A `BayesTools_table` with heterogeneity statistics

- `component`: The public heterogeneity-component name

For decomposed
[`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
models with multiple selected components, a named list of such summary
objects is returned.

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
follows the approach described in the metafor documentation. For
[`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md)
ensembles whose members differ in their PET and PEESE terms, each draw
uses the typical sampling variance of its member's design (see
[`hatvalues.brma`](https://fbartos.github.io/RoBMA/reference/hatvalues.brma.md)).

For [`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
models, the method reports `tau` and `tau2` for each selected
random-effect component. A genuine variance-additive aggregate is
reported as `tau_total` and `tau2_total`; a mean-variance allocation
scale is reported as `tau_common` and `tau2_common`. These RoBMA I/O
names map to BayesTools' general random-effect names without changing
the fitted quantities. Formula-random summaries are present only when a
`random` formula is declared. A random-effect quantity that the fitted
block cannot represent on the original predictor scale (for example the
intercept SD, variance and correlations of a block with factor-specific
slopes of a standardized predictor) is `NA`, with a row note, in the row
that shows it. This table shows the SDs, variances and correlations, so
all of them can be `NA` here;
[`summary.brma`](https://fbartos.github.io/RoBMA/reference/summary.brma.md)
shows the correlations and the fitted-scale multipliers, so only its
correlation rows can be `NA`. Public random-effect correlations (`rho`)
are included in their component's table, including homogeneous
structures such as [`ar()`](https://rdrr.io/r/stats/ar.html), `cs()`,
and `car()`. A variance-additive total across multiple components has no
single correlation parameter and reports only its total SD and variance.
When a known group covariance \\R\\ is supplied, `tau` is the fitted
base multiplier of that covariance kernel. It need not equal every row's
marginal standard deviation when \\\mathrm{diag}(R)\\ is not one.
Row-specific marginal standard deviations remain available from
`predict(type = "terms.scale")` and in covariance-based prediction and
diagnostics. \\I^2\\ and \\H^2\\ use the original sampling variances and
the fitted location design to define the typical sampling variance.
Likelihood weights affect these summaries through the posterior
heterogeneity but do not replace the sampling variances \\v_i\\ by \\v_i
/ w_i\\, so rescaling all weights does not change the reference variance
(as in metafor). Relative \\I^2\\ and \\H^2\\ summaries are not reported
for general known-V covariance structures, nor for a meta-regression
whose location design leaves no residual degrees of freedom (as many
estimable coefficients as estimates), where the typical sampling
variance is undefined. For `component = "total"`, independent component
variances are summed before reporting `tau2_total`; `tau_total` is the
square root of this variance draw. When a random-formula model uses a
shared total-SD plus variance-allocation node, `component = "all"` also
includes an allocation-node table with the appropriate aggregate SD and
variance plus `tau2_prop(<block>)` rows. For gated allocations, the
aggregate is the realized model-averaged total and includes the all-off
zero branch; proportions are realized shares conditional on positive
total variance. For nested formulas such as
`random = ~ 1 | study / esid`, this table is displayed under the
user-facing path `"study/esid"`. For heterogeneous structured random
effects, such as `random = ~ har(time | study)`, the allocation table
reports the exhaustive semantic allocation family: aggregate SD and
variance, level-specific component SDs and variances, correlations, SD
multipliers, and variance multipliers. For row-shaped scale sources,
absolute level SDs and variances use the observed rows assigned to each
level; their aggregate uses the observed-source RMS. A redundant block
owner is omitted for a bare formula or unnamed one-entry list and
retained for explicitly named one-entry lists and multiple components.

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
