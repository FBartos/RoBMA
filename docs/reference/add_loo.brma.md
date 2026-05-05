# Add LOO-PSIS to brma Objects

Compute approximate leave-one-out cross-validation (LOO-CV) using Pareto
smoothed importance sampling (PSIS) for brma model objects and store the
result in the object.

## Usage

``` r
# S3 method for class 'brma'
add_loo(object, unit = "estimate", r_eff = NULL, parallel = FALSE, ...)
```

## Arguments

- object:

  a brma model object.

- unit:

  output/deletion unit. `"estimate"` computes one contribution per
  effect-size estimate. `"cluster"` computes one contribution per
  cluster and is available only for multilevel models.

- r_eff:

  optional vector of relative effective sample sizes. If not provided,
  it is computed from the log-likelihood values.

- parallel:

  Logical. If `TRUE`,
  [`loo::relative_eff()`](https://mc-stan.org/loo/reference/relative_eff.html)
  and [`loo::loo()`](https://mc-stan.org/loo/reference/loo.html) use
  `RoBMA.get_option("max_cores")`. Log-likelihood construction is
  unchanged. If `FALSE`, those computations use one core.

- ...:

  additional arguments (currently ignored).

## Value

The brma object with the LOO result stored in `object[["loo"]][[unit]]`.

## Details

With `unit = "estimate"`, LOO-CV is computed with one contribution per
effect-size estimate. For binomial and Poisson models, each pair of
counts (ai/ci or x1i/x2i) that defines a single effect size estimate is
treated as one contribution.

With `unit = "cluster"`, LOO-CV is computed with one joint contribution
per cluster. For unweighted normal models without selection this uses
the analytic cluster block covariance. Selection, data-weighted normal,
and GLMM models integrate the held-out cluster effect with Gauss-Hermite
quadrature.

For selection models, the LOO evaluates the weighted likelihood,
conditioning on the posterior omega samples.

The PSIS object is essential for model comparison via
[`loo_compare`](https://mc-stan.org/loo/reference/loo_compare.html) and
is automatically saved in the loo result. RoBMA stores target metadata
so comparisons can reject mismatched data, unit, or conditioning-depth
targets.

**Important for model comparison:** When comparing models via
[`loo_compare`](https://mc-stan.org/loo/reference/loo_compare.html), the
selection is based on expected out-of-sample predictive performance.
This evaluates how well models predict *new* observations, not how well
they fit the observed data.

## References

(Vehtari et al. 2017) (Vehtari et al. 2024)

## See also

[`loo.brma`](https://fbartos.github.io/RoBMA/reference/loo.brma.md),
[`loo`](https://mc-stan.org/loo/reference/loo.html),
[`loo_compare`](https://mc-stan.org/loo/reference/loo_compare.html),
[`pareto_k_ids`](https://mc-stan.org/loo/reference/pareto-k-diagnostic.html)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")
  fit <- bPET(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")

  fit <- add_loo(fit)
  loo_fit <- loo(fit)
  print(loo_fit)
  plot(loo_fit)
}
} # }
```
