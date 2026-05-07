# Residuals for brma Objects

Computes residuals (observed minus fitted values) from a fitted brma
object, with options for standardization.

## Usage

``` r
# S3 method for class 'brma'
residuals(
  object,
  type = "outcome",
  unit = "estimate",
  conditioning_depth = "marginal",
  bias_adjusted = FALSE,
  ...
)
```

## Arguments

- object:

  a fitted brma object

- type:

  the type of residuals to compute. Options are:

  - `"outcome"` (default): Raw residuals (observed - fitted). Available
    for all model types.

  - `"pearson"`: Pearson (semi-standardized) residuals, dividing raw
    residuals by the marginal standard error \\\sqrt{v_i + \tau^2}\\.
    Only available for normal outcome models without selection
    (weightfunction) bias adjustment.

  - `"rstandard"`: Internally standardized residuals, dividing raw
    residuals by their standard errors computed using the hat matrix.
    Only available for normal outcome models without selection
    (weightfunction) bias adjustment.

  - `"LOO-PIT"` (alias: `"rstudent"`): Leave-one-out probability
    integral transform (PIT) residuals computed via Pareto smoothed
    importance sampling. The LOO-CDF value for each observation is
    computed and transformed to standard normal quantiles via
    \\\Phi^{-1}(u_i)\\. Under a correctly specified model, these
    residuals should follow a standard normal distribution. This is the
    recommended standardized residual for Bayesian models as it properly
    accounts for estimation uncertainty and leverage. Available for all
    model types. Note: This requires that the loo has been computed
    previously (see
    [`add_loo()`](https://fbartos.github.io/RoBMA/reference/add_loo.brma.md)
    function).

- unit:

  output unit. Only `"estimate"` is implemented currently..

- conditioning_depth:

  conditioning depth for non-LOO residuals. `"marginal"` uses fixed
  effects only, `"cluster"` conditions on cluster-level random effects,
  and `"estimate"` conditions on the full estimate-level fitted value.
  LOO-PIT residuals always use the estimate-unit LOO target.

- bias_adjusted:

  whether residuals should be computed from bias-adjusted fitted values.
  Defaults to `FALSE`, which means residuals are computed as the
  difference between observed values and raw (biased) predictions
  including PET/PEESE terms. Set to `TRUE` to compute residuals from
  bias-corrected fitted values. Applies to outcome and Pearson
  residuals. Note that bias-adjusted residuals are not residuals in the
  traditional sense.

- ...:

  additional arguments.

## Value

A numeric vector of residual means, one per estimate.

## Details

Raw residuals (`type = "outcome"`) are computed as: \$\$e_i = y_i -
\hat{\mu}\_i\$\$ where \\y_i\\ is the observed effect size and
\\\hat{\mu}\_i\\ is the fitted value (prediction from the fixed
effects).

Pearson residuals (`type = "pearson"`) divide raw residuals by the
marginal standard error: \$\$r_i^{Pearson} = \frac{e_i}{\sqrt{v_i +
\tau^2}}\$\$ where \\v_i\\ is the sampling variance and \\\tau^2\\ is
the relevant heterogeneity component. Only available for normal outcome
models without selection (weightfunction) bias adjustment.

Standardized residuals (`type = "rstandard"`) use the hat matrix to
compute residual standard errors that account for the uncertainty in
estimated coefficients: \$\$z_i =
\frac{e_i}{\sqrt{\[(I-H)M(I-H)'\]\_{ii}}}\$\$ where \\H\\ is the hat
matrix and \\M\\ is the marginal variance-covariance matrix. For models
without moderators, this simplifies to the Pearson formula. Only
available for normal outcome models without selection (weightfunction)
bias adjustment.

LOO-PIT residuals (`type = "LOO-PIT"`) are the Bayesian equivalent of
studentized deleted residuals (Vehtari et al. 2017) . They are computed
via leave-one-out probability integral transformation: \$\$r_i =
\Phi^{-1}(u_i)\$\$ where \\u_i = \sum_s w\_{is} F(y_i \| \theta^{(s)})\\
is the LOO-weighted CDF value, \\w\_{is}\\ are the normalized PSIS
weights, and \\F\\ is the cumulative distribution function of the
estimate-unit predictive distribution used by LOO. Under a correctly
specified model, LOO-PIT residuals should follow a standard normal
distribution. Unlike traditional standardized residuals, LOO-PIT
residuals properly account for estimation uncertainty and leverage
without requiring a hat matrix. This is the recommended method for
standardized residuals in Bayesian meta-analysis.

For meta-regression models, fitted values incorporate moderator effects.
For models without moderators, all fitted values equal the pooled
effect.

For GLMM models (binomial or Poisson), observed effect sizes and their
sampling variances are computed from the raw frequency data using the
same formulas as
[`metafor::escalc`](https://wviechtb.github.io/metafor/reference/escalc.html)
with the default zero-cell adjustment (adding 0.5 to all cells when any
cell is zero). GLMM residuals and LOO-PIT values are therefore
approximate effect-size-scale diagnostics, not exact PIT diagnostics for
the raw count likelihood.

The residuals are computed separately for each posterior sample,
naturally propagating uncertainty in model parameters to the residuals.

## References

Vehtari A, Gelman A, Gabry J (2017). “Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC.” *Statistics
and computing*, **27**(5), 1413–1432.
[doi:10.1007/s11222-016-9696-4](https://doi.org/10.1007/s11222-016-9696-4)
.

## See also

[`predict.brma()`](https://fbartos.github.io/RoBMA/reference/predict.brma.md),
[`blup.brma()`](https://fbartos.github.io/RoBMA/reference/blup.brma.md),
[`pooled_effect()`](https://fbartos.github.io/RoBMA/reference/pooled_effect.md),
[`rstandard.brma()`](https://fbartos.github.io/RoBMA/reference/rstandard.brma.md),
[`rstudent.brma()`](https://fbartos.github.io/RoBMA/reference/rstudent.brma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")

  fit <- bPET(
    yi      = yi,
    vi      = vi,
    data    = dat.lehmann2018,
    measure = "SMD",
    seed    = 1,
    silent  = TRUE
  )

  # raw residuals (default)
  residuals(fit)

  # Pearson and internally standardized residuals
  residuals(fit, type = "pearson")
  residuals(fit, type = "rstandard")

  # LOO-PIT residuals require stored LOO
  fit <- add_loo(fit)
  residuals(fit, type = "LOO-PIT")

  # residuals from bias-adjusted predictions
  residuals(fit, bias_adjusted = TRUE)

  # check Pareto k diagnostics
  plot(loo(fit))
}
} # }
```
