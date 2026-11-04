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
  max_samples = Inf,
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
    residuals by the SD of the corresponding outcome target: marginal
    includes all applicable heterogeneity, cluster depth includes
    within- cluster heterogeneity, and estimate depth uses sampling
    variance only. Only available for normal outcome models without
    selection (weightfunction) bias adjustment.

  - `"rstandard"`: Internally standardized projection residuals, divided
    by their corresponding SDs computed using the hat matrix. Only
    available for normal outcome models without selection
    (weightfunction) bias adjustment.

  - `"LOO-PIT"` (alias: `"rstudent"`): Leave-one-out probability
    integral transform (PIT) residuals computed via Pareto smoothed
    importance sampling. The LOO-CDF value for each observation is
    computed and transformed to standard normal quantiles via
    \\\Phi^{-1}(u_i)\\. Under a correctly specified model, these
    residuals should follow a standard normal distribution. This is the
    recommended standardized residual for Bayesian models as it properly
    accounts for estimation uncertainty and leverage. Available for
    normal outcome models and, as a descriptive normal approximation on
    the effect-size scale, for binomial and Poisson GLMMs (see Details).
    Note: This requires that the loo has been computed previously (see
    [`add_loo()`](https://fbartos.github.io/RoBMA/reference/add_loo.brma.md)
    function).

- unit:

  output unit. Only `"estimate"` is implemented currently.

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

- max_samples:

  maximum posterior draws used for known-`V`
  [`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
  marginal covariance computations in Pearson and internally
  standardized residuals. Defaults to `Inf`. Finite values
  deterministically thin draws across posterior row order.

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
matrix and \\M\\ is the marginal variance-covariance matrix. Without
additional likelihood weights, the GLS projection satisfies \\HM = MH' =
HMH'\\, so the denominator equals \\M\_{ii} -
\[X(X'M^{-1}X)^{-1}X'\]\_{ii}\\, which is strictly smaller than the
Pearson denominator \\\sqrt{M\_{ii}}\\. This holds with no moderators
too, where fitting the pooled effect alone removes \\1/(1'M^{-1}1)\\
from every residual variance, so standardized and Pearson residuals do
not coincide. Only available for normal outcome models without selection
(weightfunction) bias adjustment.

Likelihood `weights` affect fitting, the fitted projection, and PSIS
deletion weights; they do not divide the observation variance used for
diagnostic scaling. Pearson residuals and the predictive CDF/moments in
[`rstudent()`](https://rdrr.io/r/stats/influence.measures.html) use the
original outcome law. Internally standardized residuals combine the
fitted weighted projection \\H\\ with that law's covariance \\M\\, as in
metafor diagnostics with custom fitting weights. RoBMA retains its
Bayesian posterior averaging and LOO-PIT transformation; these are not
identical to frequentist deleted-residual z-statistics. Each retained
row produces one residual, without multiplying or repeating it by its
observation weight. Estimate-unit LOO removes the row's entire weighted
likelihood contribution, for integer and non-integer weights alike.

LOO-PIT residuals (`type = "LOO-PIT"`) are the Bayesian equivalent of
studentized deleted residuals (Vehtari et al. 2017) . They are computed
via leave-one-out probability integral transformation: \$\$r_i =
\Phi^{-1}(u_i)\$\$ where \\u_i = \sum_s w\_{is} F(y_i \| \theta^{(s)})\\
is the LOO-weighted CDF value, \\w\_{is}\\ are the normalized PSIS
weights, and \\F\\ is the cumulative distribution function of the
estimate-unit outcome distribution under the leave-one-out fit. With
likelihood weights, PSIS removes the weighted likelihood contribution
while this CDF remains unpowered. For weighted Gaussian cluster models,
the cluster effect is integrated given the retained weighted rows. These
retained moments are paired with the integrated powered-likelihood PSIS
ratios; full-data sampled cluster effects are not reused. Numerically
unresolved retained moments raise `RoBMA_weighted_gaussian_unavailable`,
a `RoBMA_numerical_unavailable` condition. Lower and upper probabilities
are accumulated separately on the log scale, and the smaller tail is
passed directly to the normal quantile function. Therefore extreme
residuals are not truncated by probability clipping. Under a correctly
specified model, LOO-PIT residuals should follow a standard normal
distribution. Unlike traditional standardized residuals, LOO-PIT
residuals properly account for estimation uncertainty and leverage
without requiring a hat matrix. This is the recommended method for
standardized residuals in Bayesian meta-analysis.

For meta-regression models, fitted values incorporate moderator effects.
For models without moderators, all fitted values equal the pooled
effect.

For binomial and Poisson GLMMs the fitted predictive distributions are
discrete, so the LOO-PIT residuals are a descriptive normal
approximation on the effect-size scale: \\F\\ is the normal distribution
of the continuity-corrected log odds ratio or log incidence rate ratio
with its approximate sampling standard error, as metafor computes
residuals of GLMMs. They are not probability integral transforms of the
fitted count likelihood, and the functions say so in a message (plots
label the approximation).

For correlated known-`V`
[`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
models, LOO-PIT residuals use the same conditional estimate-unit target
as LOO, \\p(y_i \mid y\_{-i}, \theta)\\. Random-formula
[`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
models support marginal Pearson and internally standardized residuals
through the marginal covariance \\V + ZGZ'\\ and estimate-depth
residuals through the fitted existing-level random-effect target.
Selection models retain their declared context and original publication
event during deletion. Retained product weights cancel; best selection
also retains the best p-value among the observed outcomes in that event.
Pearson and hat-matrix standardized residuals remain unavailable for
selection models; use LOO-PIT residuals for a likelihood-aware
standardized diagnostic.

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
