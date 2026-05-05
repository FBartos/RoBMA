# DFBETAS for brma Objects

Computes DFBETAS (Difference in BETAS, standardized) for a fitted brma
object. DFBETAS measures the influence of each observation on the
estimated model coefficients. Positive values indicate that deleting the
observation yields a smaller estimate, negative values indicate that
deleting the observation yields a larger estimate.

## Usage

``` r
# S3 method for class 'brma'
dfbetas(
  model,
  type = "mods",
  standardized_coefficients = FALSE,
  transform_factors = TRUE,
  return_loo_estimates = FALSE,
  ...
)
```

## Arguments

- model:

  a fitted brma object.

- type:

  type of parameters to be summarized. Defaults to `"mods"` (for the
  effect size and meta-regression coefficients). The other options are
  `"scale"` (for the heterogeneity and scale-regression coefficients)
  and `"bias"` (for omega, PET, and PEESE publication-bias parameters).

- standardized_coefficients:

  whether to show standardized meta-regression coefficients. Defaults to
  `FALSE`. When set to `TRUE`, standardized meta-regression coefficients
  are returned for the intercept and continuous predictors. These
  coefficients correspond to the standardized scale on which prior
  distributions are specified by default (i.e.,
  `standardize_continuous_predictors = TRUE`).

- transform_factors:

  whether to transform factors to their original names. Defaults to
  `TRUE`.

- return_loo_estimates:

  whether to return the leave-one-out coefficient estimates used to
  compute DFBETAS instead of standardized DFBETAS values. Defaults to
  `FALSE`.

- ...:

  additional arguments (currently ignored).

## Value

If `return_loo_estimates = FALSE`, a data frame with \\K\\ rows
(observations) and \\P\\ columns (parameters), containing DFBETAS
values. If `return_loo_estimates = TRUE`, returns the corresponding
leave-one-out coefficient estimates. Row names correspond to study
labels (if available) or indices.

## Details

This function computes DFBETAS values using the Leave-One-Out (LOO)
approximation based on Pareto Smoothed Importance Sampling (PSIS)
weights. Ideally, DFBETAS is defined as: \$\$DFBETAS\_{ij} =
\frac{\hat{\beta}\_j -
\hat{\beta}\_{j(-i)}}{SE(\hat{\beta}\_{j(-i)})}\$\$ where
\\\hat{\beta}\_j\\ is the estimate of the \\j\\-th coefficient using the
full data, \\\hat{\beta}\_{j(-i)}\\ is the estimate when observation
\\i\\ is omitted, and \\SE(\hat{\beta}\_{j(-i)})\\ is the standard error
of the coefficient when observation \\i\\ is omitted.

In the Bayesian context using LOO approximation:

- \\\hat{\beta}\_{j(-i)}\\ is estimated as the importance sampling
  weighted mean of the posterior samples, using PSIS weights
  \\w\_{is}\\.

- \\SE(\hat{\beta}\_{j(-i)})\\ is estimated as the importance sampling
  weighted standard deviation of the posterior samples.

This approximation allows computing influence statistics without
refitting the model \\K\\ times, making it computationally efficient.
For `type = "bias"`, fixed identification parameters (e.g., the
reference \\\omega = 1\\ interval) can have zero LOO posterior standard
deviation. These parameters are retained in the output, but their
DFBETAS values are reported as `NaN` because the standardized diagnostic
is undefined.

Note: This function requires that LOO-CV has been computed for the model
using
[`add_loo`](https://fbartos.github.io/RoBMA/reference/add_loo.brma.md).

## See also

[`add_loo`](https://fbartos.github.io/RoBMA/reference/add_loo.brma.md),
[`loo_weights.brma`](https://fbartos.github.io/RoBMA/reference/loo_weights.brma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")
  fit <- bPET(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
  fit <- add_loo(fit)

  inf <- dfbetas(fit)
  plot(inf[, 1], type = "h")
}
} # }
```
