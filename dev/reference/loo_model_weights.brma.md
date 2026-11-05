# Predictive Model Weights for brma Objects

Compute stacking or pseudo-BMA weights from LOO results already stored
in multiple `brma` objects.

## Usage

``` r
# S3 method for class 'brma'
loo_model_weights(
  x,
  ...,
  unit = "estimate",
  method = c("stacking", "pseudobma"),
  optim_method = "BFGS",
  optim_control = list(),
  BB = TRUE,
  BB_n = 1000,
  alpha = 1,
  r_eff_list = NULL,
  cores = getOption("mc.cores", 1)
)
```

## Arguments

- x:

  a `brma` model object (the first model to weight).

- ...:

  additional `brma` model objects or RoBMA-targeted `psis_loo` objects.

- unit:

  output/deletion unit used when extracting stored LOO results from
  `brma` objects.

- method:

  weighting method; `"stacking"` or `"pseudobma"`.

- optim_method, optim_control:

  optimization method and controls for stacking. See
  [`loo_model_weights`](https://mc-stan.org/loo/reference/loo_model_weights.html).

- BB, BB_n, alpha:

  Bayesian-bootstrap controls for pseudo-BMA weighting. See
  [`loo_model_weights`](https://mc-stan.org/loo/reference/loo_model_weights.html).

- r_eff_list:

  optional relative effective sample sizes. Retained for the upstream
  interface and ignored when stored `psis_loo` objects are supplied.

- cores:

  number of cores passed to
  [`loo::loo_model_weights()`](https://mc-stan.org/loo/reference/loo_model_weights.html).

## Value

A named numeric vector containing one predictive weight per model.

## Details

Each `brma` object must first be updated with
`object <- add_loo(object, unit = unit)`. All models must use the same
outcome data, deletion `unit`, and likelihood target. RoBMA validates
these targets before delegating to
[`loo::loo_model_weights()`](https://mc-stan.org/loo/reference/loo_model_weights.html).

Model names are derived from the supplied model expressions. See
[`loo_model_weights`](https://mc-stan.org/loo/reference/loo_model_weights.html)
for details about stacking, pseudo-BMA, and pseudo-BMA+.

## See also

[`add_loo`](https://fbartos.github.io/RoBMA/dev/reference/add_loo.brma.md),
[`loo.brma`](https://fbartos.github.io/RoBMA/dev/reference/loo.brma.md),
[`loo_compare.brma`](https://fbartos.github.io/RoBMA/dev/reference/loo_compare.brma.md),
[`loo_model_weights`](https://mc-stan.org/loo/reference/loo_model_weights.html)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")

  fit_bias <- RoBMA(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
  fit_nobias <- BMA(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")

  fit_bias   <- add_loo(fit_bias)
  fit_nobias <- add_loo(fit_nobias)

  loo_model_weights(fit_bias, fit_nobias)
}
} # }
```
