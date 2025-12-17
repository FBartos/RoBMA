# Extract method for Robust Bayesian Meta-Analysis Fits

`residuals.RoBMA` extract residuals based on the RoBMA model. Only
available for normal-normal models estimated using the spike-and-slab
algorithm (i.e., `algorithm = "ss"`).

## Usage

``` r
# S3 method for class 'RoBMA'
residuals(
  object,
  conditional = FALSE,
  output_scale = NULL,
  probs = c(0.025, 0.975),
  as_samples = FALSE,
  ...
)
```

## Arguments

- object:

  a fitted RoBMA object

- conditional:

  show the conditional estimates (assuming that the alternative is
  true). Defaults to `FALSE`. Only available for `type == "ensemble"`.

- output_scale:

  transform the meta-analytic estimates to a different scale. Defaults
  to `NULL` which returns the same scale as the model was estimated on.

- probs:

  quantiles of the posterior samples to be displayed. Defaults to
  `c(.025, .975)`

- as_samples:

  whether posterior samples instead of a summary table should be
  returned. Defaults to `FALSE`.

- ...:

  additional arguments

## Value

`pooled_effect` returns a list of tables of class 'BayesTools_table'.

## See also

[`predict.RoBMA()`](reference/predict.RoBMA.md)

## Examples

``` r
if (FALSE) { # \dontrun{
require(metafor)
dat <- escalc(measure = "OR", ai = tpos, bi = tneg, ci = cpos, di = cneg,
              data = dat.bcg)

# fit meta-regression
robma_dat <- data.frame(
  logOR  = dat$yi,
  se     = sqrt(dat$vi),
  ablat  = dat$ablat,
  alloc  = dat$alloc
)

fit <- NoBMA.reg(~ ablat + alloc, data = robma_dat,
                 seed = 1, algorithm = "ss", parallel = TRUE)

residuals(fit)

} # }
```
