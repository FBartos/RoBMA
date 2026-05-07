# Estimated Marginal Means for brma Objects

Computes estimated marginal means for a fitted `brma` object with
moderators using
[`BayesTools::as_marginal_inference()`](https://fbartos.github.io/BayesTools/reference/as_marginal_inference.html).

## Usage

``` r
# S3 method for class 'brma'
marginal_means(
  object,
  null_hypothesis = 0,
  normal_approximation = FALSE,
  n_samples = 10000,
  output_measure = NULL,
  transform = NULL,
  bf = NULL,
  ...
)
```

## Arguments

- object:

  a fitted `brma` object with moderators.

- null_hypothesis:

  point null hypothesis used for inclusion Bayes factors. Defaults to
  `0`.

- normal_approximation:

  whether prior and posterior density at the null should be approximated
  with a normal distribution. Defaults to `FALSE`.

- n_samples:

  number of samples/grid points used by BayesTools for marginal prior
  densities. Defaults to `10000`.

- output_measure:

  effect-size measure for location/effect predictions. Defaults to the
  fitted measure. Supported conversions are among `"SMD"`, `"COR"`,
  `"ZCOR"`, and `"OR"`; `"RR"`, `"HR"`, `"IRR"`, `"RD"`, and `"GEN"` can
  only be returned on their fitted measure. Use `transform = "EXP"` for
  ratio-scale output from log-scale measures.

- transform:

  optional display transformation. Currently `"EXP"` exponentiates
  log-scale measures `"OR"`, `"RR"`, `"HR"`, and `"IRR"`.

- bf:

  whether inclusion Bayes factors should be shown by default in
  summaries. Defaults to `TRUE` for RoBMA/BMA objects and `FALSE` for
  single-model `brma` objects.

- ...:

  additional arguments (currently ignored).

## Value

A list of class `marginal_means.brma` containing the BayesTools
`marginal_inference` object and parameter metadata.

## See also

[`summary()`](https://rdrr.io/r/base/summary.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html),
[`summary.brma()`](https://fbartos.github.io/RoBMA/reference/summary.brma.md),
[`regplot()`](https://fbartos.github.io/RoBMA/reference/regplot.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE) &&
    requireNamespace("metafor", quietly = TRUE)) {
  data(dat.bcg, package = "metadat")
  dat <- metafor::escalc(
    measure = "RR",
    ai      = tpos,
    bi      = tneg,
    ci      = cpos,
    di      = cneg,
    data    = dat.bcg
  )

  fit <- brma(yi = yi, vi = vi, mods = ~ alloc, data = dat, measure = "RR")
  mm  <- marginal_means(fit)
  summary(mm)
  plot(mm, parameter = "alloc")
}
} # }
```
