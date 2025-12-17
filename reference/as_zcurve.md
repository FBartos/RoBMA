# Transform RoBMA object into z-curve

This function transforms the estimated RoBMA model into a z-curve object
that can be further summarized and plotted. Only available for
normal-normal models estimated using the spike-and-slab algorithm (i.e.,
`algorithm = "ss"`). See Bartoš and Schimmack (2025) and
[[`vignette("ZCurveDiagnostics", package = "RoBMA")`](articles/ZCurveDiagnostics.md)](doc/ZCurveDiagnostics.md)
for more detail.

## Usage

``` r
as_zcurve(x, significance_level = stats::qnorm(0.975), max_samples = 1000)
```

## Arguments

- x:

  A RoBMA object

- significance_level:

  Significance level used for computation of z-curve estimates.

- max_samples:

  Maximum number of samples from the posterior distribution that will be
  used for estimating z-curve estimates.

## Value

`as_zcurve` returns a list of tables of class 'zcurve_RoBMA'.

## Examples

``` r
if (FALSE) { # \dontrun{
# using the example data from Anderson et al. 2010 and fitting the default model
# (note that the model can take a while to fit)
fit <- RoBMA(r = Anderson2010$r, n = Anderson2010$n,
             study_names = Anderson2010$labels, algorithm = "ss")

zcurve_fit <- as_zcurve(fit)
summary(zcurve_fit)
} # }
```
