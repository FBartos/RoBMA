# Variance Inflation Factors for brma Objects

Computes variance inflation factors (VIF) and generalized VIF (GVIF) for
a fitted brma meta-regression model. Also optionally returns the
posterior correlation matrix of regression coefficients.

## Usage

``` r
# S3 method for class 'brma'
vif(object, posterior_correlation = TRUE, ...)
```

## Arguments

- object:

  a fitted brma object with moderators

- posterior_correlation:

  logical; whether to also compute and return the posterior correlation
  matrix of regression coefficients. Defaults to `TRUE`.

- ...:

  additional arguments (currently ignored)

## Value

An object of class `vif.brma` containing:

- vif:

  A data frame with columns `term`, `df`, `GVIF`, and `GVIF^(1/(2*df))`
  (= \\GVIF^{1/(2 \cdot df)}\\). For single-df terms, GVIF equals the
  standard VIF.

- posterior_correlation:

  A correlation matrix of posterior regression coefficient samples when
  requested and available; otherwise `NULL`.

## Details

VIF is computed from the correlation matrix derived from the coefficient
variance-covariance matrix \\(X'WX)^{-1}\\. For standard meta-regression
models, \\W = \mathrm{diag}(w_i/(v_i + \hat\tau^2))\\ with \\\hat\tau\\
equal to the posterior mean heterogeneity. For scale and multilevel
models, the coefficient covariance is averaged across posterior
heterogeneity draws, using observation-specific \\\tau_i\\ and
block-structured multilevel covariance where applicable.

A VIF of 1 indicates no collinearity; values above 5 or 10 are commonly
considered problematic.

For multi-column terms, such as factor contrasts, the Generalized VIF
(GVIF) of Fox and Monette (1992) is reported. GVIF captures the joint
inflation for all coefficients belonging to the same term. To enable
comparison across terms with different degrees of freedom, \\GVIF^{1/(2
\cdot df)}\\ is also reported; this value can be compared against the
usual VIF thresholds (after squaring).

When `posterior_correlation = TRUE`, the function also returns the
posterior correlation matrix of the regression coefficients. This
Bayesian diagnostic complements VIF: while VIF diagnoses the *potential*
for collinearity problems (a data property), the posterior correlation
shows the *realized* identification given the data and priors.
Informative priors can mitigate collinearity, reducing posterior
correlations even when VIF is high.

## References

Fox J, Monette G (1992). “Generalized collinearity diagnostics.”
*Journal of the American Statistical Association*, **87**(417), 178–183.
[doi:10.1080/01621459.1992.10475190](https://doi.org/10.1080/01621459.1992.10475190)
.

## See also

[`regplot()`](https://fbartos.github.io/RoBMA/reference/regplot.md),
[`summary.brma()`](https://fbartos.github.io/RoBMA/reference/summary.brma.md)

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
  fit <- brma(
    yi      = yi,
    vi      = vi,
    mods    = ~ ablat + year,
    data    = dat,
    measure = "RR",
    seed    = 1,
    silent  = TRUE
  )

  vif(fit)
}
} # }
```
