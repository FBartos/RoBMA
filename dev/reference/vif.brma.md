# Variance Inflation Factors for brma Objects

Computes variance inflation factors (VIF) and generalized VIF (GVIF) for
a fitted brma meta-regression model. Also optionally returns the
posterior correlation matrix of regression coefficients.

## Usage

``` r
# S3 method for class 'brma'
vif(
  object,
  posterior_correlation = TRUE,
  max_samples = Inf,
  ...,
  probs = NULL,
  include_samples = FALSE
)
```

## Arguments

- object:

  a fitted brma object with moderators

- posterior_correlation:

  logical; whether to also compute and return the posterior correlation
  matrix of regression coefficients. Defaults to `TRUE`.

- max_samples:

  maximum posterior draws used for the VIF/GVIF computation in every
  model class. Defaults to `Inf`. Finite values deterministically thin
  draws across posterior row order. The optional posterior coefficient
  correlation still uses all draws.

- ...:

  additional arguments (currently ignored)

- probs:

  optional credible-interval probabilities for GVIF (ordinary VIF for
  one-df terms). Defaults to `NULL`, which retains the point-summary
  table. Requested empirical quantiles appear as `CI_<probability>`
  columns. These intervals describe GVIF, not its adjusted root.

- include_samples:

  whether to include a `samples` data frame with posterior `draw` row
  indices before thinning, `term_index` (the row of the `vif` table,
  also distinguishing repeated display labels), `term`, `df`, `GVIF`,
  and `GVIF^(1/(2*df))`. Defaults to `FALSE`.

## Value

An object of class `vif.brma` containing:

- vif:

  A data frame with columns `term`, `df`, `GVIF`, and `GVIF^(1/(2*df))`
  (= \\GVIF^{1/(2 \cdot df)}\\). For single-df terms, GVIF equals the
  standard VIF.

- posterior_correlation:

  A correlation matrix of posterior regression coefficient samples when
  requested and available; otherwise `NULL`.

- samples:

  The optional long data frame of draw-wise diagnostics when
  `include_samples = TRUE`. Printing and data-frame coercion of the
  result retain the summary-table view; extract this element explicitly
  for draws.

## Details

VIF is computed from the correlation matrix derived from the coefficient
variance-covariance matrix \\(X'WX)^{-1}\\. For standard meta-regression
models, \\W = \mathrm{diag}(w_i/(v_i + \tau^2))\\ at each posterior
heterogeneity draw. Scale and multilevel models use observation-specific
\\\tau_i\\ and block-structured covariance where applicable. VIF/GVIF
and its adjusted root are computed in each draw and averaged separately.
These posterior means are generally different from applying VIF to an
averaged covariance or to a posterior-mean heterogeneity. The mean
adjusted GVIF is likewise generally different from transforming the mean
GVIF. The diagnostic retains its GLS design-based meaning; it is not
calculated from the empirical posterior covariance of the coefficients.

A VIF of 1 indicates no collinearity; values above 5 or 10 are commonly
considered problematic.

The members of a
[`RoBMA()`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.md)
ensemble can differ in their bias terms (PET adds a `PET` column, PEESE
a `PEESE` column, and the null and weight-function members add none), so
each posterior draw uses the design of its member. A location term
summarizes all selected posterior draws; `PET` and `PEESE` summarize
only draws in their own class (`NA` when no selected draw belongs to
it). Their samples are `NA` outside that class, not zero, and their
credible intervals are conditional on class membership. Location terms
keep the specified location design, including terms with spike-and-slab
priors.

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

Data-only objects created with `only_data = TRUE` are rejected because
VIF depends on fitted heterogeneity information, not just the parsed
data.

For
[`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
known-`V` models, VIF is computed from the full marginal GLS covariance.
Formula random effects are marginalized through BayesTools'
random-effect covariance metadata rather than treated as conditioned
fitted effects.

For known-`V` models, the VIF requires applying the inverse of the
marginal covariance. A covariance that is too ill-conditioned for that
(a reciprocal condition number of about `2e-10` or less in any posterior
draw used) is unavailable: the call stops with an error of class
`RoBMA_known_v_precision_unavailable`, whose parent class is
`RoBMA_numerical_unavailable`, instead of returning a rounded number.
The condition carries the draw positions (`rows`) and their
`reciprocal_condition`.

For generalized linear mixed models, this is a working effect-size-scale
diagnostic based on the fitted model's approximate sampling variances.
It is not a VIF derived directly from the binomial or Poisson
information matrix. When the marginal covariance is fixed, the
calculation reduces to the conventional fixed-covariance VIF used by
metafor.

## References

Fox J, Monette G (1992). “Generalized collinearity diagnostics.”
*Journal of the American Statistical Association*, **87**(417), 178–183.
[doi:10.1080/01621459.1992.10475190](https://doi.org/10.1080/01621459.1992.10475190)
.

## See also

[`regplot()`](https://fbartos.github.io/RoBMA/dev/reference/regplot.md),
[`summary.brma()`](https://fbartos.github.io/RoBMA/dev/reference/summary.brma.md)

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
