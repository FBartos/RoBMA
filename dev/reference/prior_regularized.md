# Regularized Regression Priors

Set ridge, lasso, or horseshoe priors for moderator or scale regression
coefficients. A named group in `prior_mods` or `prior_scale` can share
one shrinkage scale across its terms. The intercept retains its own
prior. For a centered factor, lasso shrinks the length of the fitted
contrast vector and horseshoe assigns one local scale to that vector.
Treatment and independent contrasts have coordinate-wise shrinkage. The
scale is in fitted coefficient units, so the same scale can give
different prior spreads for factor-level effects under different
contrasts. An interaction inherits its factors' fitted codings. To
change the coding of a factor in the interaction, specify it on that
factor's main effect.

## Usage

``` r
prior_regularized(
  family,
  scale,
  terms = NULL,
  slab_scale = NULL,
  slab_df = 4,
  contrast = NULL
)
```

## Arguments

- family:

  One of `"ridge"`, `"lasso"`, and `"horseshoe"`.

- scale:

  A positive number fixing the global scale, or a proper scalar
  continuous prior with positive support from which to learn it. The
  scale is in fitted coefficient units. For a single lasso coefficient
  it is the Laplace scale. For a centered factor it scales the Euclidean
  norm in the group-lasso density.

- terms:

  Optional character vector of fixed-effect term labels. Place a group
  under its own name in the formula prior list. Omitting `terms` selects
  the non-intercept terms without a prior of their own. Separate named
  groups have independent scales, even if they use the same scale prior.
  A multi-column term uses its fitted coefficient basis.

- slab_scale:

  Positive scale of the Student-t slab for the regularized horseshoe.
  The default `NULL` selects the classical horseshoe.

- slab_df:

  Positive degrees of freedom of the Student-t slab.

- contrast:

  Optional factor coding, one of `"treatment"`, `"independent"`,
  `"meandif"` or `"orthonormal"`. With `NULL`,
  [`JAGS_formula()`](https://fbartos.github.io/BayesTools/reference/JAGS_formula.html)
  uses treatment coding by default and a calling model can supply its
  factor coding.

## Value

With numeric scale and omitted terms, ridge and lasso return ordinary
scalar priors accepted by
[`prior_spike_and_slab()`](https://fbartos.github.io/BayesTools/reference/prior_spike_and_slab.html)
and
[`prior_mixture()`](https://fbartos.github.io/BayesTools/reference/prior_mixture.html).
Otherwise the result declares a coefficient group for
[`JAGS_formula()`](https://fbartos.github.io/BayesTools/reference/JAGS_formula.html).
[`is.prior.regularized()`](https://fbartos.github.io/BayesTools/reference/is.prior.html)
identifies group declarations.
[`is.prior.regularized_term()`](https://fbartos.github.io/BayesTools/reference/is.prior.html)
identifies coefficient priors bound to a group.

## See also

[`prior_specification()`](https://fbartos.github.io/RoBMA/dev/reference/prior_specification.md),
[`BayesTools::prior_regularized()`](https://fbartos.github.io/BayesTools/reference/prior_regularized.html)

## Examples

``` r
global <- prior("cauchy", list(0, 0.2), list(0, Inf))
moderator_priors <- list(
  shrinkage = prior_regularized("horseshoe", terms = c("age", "dose"),
                               scale = global, slab_scale = 2, slab_df = 4)
)
# Jointly shrink the contrasts of a design factor in a meta-regression.
design_priors <- list(
  design_group = prior_regularized("lasso", terms = "design", scale = 0.3,
                                   contrast = "meandif")
)
```
