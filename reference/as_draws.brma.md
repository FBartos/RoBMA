# Convert brma Objects to posterior Draws Formats

Provides an interface to the posterior package for `brma` objects. These
functions convert the MCMC samples from a fitted `brma` model to various
draws formats supported by the posterior package, enabling the use of
posterior's rich set of diagnostics and summary functions.

## Usage

``` r
as_draws(x, ...)

as_draws_array(x, ...)

as_draws_df(x, ...)

as_draws_list(x, ...)

as_draws_matrix(x, ...)

as_draws_rvars(x, ...)

# Default S3 method
as_draws(x, ...)

# Default S3 method
as_draws_array(x, ...)

# Default S3 method
as_draws_df(x, ...)

# Default S3 method
as_draws_list(x, ...)

# Default S3 method
as_draws_matrix(x, ...)

# Default S3 method
as_draws_rvars(x, ...)

# S3 method for class 'brma'
as_draws(x, include_auxiliary = FALSE, ...)

# S3 method for class 'brma'
as_draws_array(x, include_auxiliary = FALSE, ...)

# S3 method for class 'brma'
as_draws_df(x, include_auxiliary = FALSE, ...)

# S3 method for class 'brma'
as_draws_list(x, include_auxiliary = FALSE, ...)

# S3 method for class 'brma'
as_draws_matrix(x, include_auxiliary = FALSE, ...)

# S3 method for class 'brma'
as_draws_rvars(x, include_auxiliary = FALSE, ...)
```

## Arguments

- x:

  an object to convert. The `brma` methods expect a fitted `brma`
  object; default methods forward non-`brma` objects to the
  corresponding posterior conversion function.

- ...:

  additional arguments passed to the corresponding posterior function.

- include_auxiliary:

  logical; whether to include the hidden intermediate parameters (see
  Details) in addition to the parameters shown by the default
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) of the fit.
  Formula-random parameters retain their RoBMA semantic names in either
  case.

## Value

An object of the corresponding posterior draws class.

## Details

By default, the draws hold exactly the parameters that the default
[`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) of the fit show (the
location, heterogeneity, and moderator coefficients, the multilevel
allocation `rho`, and the publication-bias parameters), together with
the model-inclusion indicators behind the reported inclusion
probabilities. The hidden intermediate parameters are added by
`include_auxiliary = TRUE`: the latent cluster and estimate effects
(`gamma[k]` and the GLMM `theta[i]`), the selection-model latent
variables (`eta`), the GLMM nuisance base rates and rates (`pi[i]` and
`phi[i]`), the random-effect variances, component SDs derived from a
variance allocation, and total variance of formula-random models (for
example `tau2` or `tau2_total`), and the other backend variables listed
below.

These functions are S3 generics. Their default methods forward to the
corresponding posterior generics so attaching RoBMA preserves the usual
posterior behavior for non-`brma` objects.

The following conversion functions are available:

- `as_draws`: converts to the default draws format

- `as_draws_array`: converts to a 3-D array (iteration x chain x
  variable)

- `as_draws_df`: converts to a data frame with columns for iterations,
  chains, and variables

- `as_draws_list`: converts to a list of lists

- `as_draws_matrix`: converts to a 2-D matrix (draw x variable)

- `as_draws_rvars`: converts to random variable objects

These methods require the posterior package to be installed. For `brma`
methods, conversion is performed by first extracting the MCMC samples as
a `mcmc.list` object and then using the corresponding posterior
conversion function. By default, backend-private latent effects,
known-`V` dependency factors, random-effect simulation and Cholesky
factors, and prior-parameterization variables are omitted. Stable fitted
coordinates such as model indicators remain available. Ordinary
spike-and-slab backend slab and inclusion coordinates are omitted by
default; `include_auxiliary = TRUE` retains them. Semantic random-effect
indicators remain available as described below. Formula-random
coordinates are replaced by their semantic RoBMA quantities using the
same names as summaries, plots, and hypotheses, such as `tau_total`,
`tau_common`, `rho(...)`, and `tau2_prop(...)`. A `rho(...)` column can
contain `NA` draws: an original-scale correlation is undefined in draws
where one of its standard deviations is zero, e.g. where an inclusion
gate or a point prior sets it to zero. Summaries and hypothesis tests of
the correlation use its defined draws; summaries computed from the
exported draws need to handle the `NA` values. Random-effect SD switches
are exposed as exact columns named after the standard deviation they act
on, with an `_indicator` suffix, pairing with the semantic SD column the
same way `mu_x_indicator` pairs with `mu_x`. A heterogeneity-component
inclusion gate is 0/1, such as `study: tau_indicator` and
`observation: tau_indicator`, and carries the inclusion state the
summary reports under `Heterogeneity:`. An allocation-SD prior
indicator, such as `tau_total_indicator`, records which mixture
component allocated the total standard deviation and therefore takes one
value per component, matching the order of that prior's components.
Ordered predictor terms return their fitted totals and normalized
sampled allocation shares by default, replacing their deterministic
increment coordinates. Shares are normalized from the sampled Gamma
primitives and remain defined when a total is zero or its inclusion
event is off. Totals with point priors are exact constants for every
fitted slice; fixed allocations do not create share columns. The labels
come from the fitted BayesTools ordered specification, including reached
factor levels and interaction slices. These quantities remain on the
fitted predictor scale. Ordered-total mixture indicators and inclusion
probabilities remain available, while nested slab and pseudo-prior
coordinates are omitted. With `include_auxiliary = TRUE`, raw ordered
increments and eligible total auxiliaries are also retained. Internal
Gamma coordinates are never exported. `as_draws_rvars()` preserves
allocation labels as scalar random variables and keeps ordinary vector
and multi-slice total indices. If an ordered sampled parameter conflicts
with an existing fitted column, conversion stops with classes
`RoBMA_ordered_draws_conflict` and `RoBMA_ordered_draws_unavailable`,
with the message
`Ordered sampled parameter '<name>' conflicts with existing fitted draws. Refit the model with the current RoBMA/BayesTools build.`
Missing required ordered sources stop with
`RoBMA_ordered_draws_unavailable` and
`BayesTools_ordered_coordinates_unavailable`, with parent class
`BayesTools_ordered_unavailable` and the message
`Ordered source draws are unavailable. Refit the model with the current RoBMA/BayesTools build.`
Entirely fixed numeric-point terms need only the existing chain row
count to materialize their declared totals. Both RoBMA ordered refit
errors also inherit `RoBMA_refit_required` and
`BayesTools_refit_required`. Structural point-prior coordinates are
returned as exact constant columns, with chain timing taken from the
fitted draw geometry. Private backend anchors are never exposed. Use
`include_auxiliary = TRUE` to skip RoBMA auxiliary filtering; added
backend variables retain their coordinate names, while public
random-effect parameters retain their RoBMA names. `brma_samples`
objects have separate methods documented at
[`as_draws.brma_samples`](https://fbartos.github.io/RoBMA/reference/as_draws.brma_samples.md).

## See also

[`draws`](https://mc-stan.org/posterior/reference/draws.html),
[`brma`](https://fbartos.github.io/RoBMA/reference/brma.md),
[`as_draws.brma_samples`](https://fbartos.github.io/RoBMA/reference/as_draws.brma_samples.md)
