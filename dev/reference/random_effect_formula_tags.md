# Random-Effect Formula Structure Tags

[`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
random-effect formulas use BayesTools random-effect formula structure
tags to describe the covariance structure of random coefficients. These
tags are formula specials parsed by
[`BayesTools::random_effects_formula()`](https://fbartos.github.io/BayesTools/reference/random_effects_formula.html),
not ordinary user-facing functions.

## Details

Use tags as bare names inside the `random` formula supplied to
[`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md),
for example `random = ~ cs(comp | study)`. They are parsed from the
formula call and are not called as normal R functions, so users should
not write `RoBMA::cs()` or `BayesTools::cs()`. RoBMA forwards the
formula to
[`BayesTools::random_effects_formula()`](https://fbartos.github.io/BayesTools/reference/random_effects_formula.html),
stores the parsed BayesTools random-effect metadata, and uses that
metadata for JAGS compilation, prior matching, prediction, and summary
labels.

BayesTools distinguishes random-coefficient formulas from structured
index specifications. For `id()`,
[`diag()`](https://rdrr.io/r/base/diag.html), and `us()` / `un()`, the
left side of `|` is a coefficient formula. It can contain an intercept,
continuous or factor slopes, and interactions. As in an ordinary
random-slope formula, `1`, `0`, and `-1` include or remove the random
intercept. Plain `(expr | group)` defaults to `us(expr | group)`, while
`expr || group` is the diagonal shorthand.

Factor coding inside a random-coefficient block is not selected by
removing the intercept. BayesTools reuses an already established
fixed-factor contrast by default and lets a
[`random_block`](https://fbartos.github.io/RoBMA/dev/reference/prior_random.md)
set an explicit block-specific basis through `contrasts`. For example,
`us(0 + group | study)` together with
`random_block(contrasts = c(group = "independent"))` creates one
correlated random coefficient per `group` level and no random intercept.
Without the independent contrast override, `0 + group` preserves the
resolved contrast basis instead of forcing raw level indicators.

In contrast, `cs()`, `hcs()`, `ar1()` /
[`ar()`](https://rdrr.io/r/stats/ar.html), `har()`, and `car()` own a
structure-specific index basis. Their left side names index columns
rather than coefficients. These tags reject explicit `1`, `0`, and `-1`,
and their blocks do not accept `random_block(contrasts = ...)`.

Supported structure tags:

- `id()`: independent random-coefficient columns with one shared
  standard deviation.

- [`diag()`](https://rdrr.io/r/base/diag.html): independent
  random-coefficient columns with one standard deviation per column.
  `||` is shorthand for this structure.

- `us()` / `un()`: one standard deviation per random-coefficient column
  and an unrestricted correlation matrix. Plain `|` syntax defaults to
  this structure.

- `cs()`: compound symmetry over one or more discrete index columns;
  index levels share one standard deviation and one common pairwise
  correlation.

- `hcs()`: heterogeneous compound symmetry over one or more discrete
  index columns; levels have separate standard deviations and one common
  pairwise correlation.

- `ar1()` / [`ar()`](https://rdrr.io/r/stats/ar.html): homogeneous
  autoregressive order-1 covariance over exactly one discrete index;
  levels share one standard deviation and correlations decay as
  \\\mathrm{cor}^{\|i-j\|}\\ over the resolved order.

- `har()`: heterogeneous autoregressive order-1 covariance; index levels
  have separate standard deviations and AR(1) correlations.

- `car()`: continuous-time autoregressive covariance; use exactly one
  finite numeric/integer coordinate or an ordered factor with numeric
  level labels, as in `car(time | study)`. Correlations are
  \\\mathrm{cor}^d\\ based on actual absolute distances, and
  \\\mathrm{cor}\\ is constrained to \\\[0, 1)\\.

- `random()` / `re()`: explicit random-effect block wrapper. By default
  this uses an unstructured covariance, and it can be useful when
  specifying an explicit `name` for prior matching and summaries.

Under [`diag()`](https://rdrr.io/r/base/diag.html), `us()`, and plain
bars, the contrast columns of a mean-difference or orthonormal factor
term are one exchangeable unit: they share one standard deviation and
are uncorrelated with each other, while `us()` keeps the unit's
correlations with the other columns. A `us()` block that is a single
such unit therefore has no correlation parameter, and a `us()` block
cannot hold two such terms with three or more levels.

The discrete structured tags `cs()`, `hcs()`, `ar1()`, and `har()`
accept factor, character, numeric/integer, or logical index columns.
Existing factor level order is retained; otherwise sorted unique values
determine the levels and, for AR(1), their order. `cs()` and `hcs()` may
combine multiple columns with `+`; BayesTools forms the observed
interaction in lexicographic order. `ar1()` and `har()` require exactly
one index column. AR1/HAR warns when numeric labels in a character or
unordered factor resolve to a non-increasing numeric order, such as
`1, 10, 2`, while preserving that order. Use numeric values for numeric
ordering, or an ordered factor with explicit levels to confirm another
intended order.

`hcs(index | group)` and `us(0 + index | group)` are different models
even when both produce one standard deviation per displayed level. HCS
constrains every pair to the same correlation. US estimates an
unrestricted correlation matrix among the selected coefficient columns
and its factor basis must be specified through contrast metadata when it
matters. Alias tags are normalized before model compilation: `un()` to
`us()`, [`ar()`](https://rdrr.io/r/stats/ar.html) to `ar1()`, and
unqualified `random()` / `re()` to `us()` unless a covariance override
changes the structure. In these wrappers,
`covariance = "cs", hom = FALSE` resolves to `hcs()`, and
`covariance = "ar1", hom = FALSE` resolves to `har()`. Any explicit
covariance structure in
[`prior_random()`](https://fbartos.github.io/RoBMA/dev/reference/prior_random.md)
must match the resolved structure. Direct `cs()` / `ar1()` tags retain
their homogeneous meaning.

[`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
also accepts metafor-style known group covariance through its `R` and
`Rscale` arguments. This covariance is attached to random-effect
grouping levels and is distinct from the known sampling covariance `V`.
RoBMA translates `R` internally to BayesTools random-effect group
covariance metadata; users should continue to use `R`/`Rscale` in
[`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
rather than BayesTools-native `group_covariance` arguments. Current
support is limited to random-intercept blocks such as
`random = ~ 1 | study`. Supported one-to-one known-`R` blocks in
known-`V` normal models can be compiled as marginalized diagonal
row-variance terms; other supported known-`R` blocks remain sampled.

The scalar-correlation structures `cs()`, `hcs()`, `ar1()` /
[`ar()`](https://rdrr.io/r/stats/ar.html), and `har()` use uniform
defaults over their complete valid raw-correlation intervals:
`(-1 / (K - 1), 1)` for CS/HCS and `(-1, 1)` for AR1/HAR, where `K` is
the number of resolved index levels. CAR uses `Uniform(0, 1)` because
its correlation-decay parameter is one-sided. These defaults therefore
do not silently exclude negative correlations when the structure permits
them.

A list of formulas defines top-level random-effect components:

    random = list(
      study  = ~ cs(comp | study),
      design = ~ cs(comp | design)
    )

A bare formula or unnamed one-entry list has no redundant user-facing
component prefix. An explicitly named one-entry list retains its name.
In a list with two or more formulas, missing names are generated as
`component 1`, `component 2`, and so on. These component names are used
by RoBMA summaries and by component-specific scale models. Explicit
`name =` arguments inside `random()` or `re()` define the final
random-effect block name used by
[`prior_random()`](https://fbartos.github.io/RoBMA/dev/reference/prior_random.md)
block overrides.

### Public parameter names

BayesTools stores concrete coordinates, selectable semantic quantities,
and aliases in one versioned parameter map. Its canonical random-effect
names use `(formula) owner: quantity(arguments)`. RoBMA maps these to
meta-analytic I/O names: a sole random intercept is `tau` for a bare
block and `study: tau` for a named block, while a non-intercept
coefficient remains explicit, for example `study: tau(x)`. The
owner-free `tau` shorthand is accepted only when it resolves uniquely.
Bare formulas and unnamed one-entry lists omit a redundant owner;
explicitly named one-entry lists and multi-component models retain it.

Public correlations are named `rho`. Automatic total-variance
allocations expose `tau_total`, `tau2_total`, and `tau2_prop(...)`.
Mean-variance allocations within heterogeneous blocks expose
`tau_common`, `tau2_common`, `tau2_mult(...)`, and `tau_mult(...)`.
BayesTools' general names and compact covariance-construction
coordinates remain internal to RoBMA. These meta-analytic names are the
supported selectors for summaries, plotting, density estimation, and
hypotheses.

## See also

[`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md),
[`random_effect_prior_specification`](https://fbartos.github.io/RoBMA/dev/reference/random_effect_prior_specification.md),
[`prior_random`](https://fbartos.github.io/RoBMA/dev/reference/prior_random.md),
[`random_covariance`](https://fbartos.github.io/RoBMA/dev/reference/prior_random.md),
[`BayesTools::random_effects_formula()`](https://fbartos.github.io/BayesTools/reference/random_effects_formula.html)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- brma.mv(
  yi      = yi,
  V       = V,
  random  = ~ cs(outcome | study),
  data    = dat,
  measure = "GEN",
  prior_unit_information_sd = 1
)

fit_network <- brma.mv(
  yi = yi,
  V  = V,
  random = list(
    study  = ~ cs(comp | study),
    design = ~ cs(comp | design)
  ),
  data    = dat,
  measure = "OR"
)

group_prior <- prior_factor(
  distribution = "normal",
  parameters   = list(mean = 0, sd = 1),
  contrast     = "independent"
)
random_prior <- prior_random(
  study = random_block(
    contrasts = c(group = "independent")
  )
)
fit_group <- brma.mv(
  yi = yi, V = V, ni = ni,
  mods = ~ 0 + group,
  random = ~ us(0 + group | study),
  prior_mods = group_prior,
  prior_heterogeneity = random_prior,
  data = dat, measure = "GEN",
  prior_unit_information_sd = 1
)
} # }
```
