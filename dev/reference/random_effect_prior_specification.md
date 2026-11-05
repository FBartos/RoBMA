# Random-effect prior specification

Random-effect prior distributions are used by
[`brma.mv`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
when a `random` formula is supplied. The formula and the prior have
deliberately separate roles. The formula defines the random-effect
design, grouping variables, and covariance family. The
`prior_heterogeneity` argument defines the prior scale and can
additionally customize correlations, variance allocation, monitoring,
prediction, and computational parameterization.

RoBMA constructs a complete random-effect prior when
`prior_heterogeneity` is omitted. An ordinary positive prior
distribution changes only the base random-effect standard-deviation
prior. A
[`prior_random`](https://fbartos.github.io/RoBMA/dev/reference/prior_random.md)
object can either customize selected non-scale settings while inheriting
RoBMA's scale defaults, or explicitly own the complete scale
architecture.

## Arguments

- random:

  an optional formula or list of formulas specifying BayesTools
  random-effect terms for
  [`brma.mv`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md).
  See
  [`random_effect_formula_tags`](https://fbartos.github.io/RoBMA/dev/reference/random_effect_formula_tags.md)
  for the formula syntax, component naming rules, and public semantic
  parameter names.

- prior_heterogeneity:

  prior distribution for residual heterogeneity. In random-formula
  [`brma.mv`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
  models, omission constructs all defaults; an ordinary positive prior
  distribution replaces the base or aggregate random-effect SD prior
  while retaining automatic allocation and correlation defaults; and a
  [`prior_random`](https://fbartos.github.io/RoBMA/dev/reference/prior_random.md)
  object supplies partial or complete random-effect settings as
  described below. With a `scale` model, the row-shaped SD source of a
  random block must be bound to the block that its `scale` model
  targets; a source bound to another block stops with an error of class
  `RoBMA_prior_scale_source_target` (also of class `RoBMA_prior_error`).

- scale:

  an optional formula or named list of formulas specifying
  location-scale predictors. With random-effect formulas, `scale` models
  the row-wise total random-effect standard deviation consumed by a
  random-effect component.

- prior_scale:

  prior distribution for scale-regression coefficients. A single prior
  applies to all terms; a named list can specify term-specific priors.
  If omitted or `NULL`, default priors are used.

- prior_unit_information_sd:

  numeric. The unit information standard deviation used to construct the
  default base heterogeneity prior distribution.

- rescale_priors:

  numeric. A scaling factor for supported prior distributions. It
  rescales the base heterogeneity prior distribution before RoBMA
  constructs the random-effect prior specification.

## Details

### Division of responsibility

The random-effect setup has three layers:

|  |  |  |
|----|----|----|
| Input | Controls | Does not redefine |
| `random` | blocks, grouping, coefficient/index basis, covariance family | prior scale or prior family |
| ordinary `prior_heterogeneity` | base or aggregate random-effect SD prior | formula structure or automatic allocation |
| `prior_random(...)` | block overrides, covariance priors, allocation, and policies | the structure parsed from `random` |

Plain `(expr | group)` syntax and `us()` / `un()` define unstructured
random-coefficient blocks. `||` and
[`diag()`](https://rdrr.io/r/base/diag.html) define independent
random-coefficient blocks, while `id()` gives independent coefficients
one shared SD. In these families, `1`, `0`, and `-1` control the random
intercept; continuous slopes, factor slopes, and interactions are
supported. Factor coding comes from the resolved contrast metadata, not
from intercept syntax. Use `random_block(contrasts = ...)` when a
block-specific factor basis is required.

The `cs()`, `hcs()`, `ar1()` /
[`ar()`](https://rdrr.io/r/stats/ar.html), `har()`, and `car()` families
instead own a structure-specific index basis. Their left side supplies
index variables rather than a coefficient formula, so these structures
reject intercept controls and block contrast overrides. See
[`random_effect_formula_tags`](https://fbartos.github.io/RoBMA/dev/reference/random_effect_formula_tags.md)
for the complete formula and variable-type contract.

Nested grouping syntax, such as `random = ~ 1 | study / esid`, creates
multiple blocks. A formula list creates top-level components:


    random = list(
      study   = ~ 1 | study,
      outcome = ~ 1 | outcome
    )

A bare formula or unnamed one-entry list has no redundant public
component prefix. An explicitly named one-entry list retains its name.
In lists with two or more entries, missing names become `component 1`,
`component 2`, and so on. Component names are used by summaries and by
component-specific `scale` models. Block overrides in
[`prior_random()`](https://fbartos.github.io/RoBMA/dev/reference/prior_random.md)
must match the resolved block names, which can be inspected with
`print_prior(object, parameter = "random")`.

### Default standard-deviation prior

RoBMA first constructs the same base heterogeneity SD prior described in
[`prior_specification`](https://fbartos.github.io/RoBMA/dev/reference/prior_specification.md).
Under the standard UISD defaults, \$\$\sigma\_{\mathrm{base}} \sim
\mathrm{Normal}^{+}(0, \mathrm{UISD}/4).\$\$ The scale multiplier can be
changed through `RoBMA.options(default_UISD.heterogeneity = ...)`. The
same `measure`, `ni`, `prior_unit_information_sd`, informed-prior, and
`rescale_priors` rules used for ordinary heterogeneity apply before the
prior is assigned to the random-effect structure.

Supplying an ordinary positive prior through `prior_heterogeneity`
replaces this base prior only. RoBMA still determines whether it applies
directly to one SD, to the RoBMA quantity `tau_total` for a
total-variance allocation, or to `tau_common` for a mean-variance
allocation.

### Structure-specific defaults

Let \\K\\ denote the resolved number of random-coefficient columns for a
coefficient structure or the number of index levels for an index
structure. The formula family determines the number of marginal SDs and
the default correlation prior:

|  |  |  |  |
|----|----|----|----|
| Structure | Basis | Marginal SDs | Default correlation prior |
| `id()` | coefficients | one shared SD | none |
| [`diag()`](https://rdrr.io/r/base/diag.html) or `||` | coefficients | one per column or factor unit | none |
| plain `|`, `us()` / `un()` | coefficients | one per column or factor unit | `LKJ(1)` matrix |
| `cs()` | index levels | one shared SD | `Uniform(-1 / (K - 1), 1)` on raw `cor` |
| `hcs()` | index levels | one per level | `Uniform(-1 / (K - 1), 1)` on raw `cor` |
| `ar1()` / [`ar()`](https://rdrr.io/r/stats/ar.html) | ordered index | one shared SD | `Uniform(-1, 1)` on raw `cor` |
| `har()` | ordered index | one per level | `Uniform(-1, 1)` on raw `cor` |
| `car()` | numeric coordinates | one shared SD | `Uniform(0, 1)` on raw `cor`, supported on `[0, 1)` |

A mean-difference or orthonormal factor term (the default coding of
factor moderators in model-averaging constructors) is one exchangeable
unit: its contrast columns share one SD, `tau(f)`, and are uncorrelated,
so the model of the level deviations does not depend on the arbitrary
contrast basis. Under `us()` the unit keeps its correlations with the
other coefficients, reported per level (`rho(intercept, f[a])`) and as
the multiple correlation `rho_total(intercept, f)`; a `us()` block can
hold at most one such term with three or more levels. SD-component
allocations count units, not columns: \\K\\ below is the number of SD
components.

Correlation priors are added after the formula basis and \\K\\ are
resolved, and only when the block has more than one random coefficient.
The CS/HCS lower bound is the exact positive-definiteness bound, so the
default includes all admissible negative correlations. `LKJ(1)` is
uniform over correlation matrices and does not favor positive
correlations. In particular, `hcs(index | group)` has one common
pairwise correlation, while `us(0 + index | group)` can estimate a
different correlation for every pair of resolved coefficient columns
outside an exchangeable factor unit.

### Automatic variance allocation

RoBMA uses symmetric `Dirichlet(1, ..., 1)` priors whenever a variance
scale must be divided. The allocation is hierarchical when the formula
contains multiple top-level components, nested blocks, or heterogeneous
blocks:

- A single homogeneous block receives the base SD prior directly.

- A split across blocks or components uses `scale = "total_variance"`.
  For child weight \\w_j\\, \$\$\sigma_j =
  \sigma\_{\mathrm{total}}\sqrt{w_j}, \qquad \sum_j \sigma_j^2 =
  \sigma\_{\mathrm{total}}^2.\$\$ The public aggregate quantities are
  `tau_total` and `tau2_total`; components are `tau2_prop(...)`.

- A heterogeneous block, such as DIAG, US, HCS, or HAR with \\K \> 1\\,
  adds an SD-component allocation with `scale = "mean_variance"`:
  \$\$\sigma_j = \sigma\_{\mathrm{common}}\sqrt{K w_j}, \qquad
  K^{-1}\sum_j \sigma_j^2 = \sigma\_{\mathrm{common}}^2.\$\$ The public
  aggregate quantities are `tau_common` and `tau2_common`; components
  are `tau2_mult(...)` and `tau_mult(...)`.

Without component-specific scale models, the base prior controls a
shared total or average variance scale. Component-specific or partial
`scale` models instead create independent roots for the targeted
components; every unscaled component retains its own base SD prior.
Splits remain within each component, so their combined variance is not
governed by one global total prior.

### Public random-effect parameter names

BayesTools stores concrete coordinates, semantic quantities, and aliases
in one versioned parameter map. RoBMA summaries, plots, density
estimation, and hypotheses use its semantic catalog view. BayesTools
canonical random-effect names have the form
`(formula) owner: quantity(arguments)`. RoBMA consistently maps these to
meta-analytic names: a sole intercept is `tau` for a bare block and
`study: tau` for a named block, while a non-intercept coefficient
remains explicit, for example `study: tau(x)`. The unique owner-free
shorthand `tau` is accepted when only one component matches. Bare
formulas and unnamed one-entry lists omit a redundant owner; explicitly
named one-entry lists and multi-component models retain owners for
disambiguation. Other examples are
`rho(group[sensitivity],group[specificity])`, `tau_total`, `tau_common`,
and `tau2_prop(study)`. A custom
[`random_variance_allocation()`](https://fbartos.github.io/RoBMA/dev/reference/prior_random.md)
always requires a stable internal `name`; its `display_name` and
`component_names` separately control these public owner and component
labels.

Formula-random correlations are public `rho` quantities. BayesTools'
general `cor` quantities, compact scalar correlation coordinates, and
LKJ primitives remain internal to RoBMA. In the specialized `cluster`
interface, `rho` remains the public variance-allocation and induced
intracluster correlation parameter described in
[`prior_specification`](https://fbartos.github.io/RoBMA/dev/reference/prior_specification.md).

### Ways to customize random-effect priors

There are three increasingly explicit customization levels:

1.  Supply an ordinary positive prior distribution to
    `prior_heterogeneity`. This changes the base SD prior but preserves
    all automatic allocations and default correlation priors.

2.  Supply a partial
    [`prior_random()`](https://fbartos.github.io/RoBMA/dev/reference/prior_random.md)
    object containing only correlation, contrast, monitoring, new-level
    prediction, or parameterization settings. RoBMA supplies its base SD
    prior and every required variance allocation.

3.  Supply a
    [`prior_random()`](https://fbartos.github.io/RoBMA/dev/reference/prior_random.md)
    object containing scale architecture. The object is complete and
    RoBMA does not add or merge scale defaults.

A
[`prior_random()`](https://fbartos.github.io/RoBMA/dev/reference/prior_random.md)
object is partial only when all scale-defining fields are absent. The
following fields make it complete: top-level or block-specific `sd`;
`random_covariance(sd = ...)`; block `sd_source`; block `terms`
containing term-specific SD overrides; or `allocation`. Once any of
these is present, the user owns the entire scale architecture and must
provide every required SD or allocation source.

For example, the following partial specification changes only the factor
basis of block `study`:


    random_prior <- prior_random(
      study = random_block(contrasts = c(group = "independent"))
    )

Combined with `random = ~ us(0 + group | study)`, this creates one
correlated random coefficient per `group` level. RoBMA supplies the
UISD-scaled common-SD prior and mean-variance Dirichlet allocation;
BayesTools adds the default `LKJ(1)` correlation prior after resolving
the coefficient count. No empty
[`random_covariance()`](https://fbartos.github.io/RoBMA/dev/reference/prior_random.md)
or manual allocation is needed.

A correlation-only customization is partial as well. Use
`random_covariance(cor = prior_lkj(...))` for US/UN and
`random_covariance(cor = ..., cor_scale = ...)` for CS/HCS, AR1/HAR, or
CAR. An explicitly supplied scalar `cor` prior defaults to the Fisher-z
scale; `cor_scale = "cor"` instead places it directly on the raw
correlation, and `cor_scale = "logit"` uses a logit transform of the
structure's admissible interval. The `structure` in the formula remains
authoritative, and a conflicting covariance override is rejected.

Complete specifications can define global SD defaults, block-specific
[`random_block`](https://fbartos.github.io/RoBMA/dev/reference/prior_random.md)
settings, term-specific SDs, external SD sources, and allocations built
with
[`random_variance_allocation`](https://fbartos.github.io/RoBMA/dev/reference/prior_random.md).
Use a complete specification when the desired scale hierarchy differs
from RoBMA's automatic total/mean-variance allocation; do not add an SD
merely to customize a correlation or contrast.

The `parameterization` argument and block-specific overrides select
`"noncentered"`, `"centered"`, or `"auto"` backend representations
without changing the specified standard-deviation, allocation, or
correlation prior distributions. Explicit centering requires strictly
positive backend-owned standard-deviation coordinates; BayesTools
rejects unsupported centered combinations and lets `"auto"` fall back to
the noncentered representation. Monitoring and new-level prediction
settings similarly affect saved quantities and prediction policy, not
the statistical prior. Keep latent effects monitored when bridge
sampling or conditional prediction is needed. See
[`prior_random`](https://fbartos.github.io/RoBMA/dev/reference/prior_random.md)
for every constructor argument and complete examples.

### Location-scale random-effect models

When `scale` is used together with `random`, the scale model targets the
row-wise total random-effect standard deviation. With one top-level
random-effect component, a single `scale` formula is unambiguous. With
multiple top-level random-effect components, `scale` must be a named
list whose names identify the targeted components or concrete blocks.
Each target has its own row-wise SD root and its allocations split
within that component; untargeted components retain independent base SD
roots. The `prior_scale` argument controls the scale-regression
coefficients. A custom
[`prior_random()`](https://fbartos.github.io/RoBMA/dev/reference/prior_random.md)
containing explicit SD architecture must agree with these row-shaped
scale sources; incompatible sources are rejected rather than silently
replaced.

### Inspecting the resolved priors

The completed prior depends on the resolved formula columns and index
levels. Use `only_priors = TRUE` to perform data, formula, and prior
resolution without MCMC, then inspect the result with
`print_prior(object, parameter = "random")`. The printed object shows
the resolved block names, SD sources, allocation hierarchy, correlation
priors, contrasts, and policies that will be used for fitting. Its
mathematical `sigma_total` and `sigma_common` labels describe the prior
architecture; fitted RoBMA selectors use `tau_total` and `tau_common`,
respectively.

## See also

[`prior_specification`](https://fbartos.github.io/RoBMA/dev/reference/prior_specification.md),
[`random_effect_formula_tags`](https://fbartos.github.io/RoBMA/dev/reference/random_effect_formula_tags.md),
[`brma.mv`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md),
[`prior_random`](https://fbartos.github.io/RoBMA/dev/reference/prior_random.md),
[`random_block`](https://fbartos.github.io/RoBMA/dev/reference/prior_random.md),
[`random_covariance`](https://fbartos.github.io/RoBMA/dev/reference/prior_random.md),
[`random_variance_allocation`](https://fbartos.github.io/RoBMA/dev/reference/prior_random.md)

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- data.frame(yi = c(0.1, 0.2, 0.3, 0.4),
                  ni = rep(40, 4), study = rep(c("a", "b"), each = 2),
                  group = factor(rep(c("control", "treatment"), 2)))
V <- diag(rep(0.04, nrow(dat)))
group_prior <- prior_factor("normal", list(0, 1), contrast = "independent")
random_prior <- prior_random(
  study = random_block(contrasts = c(group = "independent"))
)

priors <- brma.mv(
  yi = yi, V = V, ni = ni,
  mods = ~ 0 + group,
  random = ~ us(0 + group | study),
  prior_mods = group_prior,
  prior_heterogeneity = random_prior,
  data = dat, measure = "GEN",
  only_priors = TRUE
)
print_prior(priors, parameter = "random")
} # }
```
