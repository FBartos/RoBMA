# Model Interfaces and BayesTools Integration

Use this guide for model constructors, random-effect prior completion,
parameter metadata, summaries, selectors, and public result interfaces.

## Formula ownership

BayesTools owns formula semantics and compiled metadata. When its source is
available, consult its `jags-formula.md` guide for the detailed shared contract.
RoBMA consumes the resolved column counts, bases, levels, covariance ownership,
and labels. Do not introduce a separate parser or reconstruct a basis with
`stats::model.matrix()`, fixed formulas, or posterior draws.

Factor contrasts belong to `prior_factor()` and stored metadata, independently
of intercept syntax. Thus `mods = ~ 0 + group` with an independent factor
prior means a structural zero intercept and one coefficient per group level.
In random-coefficient blocks, `random_block(contrasts = ...)` overrides the
default reuse of a fixed factor's concrete contrasts; `0 + group` alone does
not force level indicators.

Keep random-coefficient structures (`id()`, `diag()`, `us()`/`un()`) separate
from structure-owned index specifications (`cs()`/`hcs()`,
`ar1()`/`ar()`/`har()`, `car()`). Coefficient formulas support intercept controls,
slopes, and interactions. Plain bars mean US; double bars mean diagonal.
ID shares one SD, diagonal uses independent column-specific SDs, and US adds
an unrestricted correlation matrix. Index structures reject intercept controls
and block contrast overrides. CS/HCS combines discrete columns by observed
interaction; AR1/HAR orders one discrete column by factor levels or sorted
values; CAR uses one finite numeric coordinate or ordered factor with numeric
labels and actual distances. Discrete indices may be factor, character,
numeric/integer, or logical. HCS has level-specific SDs and one common pairwise
correlation; it is not US in a level-indicator basis. Never relabel
coefficient-basis covariance as level-basis covariance.

Treat `prior_random()` with no global/block SD, SD source, term-specific SD
override, covariance-owned SD, or variance allocation as a partial override.
Complete its scale with RoBMA's ordinary UISD/allocation rules while preserving
the user's contrasts, covariance priors, and policies. When any scale
architecture is supplied, do not merge additional scale defaults. Correlation
defaults remain BayesTools-owned: omitted US/UN uses `LKJ(1)` and omitted scalar
correlations use the complete structure-specific raw interval.

## Parameter metadata and naming

BayesTools stores one versioned `parameter_map()` with linked coordinate,
quantity, and alias tables. `parameter_coordinates()` exposes backend rows
keyed by `coordinate_name`; `parameter_catalog()` exposes semantic quantities
keyed by `canonical_name`. Resolve public selectors through the catalog and
obtain draws through `parameter_draws()`.

Internal LKJ primitives, compact scalar-correlation coordinates, allocation
weights, and covariance-construction dependencies remain coordinate-only.
Monitoring them does not make them public aliases or allow them through
plotting, density, or hypothesis gates.

Random-effect density and marginal-diagnostic fast paths consume BayesTools'
`random_effects_marginal_update_plan()`. Route metadata-declared affine families
through exact covariance algebra and retain the generic evaluator for factor,
Markov, and unsupported families. Do not infer affineness from posterior draws
or evaluated covariance candidates.

BayesTools canonical names are `(formula) owner: quantity(arguments)`.
Parentheses hold coefficient/parameter names; square brackets hold factor/index
levels. RoBMA maps `sd`/`var`/`cor` at its I/O boundary to `tau`/`tau2`/`rho`,
including `tau_total`, `tau2_total`, `tau_common`, `tau2_common`,
`tau2_prop(...)`, `tau_mult(...)`, and `tau2_mult(...)`. Simplification removes
only a sole `intercept` argument and omits an owner only when resolution stays
unique. Keep other arguments explicit, for example
`study: rho(group[sensitivity],group[specificity])`.

A bare random formula or unnamed one-entry list omits a redundant component
prefix; a named one-entry list retains its name. Lists with multiple entries
generate missing names as `component 1`, `component 2`, etc. Generated
allocations keep a stable internal `name`, use `display_name = ""` when no
public owner is needed, and carry public `component_names` separately.

Ordinary and specialized `brma(..., cluster = ...)` models use `tau`/`tau2`
(plus specialized `rho`/`I2`). `brma.mv()` uses `tau`/`tau2` for one component,
`tau_total`/`tau2_total` for a genuine additive aggregate, and
`tau_common`/`tau2_common` for a mean-variance allocation scale.

Printed summaries and their data-frame exports follow the univariate layout.
Without `mods` or `scale`, combine location and random estimates in `Estimates`.
With predictors, random estimates join `Common Estimates`; location coefficients
use `Meta-Regression` without scale predictors and `Location` with them. Scale
coefficients retain their own table. Without location moderators, `mu` stays in
the common table. Apply the same layout to conditional estimates while retaining
the separate raw summary tables and their parameter metadata.
Merge random-effect inclusion rows into the displayed `Component Inclusion`
table with `Random: ` labels, before `Publication Bias`. Preserve bounded BF
markers and diagnostics when combining and ordering rows. Data-frame exports
follow this layout; raw inclusion subtables remain available.

## Classes and public arguments

`brma` is the base class; `brma.norm`, `brma.glmm`, and `brma.mv` specialize the
likelihood/model. `RoBMA` extends `brma` for product-space model averaging.
Wrappers prepend their class, for example `c("BMA.norm", "RoBMA", "brma")` and
`c("BMA.glmm", "RoBMA", "brma.glmm", "brma")`.

Primary constructors are `brma()`, `brma.norm()`, `brma.glmm()`, `brma.mv()`,
`bselmodel()`, `bPET()`, `bPEESE()`, `RoBMA()`, `BMA()`, `BMA.norm()`, and
`BMA.glmm()`.

Follow existing S3 result families. For method-specific results, prefer
`<generic>.<input_class>` with matching print/summary methods and retain useful
underlying classes. In BayesToolsVerse, the shared public API guide owns
diagnostic-message and tabular-coercion conventions.

Use `component` to distinguish post-fit model parts. Normalize `"mods"` and
`"location"` through shared helpers; `"bias"` selects publication-bias
parameters. Retain released `parameter_mods` and `parameter_scale` plotting
arguments through 4.x. Reserve `type` for output or prediction kind.

If shared implementation, documentation, and tests conflict, repair the
BayesTools-owned behavior there; avoid downstream compatibility workarounds.
