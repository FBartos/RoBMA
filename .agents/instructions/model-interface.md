# Model Interfaces and BayesTools Integration

Use this guide for model constructors, random-effect prior completion,
parameter metadata, summaries, selectors, and public result interfaces.

## Formula ownership

Named alternative-list intercepts supply the complete prior, while ordinary
term alternative/null entries remain separate components. Reject named null-list
intercepts (including NULL/FALSE) and explicit global declarations left without
a consumer by that override. Resolve supplied argument presence before defaults
or omission filtering. Component scale retains actual global SD and inclusion/
omission-gate consumers, including gate-only roots with local SD sources; unused
scale-intercept defaults need no UISD. Complete or structural-zero location intercepts leave an
omitted effect default unused; missing slope/global SD consumers retain UISD.
Scale intercepts retain the positive SD contract because their model uses
log(intercept). Structural-zero mods retain
their existing point-zero validation and explicit-global compatibility.

Ordered term containers pass through when they are the sole retained component.
Their mixture/inclusion declaration belongs on total; combining a container with
a separate term component refuses with the prior_ordered(total=...) and named
NULL omission remedy. Constructor/binding support does not establish downstream
ordered inclusion/export compatibility with an unaccepted dependency.

Marginal-means BF rows localize declared BayesTools formula-measure and
ordinate refusals to NA with the original measure reason, structured condition
diagnostics, and matching unavailable posterior-density source. Missing-source,
refit, other transform, invalid-context and programming errors remain strict.
An explicit late KDE summary BF request computes skipped BFs from stored
conditional cells at the stored null; existing computed results are retained.

Expanded hypothesis results retain positional statement identity in the
producer's uniform statement-major quantity loop. Identical statement text
does not merge repeated statement numbers or row-owned diagnostics.

Selected and prior-only formula draw carriers retain the actual full prior
declarations with the BayesTools design and scaling. Numeric point readers
require a finite numeric scalar; expression points use BayesTools' owned
monitor/replay route and do not become numeric constants by their wrapper.

BayesTools owns formula semantics and compiled metadata. When its source is
available, consult its `jags-formula.md` guide for the detailed shared contract.
RoBMA consumes the resolved column counts, bases, levels, covariance ownership,
and labels. Do not introduce a separate parser or reconstruct a basis with
`stats::model.matrix()`, fixed formulas, or posterior draws.
Formula-design 6, unscale-design 2 and coefficient-transform 3 require refitting
older formats. Select active targets/recipes before finite C access. Dynamic
numeric targets use original fitted state; prior laws use public BayesTools APIs.
Static qCMDE/IWMDE refuse state-dependent raw maps or stochastic formula products.
Declared measure_unavailable entries use BayesTools_formula_measure_unavailable;
Their optional structured model/log diagnostics survive replacement draw
wrappers and hypothesis refusal conversion. Malformed metadata, refit,
programming and other numerical failures remain strict. RoBMA's empirical
product-space model summaries retain their natural zero meaning. linear_weight_space distinguishes
coefficient recipes from fitted formula_contribution rows. Constant products fold
once into raw fitted weights; contribution weights never pass through C/A again.
Fitted and prior-only plots supply full ordinary multiplier laws. Current compiler
copies complete only through supported public prior/plot entry points; empty
scaling intentionally selects fitted coordinates.

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

Ordered `as_draws()` exports consume `JAGS_ordered_parameter_spec()` source
coordinates and sampled-parameter labels. Read only declared sampled or
structural sources in one targeted materialization, separately from the public
draw schema. Export fitted totals and Gamma-normalized sampled shares, including
zero-total states; never recover shares from increments. Omit deterministic
increments and nested total components by default, retaining declared event
indicators and inclusion probabilities. Auxiliary exports retain raw increments
without exposing unrelated internals. For rvars, map only declared allocation
labels through reversible opaque scalar names; preserve ordinary vector indices.
Entirely fixed numeric-point terms with no available sources reuse the public
chain solely for row count. Guard that exception with declared point totals,
no total event node and fixed allocations; never substitute for missing
stochastic, discrete, mixture or expression sources. A named input-only row
column may supply geometry when the public chain has zero coordinates.
Missing required ordered sources, including partial source loss, raise
`RoBMA_ordered_draws_unavailable`, `BayesTools_ordered_coordinates_unavailable`
and `BayesTools_ordered_unavailable` through the shared refit helper, retaining
`RoBMA_refit_required` and `BayesTools_refit_required` with the supported refit
message.

Internal LKJ primitives, compact scalar-correlation coordinates, allocation
weights, and covariance-construction dependencies remain coordinate-only.
Monitoring them does not make them public aliases or allow them through
plotting, density, or hypothesis gates.

Random-effect density and marginal-diagnostic fast paths consume BayesTools'
`random_effects_marginal_update_plan()`. Route metadata-declared affine families
through exact covariance algebra and retain the generic evaluator for factor,
Markov, and unsupported families. Do not infer affineness from posterior draws
or evaluated covariance candidates.
Covariance q-grid shortcuts certify fixed-mean invariance from the complete
replacement-coordinate dependency metadata, using at most one candidate per row;
dependent or unknown mean updates retain the generic candidate evaluator.

Allocated component SD density targets in predictor-scaled blocks require the
standardized view or KDE for original-scale requests, following the structured
scaled-component capability. The restriction is block-wide, including centered
intercepts with accidental draw parity; direct one-to-one transforms remain supported.

Post-fit mixture likelihoods keep random-source indicators with their persisted
compiled priors, including externally bound SD sources; localize only indicators
whose consuming priors were localized. Compute internally generated point-only
linear targets from declared point-prior locations using the same matrix
transformation as target draws and atom metadata. User-supplied null values
retain exact equality semantics.

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
table with `Heterogeneity: ` labels, matching the univariate `Heterogeneity`
row, before `Publication Bias`. Model-overview components in `summary_models()`
keep their own `Random: ` labels. Preserve bounded BF
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


Contribution-space hypothesis combination specs fold known constant multipliers
once into fitted weights; coefficient-space specs retain their existing units.
Stochastic contribution products and dynamic original combinations retain the
typed static qCMDE/IWMDE refusal and row reason. Point/transitive/group consumers
transport producer log ordinates and canonical logs with aligned diagnostics;
plain exports strip hidden column carriers.
Model and conditional prior contexts use schema 2 with validated model log owners;
primitive contexts use schema 1. Coefficient transformation descriptors remain schema 3.
Variance point numerical diagnostics retain SD evaluation provenance and apply
the common Jacobian to requested variance log heights. Point-region results
reuse the actual point producer's diagnostics in statement order. Variance
unit-part composition detects only the admitted real point logBF -Inf together
with its directional unit infinity before arithmetic. It retains canonical NA
and the existing undefined warning for that numerical-zero pair, preserving
independent region-prior diagnostics and legitimate boundary infinities. Scalar and
exp-affine discovery catches only declared formula-measure refusals; requested
hypotheses, missing owners, refit conditions and programming errors stay strict.
