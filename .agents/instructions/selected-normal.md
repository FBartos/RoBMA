# Selected-Normal and Posterior Routing

Use this guide when changing posterior evaluation, publication-bias selection
kernels, PIT/CDF/PDF/RNG paths, or effect-direction handling.

## Effect Direction

JAGS returns `mu`, `theta`, and random effects on the original fitted scale,
including when `effect_direction = "negative"`. Do not flip those posterior
samples again. `tau` is nonnegative and is never flipped.

PET/PEESE sampling-bias offsets follow effect direction: add them for positive
effects and subtract them for negative effects. Reuse
`.evaluate.brma.bias_offset()`.
Their coefficients are zero in inactive product-space branches, so applying
these offsets to all rows preserves branch routing.

Ordinary normal CDF code owns its direction transformation. Selected-normal
code receives the sign through `.selection_context()`. Callers must not add
another flip.

## Product-Space Routing

Predictive variance consumes the kernel's direct centered variance. Raw second
moments remain available to their callers, but are not subtracted from squared
large means. Outer sampling mixtures and scalar moment integrals use a Gaussian
location center and retain the existing numerical error and transition checks.

For model-averaged publication-bias fits, posterior `bias_indicator` identifies
the active bias branch. Never infer branch identity from `omega` values.

Build one `.selection_context()` for a posterior sample set and pass it to
selected-normal PDF, CDF, RNG, moment, threshold, and plotting helpers. The
context contains branch routing and numerical kernel inputs, including
`bias_indicator`, `use_normal`, and `kernel_mode`.

`.extract_use_normal()` is the supported convenience extractor:

- `TRUE`: the posterior row came from a non-selection branch and uses the
  ordinary normal kernel mode.
- `FALSE`: the row uses its selected-normal kernel.

Do not pass loose `use_normal`, `omega`, or sign arguments through higher
level callers when the context already owns them.

Retained-intercept qCMDE conditioning may operate within fixed product-space
branches. Resolve each branch's focal prior and preserve its indicators and
random-effect gates. Zero retained scales use the ordinary conditional evaluator;
the information cap remains based on the complete requested posterior row set.

Post-fit native kernels may coalesce adjacent, exactly equal step weights on
a validated partition, remapping observed bins to their union intervals. Preserve
the caller's context and branch identity; the compact partition only selects an
equivalent numerical kernel and retains its existing error diagnostics.

## Vector Selection Models

Estimate deletion under sampling conditioning retains the declared remaining
sampling errors and random-source context in both selected and actual normal
branches. Gaussian cancellation conditions sampling and integrated random
coordinates separately. Cluster joint deletion keeps its distinct law, including
singleton clusters; no prediction augmentation policy enters these score keys.

Weightfunction priors carry `selection_model()` with these defaults:
`estimate_random_effects = "integrate"`, `other_random_effects = "condition"`,
`known_sampling_variance = "integrate"`, `weight_rule = "product"` and
`group = NULL`. Each source field accepts `"condition"` or `"integrate"`;
the other weight rule is `"best"`. Constructors take one
`selection = selection_model(...)` input for generated weightfunction priors.
Explicit priors retain their own complete model.

Omitted or NULL `random` in any .mv constructor means a fixed-effect model:
there is no implicit `tau` or random-effects mixture, and both random-effect
selection axes are inapplicable. The whole sampling-error axis still applies.
Heterogeneity priors and scale formulas require an explicit random structure.
Specialized univariate models retain their implicit heterogeneity and clustered
`tau`/`rho`/`I2` behavior.

BayesTools' `random_effects_level_roles()` resolves estimate/other roles from
declared grouping metadata. In .mv models, allow zero or one declared factor
whose grouping levels are one-to-one with retained estimate rows; more than one
is a hard error. All other terms use `other_random_effects`. Keep each family
intact; do not split coefficient bases or infer roles from numerical loadings.
Non-diagonal known group `R` does not change an estimate factor's role. Preserve
its full covariance and dependencies; only diagonal cases use scalar algebra.

For unclustered univariate models, ordinary heterogeneity is the estimate
source. In the specialized clustered interface, within-cluster heterogeneity
is the estimate source and the shared cluster effect is the other source.
Preserve `tau`/`tau2`, `rho` and `I2` summaries. Default estimate variation is
integrated; explicitly conditioning on it defines a different target.

Known sampling variation is the complete error `e ~ N(0, V)`, equivalently
`e = S z`, `z ~ N(0, R)` and `V = S R S` for positive sampling standard
errors. `condition` retains the complete realization during selection
normalization and contributes zero sampling covariance to candidate variation.
`integrate` contributes full `V`. This applies to `vi`, `sei`, univariate and
diagonal inputs too. No residual sampling innovation is always integrated.
Keep `V` authoritative: numerical factorizations do not define statistical
source roles, and sampling coordinates never become estimate random effects.

Sampling a backend node is not a statistical conditioning decision. Retained
sources enter the candidate mean, with their population laws outside selection
normalization. Integrated sources enter candidate covariance with their full
dependencies. Backend auxiliary Gaussian variables may represent either role;
consume resolved source metadata rather than infer conditioning from monitoring.
Product normalizers factor only for conditionally independent
rows. Best uses the weight at the smallest actual p-value, including nonmonotone
weights; thresholds use original sampling SEs.
Non-unit observation `weights` require sampling `integrate` and the supported
independent-product weighted path. Unit weights are identical to omission and
remain valid with sampling `condition`.

Fitted integrated latent prediction uses the normalized powered Gaussian
augmentation documented in the predictive-target contract, conditional on the
unchanged supplied global and retained-source posterior row. Publication weights
at the observed response are constant in that conditional; likelihood powers
instead change its conditioning variance to `sei²/weights`. Population response
sampling keeps the raw variance and existing PET/PEESE direction. Existing
weighted admission refusals remain authoritative.

Publication grouping is active only for `weight_rule = "best"`. Product
branches ignore `group`, do not bind its column or check its partition, and
accept the default `group = NULL` for every supported model structure.
Product integration follows the full integrated covariance dependencies.
For best branches, bind publication groups after the common row selection,
using the explicit column, supported constructor cluster, or an existing
unambiguous singleton fallback. Active selection branches share source settings;
only best branches must share a publication partition. Best requires integrated
dependencies to remain within each group; conditioned sources may connect groups.
Ordinary `metafor::vcalc()` outputs are covariance matrices; they do not provide
publication identities. Do not infer publications from random groups or covariance
blocks. Partial vectors, deletion, prediction and zplots retain the full original
selection event. Under sampling `condition`, estimate deletion integrates the
deleted sampling error conditionally on the retained errors of remaining
estimates. Fixed choices add no parameters.

Keep `known_v_factor()` as an explicit exact computational representation,
and accept ordinary covariance matrices directly. Preserve row identity and
bind publication groups for best rules separately from covariance input.
Equivalent covariances with the same source settings and, for best rules,
publication partition define the same selection model. Do not attach conditioning
meaning or decomposition warnings to diagonal-plus-factor splits.

All-conditioned sources with positive weights almost surely reduce to the
ordinary Gaussian law and provide no likelihood information about the weights.
Zero acceptance on a positive-probability set of retained contexts is invalid;
never discard those contexts. Singular candidate covariance is not itself an
invalid observed model. Do not add jitter or residual noise. Ordinary one-shot
population publication selection integrates all sources. Selection source
settings and prediction `conditioning_depth` remain independent.

The unreleased statistical exact/approximate API and superseded field names
are replaced outright, with no aliases or object adapters. Existing scenario
labels are not selection specifications: inspect each call's explicit settings
and current defaults. Old scenario evidence does not certify a changed target.
Changed roles/workloads require new evidence before equivalence claims.
Integration diagnostics remain separate from model sensitivity. Relative
weights identify neither absolute publication probability nor missing count.

Dense product normalizers may use common/group factor quadrature with a checked
covariance approximation. Monotone weights permit covariance-order bounds;
two-sided and other nonmonotone weights require a perturbation bound. The
absolute bound uses the step weights' total variation; strictly positive
weights additionally permit a relative Gaussian covariance bound including
emitted SD/loading reconstruction error. Combine relative covariance and
quadrature errors with their denominator cross term. Zero weights retain the
absolute bound. Never use agreement of auxiliary normalizers as a covariance
bound for nonmonotone weights. Covariance error,
quadrature change, and exterior Gaussian mass share the existing error budget;
unresolved evaluations retain QMC and its diagnostics.

Deterministic acceptance also requires selection-transition resolution in the
actual standardized Gaussian axes: changed finite cutoff centers use the
supplied signed mean, selection SE and loading; their bracketing node gaps must
not exceed residual SD divided by absolute loading. Coupled axes use maximum
spacing over |z|<=3 against each axis's minimum applicable width. Zero loadings,
equal adjacent weights and analytically cancelled factors add no constraint;
singleton child variance absorbed exactly adds no child axis. Exterior factors
remain owned by existing tail gates. Check the returned candidate, preserving
actual previous-rule values/offsets without requiring resolved comparands.
Precheck the finest available/affordable rule and retain existing fallback or
unavailability, truthful measured diagnostics and unchanged budgets.

Known limits of this rule: it checks each axis, not the product of coincident
transitions, which sharpens with their number (80 identical rows with loadings
(.496, .542), residual SD .247, mean 2.13 below the cutoff and weights (1, .01)
are accepted deterministically with relative error 1.04e-2). The randomized
fallback is not calibrated for such products either: forced onto it, the same
rows passed both its MCSE and change diagnostics for 11 of 40 design seeds, all
with errors of 2-6%. Replacing the rule needs a criterion validated against
both (review 35).

Retained-latent density has separate population-forward, conditional response-
root interpolation and direct conditional routes. Sampling density uses its
actual E|Y affine axis; tails/moments use E's population axis, only when those
components are requested. The scalar native deletion ladder has its own guard.
Context numerator, denominator and outer rules remain separate, with their
common denominator and original raw mass intact. These transition checks are
additional heuristics for primitive A even inside 1/A; they do not establish
an inverse-integrand or focal-Gaussian error bound. Keep all existing Gaussian
basis, covariance, interpolation, compression, mass, exterior and change gates.
Factor-only rank-two event partition/quantile integration is outside this rule.

## Unsupported Paths

Publication-bias weight functions are unsupported for binomial and Poisson
GLMMs.

P-hacking selected-normal branches remain internal and unsupported. Do not
expose or activate them without an explicit maintainer decision and a defined
numerical contract.

## Relevant Files

- `R/evaluate.R`: posterior extraction and bias offsets.
- `R/selection-mapping.R`: selection specifications, contexts, and row routing.
- `R/selection-likelihood.R`: resolved source roles and candidate covariance.
- `R/known-v-representation.R`: authoritative sampling covariance storage.
- `R/pdf*.R`, `R/rng.R`: outcome evaluation.
- `src/selnorm/`: shared native kernels.
- `tests/testthat/helper-selection-kernel.R`: focused test fixtures.

Test positive and negative effect directions, single and product-space models,
and normal versus active selection rows. Use analytic or independent references
for tail calculations.
