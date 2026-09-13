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

Ordinary normal CDF code owns its direction transformation. Selected-normal
code receives the sign through `.selection_context()`. Callers must not add
another flip.

## Product-Space Routing

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

## Vector Selection Models

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

Bind publication groups after the common row selection, using the explicit
column or supported constructor cluster. Ordinary `metafor::vcalc()` outputs
are covariance matrices; they do not provide publication identities. Do not
infer publications from random groups or covariance blocks. Partial
vectors, deletion, prediction and zplots retain the full original selection
event. Under sampling `condition`, estimate deletion integrates the deleted
sampling error conditionally on the retained errors of remaining estimates.
Fixed choices add no parameters.

Keep `known_v_factor()` as an explicit exact computational representation,
and accept ordinary covariance matrices directly. Preserve row identity and
bind publication groups separately from covariance input. Equivalent covariances with
the same publication partition define the same selection model; do not attach
conditioning meaning or decomposition warnings to diagonal-plus-factor splits.

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

## Unsupported Paths

P-hacking selected-normal branches remain internal and unsupported. Do not
expose or activate them without an explicit maintainer decision and a defined
numerical contract.

## Relevant Files

- `R/evaluate.R`: posterior extraction and bias offsets.
- `R/selection-mapping.R`: selection specifications, contexts, and row routing.
- `R/selection-likelihood.R`: resolved source roles and candidate covariance.
- `R/known-v-representation.R`: authoritative sampling covariance storage.
- `R/pdf*.R`, `R/cdf.R`, `R/rng.R`: outcome evaluation.
- `src/selnorm/`: shared native kernels.
- `tests/testthat/helper-selection-kernel.R`: focused test fixtures.

Test positive and negative effect directions, single and product-space models,
and normal versus active selection rows. Use analytic or independent references
for tail calculations.
