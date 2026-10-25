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

The per-transition rule requires selection-transition resolution in the
actual standardized Gaussian axes: changed finite cutoff centers use the
supplied signed mean, selection SE and loading; their bracketing node gaps must
not exceed residual SD divided by absolute loading. Coupled axes use maximum
spacing over |z|<=3 against each axis's minimum applicable width. Zero loadings,
equal adjacent weights and analytically cancelled factors add no constraint;
singleton child variance absorbed exactly adds no child axis. Exterior factors
remain owned by existing tail gates. Check the returned candidate, preserving
actual previous-rule values/offsets without requiring resolved comparands.
Precheck the finest available/affordable rule and retain existing fallback or
unavailability, truthful measured diagnostics and unchanged budgets. The rule
checks each axis, not the product of coincident transitions, which sharpens
with their number (80 identical rows with loadings (0.496441, 0.542314),
residual SD 0.246872, mean -2.125758, selection SE 1, cutoff 0 and weights
(1, .01) pass it with relative error 1.04e-2). It decides the weight rule
`best`, the Smolyak sparse grid, the dense envelope ladder, the zplot density
mirrors, the funnel/zcurve cells and the sampling-conditioned scalar ladders.
The deterministic ladders of the factor and cluster kernels under `product`
weights use the shift-amplitude certificate below instead, for fitting nodes,
post-fit batches and the estimate-deletion delegation alike (`tol` is
`relative_tolerance`, default 0.005).

- Shifted rule: moving the nodes of one latent axis by `d` and multiplying its
  weights by `exp(-d x - d^2/2)` is an exact change of variables. The shifted
  rule estimates the same integral and differs from the rung by the aliasing
  error at another lattice phase. Weights are not renormalized. A shifted copy
  builds its own weights and never reads precomputed unshifted ones.
- Shifts and groups: every latent axis is its own group in every structure
  (tensor axes, the root and each leaf of a forest, the axis of a rank-one
  cluster) and is shifted by half and by a quarter of the node gap around the
  rung's own tilted mean of that axis (the gap of the interval containing it;
  the first or last gap outside the nodes, the gap to its right on a node).
  Means and shifted values come from the evaluator that produced the rung
  value. A shift recomputes the rows of its axis and of the axes below it and
  the messages of the axes above it, and reuses every other message and the
  free rows, at any depth (a leaf shift costs its own rows, a root shift its
  tree). Leaves are never shifted together, which hides leaf-leaf aliasing (a
  valid 371-row depth-one forest was accepted that way with relative error
  1.0e-2).
- Amplitude: with `e2 = expm1(q_half - q0)` and `e4 = expm1(q_quarter - q0)`
  per group, `amp = sqrt((e2/2)^2 + (e2/2 - e4)^2)` is the amplitude of the
  first harmonic of that group's aliasing error, sampled at three phases. The
  half shifts of all groups come first; when `sum |e2|/2 > tol/2` the rung is
  rejected without the quarter shifts, which cannot change that decision.
- Acceptance: the ladder starts at its first affordable rung (no precheck or
  exact start) and accepts the first rung that passes, in order, a finite
  value, the successive relative changes (factor: current and previous;
  cluster: current), the exterior-mass tail gate, the sharpness gate,
  `sum amp <= tol/2` and the combined allowance
  `tail_bound + 2 sum amp <= tol/(1+tol)`, with `tail_bound` the relative bound
  the tail gate compares. The accepted rung's own value is returned. A block
  without transitions has amplitude zero and evaluates no shifts. The factor 2
  is the measured envelope of error over amplitude on gated rungs, not a bound.
- Sharpness gate: the amplitude assumes that the first harmonic dominates. An
  edge much narrower than the node gap makes the error sawtooth-like, with
  harmonics the three phases cannot see (native blocks with edge width over
  gap of 0.02-0.1 were accepted with errors up to 2.0e-2). A rung is certified
  only when the per-transition rule resolves its transitions with every step
  width relaxed by `1/0.15`: edge width (residual SD over absolute loading) at
  least 0.15 times the node gap. Gated rungs showed error over amplitude up to
  1.76 in scans, and searches of the whole gated ladder found no accepted
  error above 77% of `tol`. The gate is conservative: sharp blocks go on to the
  next stage even when a deeper rule would be accurate.
- Refinement tier: when the ladder rejects every affordable rung, rank-one
  clusters, tensors of rank at most two and forests of depth at most one try, before Smolyak and
  the randomized fallback, the average of the finest affordable rung over the
  shifts `(m/4) h` per axis, `m = 0..3` and `h` the gap above, which keeps only
  aliasing harmonics that are multiples of 4 (leaves are averaged
  independently per leaf at each root shift, then multiplied). It is accepted
  when this average agrees within `tol/2` with the average over `(m/2) h` and
  with the same average at the previous affordable rung, which must itself pass
  the ladder's edge-width gate (a coarser comparison rung that cannot resolve
  the edges is not trusted to average out aliasing), and the tail gate holds
  for it; the result has relative change the larger check and relative MCSE 0.
  Every evaluation passes the node budget. The tier runs only when `tol` is at
  least 0.005: its gate and its checks were calibrated at the default, and the
  error of the harmonics they cannot see does not shrink with `tol` (at `tol`
  of 1e-3 and 3e-4 walls of about a thousand coincident rows were accepted with
  errors of 1.5 and 1.2 times `tol`); below it, blocks the ladder rejects go on
  to the Smolyak grid and the randomized fallback. The tier's sharpness gate
  asks for an edge width of at least `sqrt(max(1, 2 log N)) / 8` times the node
  gap on each axis, `N` the number of rows with a nonzero loading on that axis
  (all of them, including nearby and coupled rows):
  a wall of `N` coincident steps is narrower than one step by about that
  factor and keeps a harmonic-4 error (up to 44% over `tol`) that neither
  check sees. Averaging eight shifts instead does not remove it. Counting is
  conservative within each axis: nearby or crossing transitions retain the
  full penalty, while a large cluster on another axis adds no penalty here.
  Cluster refinement uses its own finest and preceding rules and shared row
  workspace; successful refinement reports the larger comparison as relative
  change, zero relative MCSE, and a deterministic result.
- Randomized fallback: for `product` weights the proposals are the prior, a
  Gaussian with the tilted mean and full covariance of the finest evaluated
  rung (tensor: weighted moments; depth-one forest: from the messages;
  rank-one cluster: its own finest rung, handed to the factor kernel) and a copy
  with standard deviations times 1.5. They require a tensor, depth-one forest or
  rank-one cluster, at least three finite rungs, a finest rung that passes the
  tail gate (moments of a rule that misses the selected mode say nothing about
  it) and a positive definite covariance; every other block, and `best`, keep
  the optimized Laplace proposals. Product weights retain the unoptimized
  boundary-start locations with unit covariance as well: optimizing every start
  can move all proposals back to a dominant central mode and erase a smaller
  but material joint tail. A further proposal covers intermediate combinations:
  in an orthonormal basis formed from loading-row directions, each coordinate
  independently mixes `N(0,1)` and the unit-variance normals centered at that
  coordinate of every boundary start, equally. Its density is a product of
  univariate mixture densities, so evaluation requires a quadratic rotation
  and a sum over coordinates and boundaries, rather than enumerating their
  exponentially many joint combinations. The basis is
  only a proposal coordinate system; no loading direction is removed from the
  model. Product weights also add a proposal that
  conditions on selection transitions. First absorb every row with zero
  loadings into its exact normalizer product `C`. For the remaining `N` rows,
  write `m = min(omega)`, `F(z) = product_i w_i(z)` for conditional row weights,
  `S(z) = sum_i (w_i(z) - m)`, and `D = E[S(z)]`, calculated from the
  univariate marginal Gaussian bin probabilities. Its density is
  `phi(z) S(z) / D`: select a row/bin with probability proportional to its
  marginal bin probability times `omega_bin - m`, draw that row's outcome
  within the bin, then draw the latent factors conditional on the outcome.
  The categorical choice and within-bin quantile share one uniform by interval
  remapping. The Gaussian proposals remain in the mixture to explore joint
  selection events. The mixture integrates the nonnegative correction
  `C (F(z) - m^N)`; the known baseline `C m^N` is added to each scramble before
  calculating the raw standard error and nested-design change. Positive
  recurrences avoid subtracting nearly equal products. Exactly identical rows
  share bin-probability evaluations; their multiplicities remain in the product,
  the proposal density and its component probabilities. The total point budget
  and stop remain `max(z * relative MCSE, relative change) <= tol`, `z = 2.36`.
  Finite MCSE values remain raw empirical standard errors. A row/bin event
  whose mass relative to the largest underflows even extended precision (a row
  hundreds of SD from its transition) has exactly zero sampling probability: it
  is left out of the draws and of the proposal density together, while the
  target keeps every row. NaN or `+Inf` proposal masses, or no event mass at
  all, instead report infinite diagnostics, so numerical unavailability cannot
  silently become an accepted zero likelihood. At a sampled point, a NaN or
  `+Inf` log target or event density is rejected the same way, a point with
  zero target contributes nothing whatever its event density, and a zero event
  density at a point with positive target is rejected rather than weighted. Exact
  categorical boundaries in valid supplied uniforms use the conditional median;
  this defines the transform on a measure-zero set without clipping interior
  probabilities.
- QMC stopping multiplier: every randomized selection normalizer is accepted
  only when `z` times its relative MCSE (and its relative change, where it
  reports one) is within `tol`, with `z = 2.36`, the 97.5% quantile of
  Student's t with 7 degrees of freedom, the default 8 scrambles, as a fixed
  constant that does not follow `scrambles`. The one constant is
  `selnorm_qmc_z` in `src/selnorm/selnorm-qmc-z.h`, used by the native stops
  (the factor and cluster fallback, the conditioned-sampling loop) and by every
  judge of a returned QMC diagnostic: the JAGS wrappers of the factor, cluster,
  dense covariance and conditioned-sampling normalizers, the dense batch's
  reuse flag, and in R `.selection_joint_factor_check_diagnostics`,
  `.selection_conditioned_sampling_diagnostics`,
  `.selection_joint_dense_loglik_check_mcse`, `.selection_joint_checked_event`
  and `.zplot_normalizer_block`. R reads it through `.selection_qmc_z()`
  (`RoBMA_selnorm_qmc_z`); never copy the number. Reported `relative_mcse`,
  `relative_change` and the values in the rejection messages stay raw, and each
  message names the rule ("accepted when 2.36 times the standard error and the
  change are within 'relative_tolerance'"; native wording
  `SELNORM_QMC_CRITERION_*` in the header, R wording
  `.selection_qmc_criterion()`), because a rejected result can show an MCSE
  below `tol`. Four checks judge their own combined error and are not scaled:
  the zplot selection projection (`.selection_factor_projection`), the zplot
  density batches (`.zplot_joint_block`), the retained-context integration
  (`.zplot_full_event_context_draw`) and the deletion-ratio QMC
  (`.selection_conditioned_sampling_deletion_loglik_chunk`). The zplot
  full-event projection (`.zplot_full_event_projection`) judges its event
  integrals through `.selection_joint_checked_event`, which is scaled.
  Deterministic rules do not use the multiplier.
- Last-evaluation memo: JAGS evaluates a node first at its current value, so
  about 17% of the node calls repeat the previous one. Each factor and cluster
  step node keeps owned copies of every input of a deterministic result (`x`,
  all parameter arrays including the quadrature tables, scalar controls, dims
  and flags) and its log density, and answers a call only when every stored
  value is equal (pointer identity never suffices: JAGS mutates in place). The
  QMC design is not an input. Only results the kernel flags `deterministic`
  (ladders, tier, sparse grid; the batch results carry the flag as an attribute
  of that name) are stored, never a randomized, non-finite or thrown one. The
  memo is a function-local map keyed by node, never deleted, and cleared as a
  whole above 4096 keys. Each entry copies the quadrature tables (about 32 KB
  for a cluster node, 5 KB for a factor node), so the map can hold roughly
  130 MB for the rest of the R session, entries of finished models included;
  a process with more step nodes than the cap clears it repeatedly and gains
  little from it. `ROBMA_SELNORM_NODE_MEMO=0`, read once per process, disables
  it. Draws are bitwise identical with and without it.
- Lean arithmetic: the nested-forest rung with its shifted copies and the
  rank-one cluster rung accumulate node products in double. A row's node
  normalizers are computed into a call-owned workspace and all validated (the
  checks of the batch normalizer routine) before any is multiplied, so a failed
  row leaves the products untouched. An evaluation (a rung or one shifted
  copy) with a failed check, a node product outside `[2^-960, 2^960]` or rows
  whose smallest positive weight cannot stay in that range (skipped up front)
  is repeated unchanged in long double. Values agree to 1e-13 relative. Tensor
  direct sums, the tier, `best` and projection paths keep long double, and no
  scratch is static (the batch entries are threaded).
- Corrected sampler: `selection.sampler = "coarse_corrected"`
  (`src/samplers/CoarseCorrectedSlice.cc`) makes one slice transition per
  eligible scalar update on a surrogate density and accepts it with the ratio
  of the exact full conditional to the surrogate, so any surrogate that is a
  deterministic function of the graph state leaves the exact target invariant;
  its quality sets the acceptance and the cost only. A dense multivariate
  child's surrogate is the exact numerator minus a rounded-anchor envelope
  normalizer (the `coarse` cache component, unit normalizer where unavailable).
  A cluster or factor child's is the exact numerator minus the log normalizer of
  ONE rule of its own schedule, rule `min(selection.coarse_max_rules, rules)`,
  at the actual parameters (`cpp_selnorm_cluster_step_surrogate_lpdf`,
  `cpp_selnorm_factor_step_surrogate_lpdf`: no change test, certificate, tier,
  sparse grid, QMC, memo or rounding), and its exact `logDensity` where that
  rule is unavailable (best-weight rule, rung beyond the node budget, non-finite
  value). The surrogate path neither reads nor stores the node memo, which holds
  exact values only. A cluster or factor child is coarse-capable when its fixed
  controls allow a rule (vector rule not fixed at a best-weight rule, kernel
  not fixed at normal, rule affordable at the declared rank) and a node is
  eligible when a child is dense or coarse-capable; the others stay exact
  terms. The vector rule and the kernel mode are dynamic parents of these
  children, because a mixture composes both from the bias indicator (a rule
  that is product in every state is still a node); a state with a best-weight
  rule has the exact density as its surrogate, so a model average with only
  best-weight branches engages without any gain and is better run with the
  default sampler. The classes' dynamic (kernel
  mode, vector rule, mean, residual SD, loadings, weights) and control
  validators split `checkParameterValue` for the sampler's fixed-control
  certificate and must stay equal to it. `test-00-selection-coarse-blocks.R`
  pins the surrogate values, fallbacks and purity, and
  `test-00-selection-coarse-exactness.R` the exactness of the corrected chains
  against QUADPACK references.

Known limits: blocks routed to the Smolyak grid (chains of rank three or more
and other supports the certificate rejects) keep the per-transition rule, so
coincident transitions remain unchecked there. The amplitude and its allowance
are measured envelopes, validated by random, adversarial and hill-climbing
searches, not a proof. The qCMDE normalizer grid anchors its error at
`2 * relative_change` (`R/selection-normalizer-grid.R`), which is not an error
bound for the accepted value. The scaled QMC stop is a standard-error
criterion, not a bound: an accepted fallback estimate can exceed
`relative_tolerance` (the `forest_84` calibration accepted an error of 0.75%
at `relative_tolerance = 0.005`). The transition-conditioned proposal addresses the
failure in which every Gaussian proposal misses a material tail contribution;
increasing a multiplier on the resulting near-zero MCSE cannot fix that
failure. Conditioning on one row does not guarantee coverage of every rare
intersection of several rows. Neither agreement between scrambles nor between
nested designs certifies that all important mass was sampled. The construction
is related to rare-event union importance sampling (Owen, Maximov and Chertkov,
2019, doi:10.1214/18-ejs1527), but its union-probability relative-error theorem
does not apply to this product correction. Check empirical uncertainty across
predeclared independent designs, alongside independent reference integrals and
mean-accuracy checks; do not treat three estimated standard errors as a bound
that every individual randomized run must satisfy. In a predeclared check of
200 fresh designs for each of four blocks under adaptive and fixed-maximum
budgets, all eight 99% mean-error intervals were within +/-0.005. One adaptive
tensor block failed the Student-reference tail-coverage screen: 12/200 errors
exceeded three estimated standard errors (94% empirical coverage), and the
2.36 multiplier covered 88%; its empirical variance was 1.04 times the mean
reported variance. Four of its 139 diagnostically accepted results exceeded
0.005 error. This is a documented limitation of the empirical criterion, not
evidence for a confidence guarantee or a reason to change the raw MCSE. On the real workloads
measured (assink2016 fits and post-fit calculations) no randomized fallback is
reached, and the few of the V7-01 fit stop at their first stage under either
multiplier. The certificate costs `1 + 2 rank` evaluations of
each rung it consults and the tier up to 32 evaluations of the finest rung,
each linear in the rows of the block, and the precheck of the finest rule is
gone, so a block that no rung certifies pays the whole ladder before the
sparse grid or the randomized fallback: 0.34 s per call for a narrow-band
block of 20 rows (0.03 s with the precheck) and 0.1 to 0.7 s, fallback
included, for walls of 380 to 2,300 coincident rows. The dense blocks of real
data (5-22 rows) are unaffected.

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
- `src/distributions/selnorm-node-memo.h`: last-evaluation memo of the factor
  and cluster step nodes.
- `tests/testthat/helper-selection-kernel.R`: focused test fixtures.

Test positive and negative effect directions, single and product-space models,
and normal versus active selection rows. Use analytic or independent references
for tail calculations.
