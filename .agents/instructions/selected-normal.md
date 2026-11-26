# Selected-Normal and Posterior Routing

Use this guide when changing posterior evaluation, publication-bias selection
kernels, PIT/CDF/PDF/RNG paths, or effect-direction handling.

JAGS 5's deviance monitor can evaluate an observed likelihood directly. A DIC
penalty that simulates replicate outcomes requires a distribution with a
sampling law. For likelihood-only nodes, verify observed deviance with
`rjags::jags.samples(..., stat = "deviance", summary = "mean")`.

See [maintainer decisions](decision-rationale.md) for why loose requests retain
the deterministic calibration cap and why this does not certify QMC accuracy.

## Dense three-outcome normalizers

For a dense three-outcome product selection block, a failed QMC normalizer
can be recovered with a bounded rank-one integral. This requires positive
off-diagonal covariances and positive weights. At or near a zero uniqueness,
an auxiliary factor can move its loadings inward to resolve the transition.
Its positive unique variances are used only for quadrature. A negative
uniqueness beyond roundoff leaves the original diagnostic in place. The
determinant divided by the sum of the principal two-row minors bounds the
smallest auxiliary covariance eigenvalue from below, including near a zero
uniqueness.
Check the perturbation between the represented factor covariance and the
supplied matrix, then combine its relative error multiplicatively with the
integral enclosure radius. The combined radius must meet the original
tolerance using the supplied work budget. This also applies to the unrestricted
Gaussian event. A truncated event or projected density has a different target.

## Effect Direction

JAGS returns `mu`, `theta`, and random effects on the original fitted scale,
including when `effect_direction = "negative"`. Do not flip those posterior
samples again. `tau` is nonnegative and is never flipped.

PET/PEESE sampling-bias offsets follow effect direction: add them for positive
effects and subtract them for negative effects. Reuse
`.evaluate.brma.bias_offset()`.
Their coefficients are zero in inactive product-space branches, so applying
the row-specific offsets preserves branch routing. A retained logical
`bias_exempt` flag zeroes the offset on that outcome only. Older fits and an
omitted indicator use the original predictors.

Ordinary normal CDF code owns its direction transformation. Selected-normal
code receives the sign through `.selection_context()`. Callers must not add
another flip.

## Product-Space Routing

Predictive variance consumes the kernel's direct centered variance. Raw second
moments remain available to their callers, but are not subtracted from squared
large means. Outer sampling mixtures and scalar moment integrals use a Gaussian
location center and retain the existing numerical error and transition checks.
When the STEP normalizer divided by the largest weight falls below 1e-280,
compute moments on either side of the Gaussian mode in coordinates anchored
at the nearest point of each bin. Combine within-bin and between-bin variance
with relative log masses so that a rare selected tail retains its contribution.

When weighted raw moments overflow, retry the STEP calculation in standard
normal coordinates and scale all bin contributions by the same power of two.
Check the recovered denominator against the log normalizer before using it.

For a continuous outcome query just inside a STEP cutoff, compute the partial
bin mass from its width in outcome units, from the represented endpoint
`z_cut * sei` to the supplied outcome query. A query on a hard support boundary
retains its exact zero tail. Fused summaries and cached plotting CDFs defer
to that calculation when standardized endpoints cannot resolve the width.

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
At a cutoff, product weights compare the signed standardized outcome with
the z cutoff. Best weights compare the signed or absolute outcome with the
cutoff in outcome units. Use the supplied z for a product density queried
on the z scale. In sampling-conditioned likelihoods, a row with zero
integrated diagonal and zero integrated loadings takes its observed value
exactly when that outcome is given. During deletion, a proposed outcome or
sampling error still varies under its conditional law.
For a positive integrated row SD below 1e-12 times the sum of absolute
reconstruction terms, the full-data likelihood refuses locations within
1e-10 times that sum of a changing cutoff in outcome units. Exempt rows,
constant steps and best-rule cuts already won by a fixed publication member
do not trigger the check. JAGS, bridge replay and deletion targets requiring
the full-data post-fit state propagate the calculation error. Fresh deleted
contexts use their own normalizer. Preserve the positive source and its prior
support.
Both rules omit exempt outcomes from the publication event. Product uses only
eligible factors; best takes the smallest eligible p-value in each publication.
An empty event has weight and normalizer one. Preserve the Gaussian density
and all covariance dependencies of the full block. When bounding, predicting
or deleting an exempt outcome, retain selection on correlated eligible
outcomes. Keep physical bin numbers in stored selection contexts and carry
the logical observation field through every subset. Fitted predictions reuse
the flags; explicit prediction rows need a logical `bias_exempt` column when
the fit has exemptions and the prediction uses publication bias.
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

Under product weights the retained rows leave a deleted row's selection event
unchanged, so its LOO-PIT tails and predictive moments are those of the scalar
selected normal at its conditional mean and variance; they are evaluated for
all rows and draws in one batch (and `rstudent()` reduces them with the PSIS
weights natively; like the univariate reducer, it skips draws whose PSIS
weight is exactly zero). Best weights condition the event on the retained rows
and keep the per-row Gaussian event evaluation. The batch uses the fused step
cell summary (paired log tails, mean, variance), which evaluates each
boundary's normal tails and density once and normalizes the tails against the
larger tail mass; a non-step cell, an unavailable mass, a tail mass below
`1e-280`, weights large enough that a bin probability lost to underflow could
reach a tail mass's last digit, or a cell whose shared moments would lose more
than four digits to cancellation (a narrow standardized interval carrying the
mass) falls back to the log-space tail and partition-moment kernels. The gate
only decides where the shared evaluation is exact enough; cells outside it get
those kernels' values unchanged, including their own limits at weights that
span beyond double precision or at subnormal queries.

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
- Loose requests do not relax the calibration of product ladders with
  transitions: the successive-change, tail, shift-amplitude, and refinement-tier
  checks use `min(relative_tolerance, 0.005)` in both the cluster and factor
  kernels. Larger requested tolerances previously accepted coarse rules with
  5% and nearly 100% errors. The default and tighter requests, constant-weight
  blocks, best rules, sparse grids and randomized fallback retain their existing
  policies; in particular QMC still uses the requested tolerance and its
  empirical diagnostic, not a claimed error bound.
- Refinement tier: when the ladder rejects every affordable rung, rank-one
  clusters, tensors of rank at most two and forests of depth at most one try, before Smolyak and
  the randomized fallback (for rank-one product clusters, the bounded product
  fallback below), the average of the finest affordable rung over the
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
  to the Smolyak grid and the randomized fallback, except rank-one product
  clusters, which go to the bounded product fallback below at every tolerance
  (best rules and higher ranks keep the randomized fallback). The tier's sharpness gate
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
- Rank-one product fallback: after the ladder/refinement declines and the
  existing exact event/constant returns, bound each exact row group using
  `P(t) = min(omega) + sum((omega - min(omega)) * bin_probability(t))`.
  Endpoint minima and endpoint/interior-peak maxima enclose each bin over a
  latent interval; finite peaks use standardized width, with ambiguous peak
  membership included. Outward mean/bin arithmetic and an explicit `64*eps`
  log allowance precede intersection with the true global weight range.
  Both Gaussian tails remain in the partition. A fixed-capacity segment tree
  maintains log totals, and a deterministic heap splits the largest absolute
  uncertainty. The actual returned midpoint is accepted only against both
  conservative endpoints at the requested tolerance. Success reports raw MCSE
  zero, relative change equal to the conservative radius, and deterministic
  true. Unresolved arithmetic, zero mass or exhausted work reports NaN log
  normalizer, Inf/Inf diagnostics, deterministic false and density -Inf;
  no QMC rescue or failed-result memo insertion. This is working-precision
  engineering evidence, not universal libm certification.
  Outward rounding of finite `-DBL_MAX` can itself produce `-Inf`; refuse that
  finite-arithmetic failure rather than reclassifying it as an infinite bin or
  zero bound, and preserve this distinction with the private rounding probe.
  Every helper scalar probability call, including preparation and latent
  segment masses, counts against `2*max_points*scrambles*active_groups*excess_bins`.
  Active leaves obey `min(SELNORM_FACTOR_NODE_BUDGET,2*max_points*scrambles)`;
  the tree has `2*next_power_of_two(leaf_cap)` slots and heap at most leaf_cap.
  `max_points` is the kernel argument: R passes `ceiling(max_points_per_scramble/2)`
  (the execution plan's `factor_max_points_per_proposal`), so the user-facing ceiling
  is `max_points_per_scramble` (rounded up to even) `*scrambles*active_groups*excess_bins`
  scalar calls and `min(SELNORM_FACTOR_NODE_BUDGET, max_points_per_scramble*scrambles)`
  leaves, half of what a reading of `max_points` as the R value would give.
  There are no group/bin arrays per leaf. Initial points still initialize QMC
  routes only. Nonfinite transition hints may be omitted while retaining full
  domain coverage; nonfinite actual bound arithmetic refuses.
- Randomized fallback: for higher-rank `product` weights the proposals are the prior, a
  Gaussian with the tilted mean and full covariance of the finest evaluated
  rung (tensor: weighted moments; depth-one forest: from the messages) and a
  copy with standard deviations times 1.5. They require a tensor or depth-one
  forest, at least three finite rungs, a finest rung that passes the
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
- Evaluation memo: JAGS evaluates a node first at its current value, so
  about 17% of the node calls repeat the previous one, and a sampler that
  restores an earlier state (a rejected correction, an enumerated indicator)
  leaves that state one evaluation back. Each factor and cluster step node
  keeps two entries, each with owned copies of every input of a deterministic
  result (`x`, all parameter arrays including the quadrature tables, scalar
  controls, dims and flags) and its log density, and answers a call only when
  every stored value of an entry is equal (pointer identity never suffices:
  JAGS mutates in place); a store replaces the entry used least recently (a
  hit counts as a use). The QMC design is not an input. Only results the
  kernel flags `deterministic` (ladders, tier, rank-one interval refinement,
  sparse grid; the batch results carry the flag as an attribute of that name)
  are stored, never a randomized, non-finite or thrown one. The memo is a
  function-local map keyed by node, never deleted, and cleared as a whole
  above 4096 keys. Each entry copies the quadrature tables (about 32 KB for a
  cluster node, 5 KB for a factor node), so the map can hold roughly 260 MB
  for the rest of the R session, entries of finished models included; a
  process with more step nodes than the cap clears it repeatedly and gains
  little from it. `ROBMA_SELNORM_NODE_MEMO`, read once per process, is `0` to
  disable it and `1` for one entry per node (the earlier memo); the test hook
  `RoBMA_corrected_sampler_test` changes the entry count at run time and reads
  its lookup and hit counters. Draws are bitwise identical with two, one and
  no entries.
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
- Range of the node products: products of many rows leave the double range. The
  nested rule keeps them in `SelNormNodeProduct` (`src/selnorm/selnorm-node-product.h`):
  a long double, bit for bit, where it has more range than double, and a mantissa
  with a binary exponent per node where it does not (Apple silicon). The cluster,
  tensor, sparse-grid and factor-box rules keep long double and have a log-scale
  or randomized fallback, which on such a platform also takes any product that
  reaches the subnormal range at any factor, also when later weights above one
  would bring it back.
- Corrected sampler: `selection.sampler = "slice.coarse"`, the default
  (`src/samplers/CoarseCorrectedSlice.cc`; the surrogate terms are in
  `CoarseSurrogate.h`, the configuration, counters and registry that all
  corrected samplers share in `CorrectedSampler.cc`), makes one slice transition per
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
  rule has the exact density as its surrogate (a model average with only
  best-weight branches falls back by rule (c) below). The classes' dynamic (kernel
  mode, vector rule, mean, residual SD, loadings, weights) and control
  validators split `checkParameterValue` for the sampler's fixed-control
  certificate and must stay equal to it. `test-00-selection-coarse-blocks.R`
  pins the surrogate values, fallbacks and purity, and
  `test-00-selection-coarse-exactness.R` the exactness of the corrected chains
  against QUADPACK references. The scalar factory skips any node with a
  `dbt_block` or `dbt_gate` marker child, whatever the factory order. A dense
  child whose kernel mode is fixed at NORMAL is not coarse-capable (its
  surrogate is its exact density), so `child_kind` decides from fixed parents
  that no surrogate of a sampler can ever differ from the exact density and the
  sampler is not created.
- Corrected block sampler (`src/samplers/CoarseCorrectedBlock.cc`,
  `RoBMA::CoarseCorrectedBlock`, a thin owner of the kernel in
  `CorrectedBlockKernel.cc` that the gated block shares): a declared block
  (`sampler.blocks`) whose
  multilevel view has a coarse-capable selection child. Its factory groups the
  marked members as BayesTools' factor slice factory does (markers read the
  same way) and is asked first: RoBMA's module loads after BayesTools', and
  JAGS asks a later module's factories first; other groups go to
  `BayesTools::FactorSlice`. Gamma-sum coordinates (`GammaGroup.cc`): members
  `eta_i ~ dgamma(a_i, r)` with fixed shapes, one fixed common rate and no
  truncation whose only children (checked on the compiled graph) are sampler
  markers, one aggregate of exactly the group, `sum()` of that aggregate and
  divisions `eta_i / sum` (every other path declines the group) are held as
  the m - 1 log-ratio coordinates `z_i = log(eta_i / eta_m)` and the total
  `T`: `T ~ Gamma(sum a, r)` is independent of the proportions and of the
  likelihood. The Jacobian of `(T, z) -> eta` is `T^(m-1) prod_i p_i`; the
  kernel adds `(m - 1) log T + sum_i log p_i` to the members' own densities.
  The corrected update draws a new `T` from its law after the surrogate sweep
  and before the full evaluation and corrects with `log Gamma(T) - log
  Gamma(T')` (the total's terms cancel in exact arithmetic, so the acceptance
  is that of the proportions, and the last full evaluation is at the state the
  update keeps, which the node memo then holds); the exact kernel draws `T`
  first (Gibbs) and sweeps. A heterogeneity allocation block (allocation SD and
  two Dirichlet gammas) thus has two directions and one gamma draw; the test
  hook's `drop_sum_jacobian` mutant fails the exactness fixtures. One update:
  the full density and the surrogate
  (each member's full density, the exact terms of the other children, the
  selected children's surrogates) at block entry, a fair forward or backward
  order draw, one slice update of the surrogate per direction with chained
  values (BayesTools' header `factor_slice_kernel.h`, `LinkingTo`), one full
  evaluation and the long-double correction; rejection and exceptions restore
  the block-entry values, an exhausted interval the direction-entry values.
  The random order makes the sweep a reversible proposal on the surrogate, so
  one correction is exact; a fixed order with one correction is not (the
  fixed-order mutant of the test hook fails the power fixture). While
  adapting, a rejected correction is followed by one exact sweep (the full
  conditional, the current axes) from the restored entry state, and the axes
  adapt to the state it leaves: started where the surrogate is poor, the chain
  would otherwise reject, not move, and learn collapsed axes that freeze it
  after adaptation too. The adaptive phase need not be invariant (its draws are
  discarded). The proposal still counts as rejected, but the adaptive chain
  leaves regions where the surrogate is poor faster than the corrected kernel
  alone would, so the adaptive rate that the low-acceptance rule (b) reads is
  optimistic; a block kept on that reading is still exact (the rule guards
  efficiency only), and the per-fit warning judges the rate after adaptation.
  The
  fallback rule is the scalar sampler's, judged on acceptance^(1/d) (d the
  number of directions; a block move composes d direction moves), the rate
  exponent that the worst-sampler rate and the warning (`sampling_judged`)
  also use; after a fallback the block runs BayesTools' fixed-order sweep of
  the full conditional with the learned axes, call for call (the test hook's
  `force_fallback` pins this against `BayesTools::FactorSlice`). A member with two finite bounds whose span is not
  finite stops at construction, as in BayesTools' sampler; every error of
  either block sampler gets the hint `sampler.blocks = FALSE`
  (`.sampler_error_hint()`), a member on a bound also the interior start. Test
  hook `RoBMA_corrected_block_test` (fixed-order mutant, evaluation trace,
  widths of the live blocks); exactness fixtures in
  `test-00-selection-corrected-block-exactness.R`.
  Precision limit: declared blocks hold their members on BayesTools'
  transformed scale (distance from a bound, on the log scale) and inherit its
  limit: a posterior narrower than about 1e-13 of a member's distance from its
  bound is not resolved, and this is not detected (see the kernel header). The
  scalar samplers do not have this limit; `sampler.blocks = FALSE` is the
  remedy.
- Corrected gated block (`src/samplers/GatedCorrectedBlock.cc`,
  `RoBMA::GatedCorrectedBlock`, asked before `BayesTools::GatedSlice`): the
  still-unclaimed gated scalar continuous nodes with a fixed support grouped by
  pointer identity of their `<active>` node (with `sampler.nested.weights`,
  the default, the shared coordinates and the fine branch's own coordinate
  form separate blocks, so no live coordinate is redrawn from its prior; fits
  whose weight priors do not qualify keep the ordinary layout). A group is claimed
  when no member depends on another, it has a coarse-capable selection child
  or a singleton step-switch row (the row cache, below) and every member
  reaches every such child only through the child's weights (no other
  parameter of such a child is a member or a deterministic descendant of one):
  selection weights, not a PET or PEESE slope (which keeps
  `GatedSlice`; under the normal kernel of its branch its surrogate is exact).
  Gate closed: every member drawn from its own distribution, redrawn up to 100
  times while a value is not strictly inside its support or a deterministic
  descendant is not finite, otherwise kept (each attempt is an independent
  prior draw, so the kernel is the restricted prior or, with a probability that
  does not depend on the state, no move); while adapting, a change of the
  children's log density is `GatedSlice`'s declaration error, and any error
  of the closed update restores the block-entry values and their
  deterministic descendants before it propagates. Gate open: one
  update of the shared block kernel, whose axes learn from these updates only;
  until a refresh first replaces them (from 100 open adaptive updates; a
  refresh that keeps the axes because an eigenvalue of the covariance is not
  positive keeps these widths too) the coordinate
  axes have widths of 2.5 prior SDs (trigamma of the shapes for gamma and beta
  members and log-ratio coordinates, 1 otherwise). `checkAdaptation()` is true
  after any update, so a branch never open while adapting does not hold up the
  adaptation. Errors start with "Corrected slice (gated block)" and get the
  `selection.sampler = "slice"` hint.
- Screened indicator (`src/samplers/ScreenedFinite.cc`,
  `RoBMA::ScreenedFinite`; the kernel `ScreenedFiniteKernel` takes a setter so
  that a joint owner can call it): a scalar `dcat` (fixed probabilities) or
  `dbern` (fixed probability) without truncation, with a fixed finite support,
  a coarse-capable selection child among the stochastic children of its view
  or a singleton step-switch row whose weights or kernel mode the node reaches
  (the bias label; the row cache, below),
  no `bt_finite_tag`/`bt_finite_witness_tag` child (a `sampler.probabilities`
  readout keeps the enumerating `BayesTools::Finite`, whose recorded
  probabilities are the exact conditionals; a surrogate proposal's
  probabilities are not), and not claimed by RoBMA's joint samplers (asked
  first). One update scores every level of nonzero prior mass with
  `s_l = log prior(l) + sum` of each stochastic child's surrogate (the exact
  density of children without one, the label-dependent prior of shared nested
  weights included), proposes by one of two moves reversible for the
  surrogate law on `exp(s)`, chosen per update by a coin that does not depend
  on the state: with probability 7/8 Liu's Metropolized Gibbs sampler (a
  level other than `k` with probability `p_k' / (1 - p_k)`, accepted on the
  surrogate with `min(1, (1 - p_k) / (1 - p_k'))`; the identity when no other
  level has mass), with probability 1/8 (`plain_proposal_probability`) the
  plain surrogate draw from
  `exp(s)` (current level included), and for `k' != k`
  evaluates the exact scores of `k` and then `k'` (so the last exact
  evaluation is the proposed state) and accepts with
  `min(1, exp((e_k' - s_k') - (e_k - s_k)))`; a self-proposal costs no exact
  evaluation. Each move with its correction is exact, so is the mixture, and
  the plain draw's holding probability makes the chain aperiodic: Liu's step
  alone swaps two equally scored levels at every update, so a chain thinned by
  two kept its initial label (tested with the tie fixture at `thin = 2` from
  both labels; at 1/8 the tied label keeps a holding probability of 1/16 per
  update). The test hook keeps either move alone (`finite_proposal = 2`
  Liu, `0` plain; `3` the default mixture). Liu's step raised the bias label's
  switches on the tutorial and Assink fits by 37% to 2x for 2-10% more time
  than the plain draw. The exact
  kernel is base::Finite's enumeration, arithmetic and single uniform (tested
  draw for draw against base::Finite, and against quadrature with the gated
  block's exact kernel after a forced fallback); its evidence record sums the
  other levels' masses directly (a difference of the cumulative total
  recorded 0 for a leave probability of 1e-30). The
  `finite_proposal = 1` mutant (no correction) fails the exactness fixtures.
  Batch scoring of the six selection labels was not built.
- Row cache of step-switch rows (`src/samplers/StepRowCache.{h,cc}`): the
  surrogate term of the singleton `dselnorm_step_switch` children of the
  corrected gated block and of a screened indicator that reaches their weights
  or kernel mode (`CorrectedBlockOptions::row_cache`, the kernel's
  `row_cache` argument; the declared-block sampler and the scalar corrected
  slice keep these rows as exact terms). In step mode with observed bin
  `o > 0` a row's density is `pw [log phi(sign y; sign m, sd) + log omega_o -
  log A]`, `A = sum_b omega_b p_b`; the fixed parameters' checks are made
  once, both Gaussian terms (normal and step mode) are a function of the
  row's observation, mean and SD and `p_b` of its mean and SD, so the cache
  keeps them for the state it last saw and computes them again when an
  evaluation finds another (the Gaussian terms at once, `p_b` at the first
  step-mode evaluation that needs them; `refresh()` brings the rows to the
  state at the start of each corrected update and of the screened test
  scores), and `value()` reads per row the kernel mode, the weights (their
  validity once per weights vector) and the mean and SD, in the order of
  `cpp_selnorm_kernel_lpdf`. The corrections are exact for any surrogate that
  is one function of the evaluated state, so every cached value and every
  choice between the cache and the node depends on the evaluated state alone:
  a screened level scores the same from every entry level (PET and PEESE
  labels move the means), which the frozen-state test checks bitwise. A row
  goes to its node's own `logDensity` when a parameter other than the mean,
  SD, weights and mode is not fixed or fails the distribution's check, its SE
  is infinite, at its mean and SD a Gaussian term or a bin probability is not
  finite or no bin probability exceeds 1e-280, its cached `A` is not above
  1e-280 (StepRowNormalizer's linear threshold; the log-scale rows), or the
  linear interval masses cannot resolve `A`: the rounding bound `sum_b
  omega_b e_b` (`e_b` from the tails the subtraction cancels and the rounding
  of the standardized bounds) exceeds `1e-10 A` (a bin of large weight whose
  probability is zero or lost its relative precision, such as a central bin
  of width 2.8e-16 with weight 1e60). Telescoping
  changes only the exact normalizer's route. Evidence: a corrected proposal of
  a sampler with cached rows is evidence when `|log acceptance ratio| > 1e-9 *
  max(1, |exact log densities|)` (`step_row_cache_mismatch`) or its exact
  proposal density is zero with a finite surrogate, in addition to the
  coarse children's flag; exact up to rounding, the cache alone neither
  falls back nor enters the record's acceptance rate. Heterogeneity, effect
  and moderator indicators move the rows' means or SDs and keep their
  samplers. Test hook `RoBMA_step_row_cache_test` (mutants: power, observed
  bin, stale cache; a perturbation of the cached probabilities; `disable`,
  the eligibility before the cache; counters; the frozen-state probe of the
  live block kernels); tests in `test-00-selection-step-row-cache.R`,
  exactness fixtures in `test-00-selection-corrected-weights-exactness.R`.
- Location block sampler (`src/samplers/LocationBlock.cc`,
  `RoBMA::LocationBlock`, declared by `R/sampler-joint-location.R` under
  `sampler.joint` "auto" or "location"): one owner of the location gates, their
  normal slabs and the standardized latents of the sampled random-effect term
  of a `.mv` fit. Its factory verifies the declaration on the compiled graph
  (gate charts by deparse, JAGS closed-class linearity of every node that
  reaches a child, members only in the children's means, a probed design with
  affine guards) and declines otherwise. The probe verifies each row of a
  design column against its own response (`column_check`: rounding of the
  row's baseline and perturbed means below 1e-10 of its response, linearity to
  1e-9 of it; a zero entry's rounding below 1e-10 of the column's largest
  response), so one row's entry is never measured against another row's
  larger response (an offset of 32768 lost a row's entries while the column
  passed, corrections reached 141); and at run time the exact correction
  compares, after its write of the final state, every likelihood mean with
  the block's: a mean the sweep moved must equal the one the block tracked,
  and a mean it did not move its entry mean, both to within 1e-6 of the row's
  conditional SD plus 1e-12 of the size of the means and design terms, so that
  an entry the probe missed is never hidden by another coefficient that moves
  the same child. A moved mean's tolerance also holds the rounding of the
  tracked computation: the tracked mean is the entry mean plus each probed
  slope times its coordinate's change over the sweep, and each slope is a
  difference of two probed means at a step of one (the design's and the bias
  slopes' at the factory's baseline, the latent's at the update's entry), so
  its error is at most `4 eps` times the size of those means, times the
  coordinate's absolute change. For the design's and the bias slopes' columns
  the size is that of the probed means themselves, the largest `|mu|` of the
  three probes (`X_size`, `column_size`, as `column_check` bounds a column);
  it does not cover a cancellation inside the graph's expression at the
  factory's baseline. For the latent it is the sum of the absolute terms of
  the entry's means (`_d_size`: the entry mean, the design terms, the
  latent's term and the bias slopes' terms at the entry label), since those
  can cancel. Without it, rows at `3000`, `6000` and `8000 / 1e-6` (which the
  `GEN` measure admits) stopped some seeds in their first adapting sweeps: the
  latent, entering at a state where an intercept of about 6,300 and the
  cluster effect cancelled, jumped by 2e4 to 5e4 with a slope probed to
  4e-13, and the tracked means differed by up to 2.9e-8 against a
  tolerance of 1.6e-8. The comparison runs while adapting too. Either is
  the error `RoBMA_location_untracked_mean` naming
  the child (the entry state restored; mutants 12, 16 and 18 of the test
  hook). An unmoved mean was once required to be bitwise unchanged, but the
  graph's recomputation of it need not be (a Linux build with GCC 10.2.1 and
  OpenBLAS 0.3.29 stopped `metafor-anova-rma-09` that way), so the
  comparison is a diagnostic of the design; the correction is exact whatever
  the size of the difference: every child outside the correction's children
  (the dirty ones, or all after a label change) with a row whose graph mean
  after the final write is not bitwise its entry mean enters the paired
  correction with the graph's own densities at the entry state (the entry
  written back once, those children evaluated, the final state written
  again; only when such a child exists) and at the final state, with the
  precision guard and the endpoint record of a selection child as for the
  others, while adapting too. (Leaving them out, as before, dropped their
  density change `delta z - delta^2 / 2` per row from the ratio: a
  coefficient whose response rounded to zero in the factory's probe, a mean
  `offset + 2^-43 b` probed at an offset of 32768, was sampled without its
  likelihood; mutant 21 leaves them out.) After adaptation the block's
  counters record the unmoved comparisons that were not bitwise equal
  (`unmoved_mean_inexact`), their largest discrepancy in the row's SDs
  (`unmoved_mean_max`; mutant 17 moves the compared entry mean by one ulp,
  mutants 22 and 23 by 0.1 and 10 times the row's tolerance, accepted and
  counted, refused), the children included this way (`unmoved_included`)
  and their largest density change in log units (`unmoved_density_max`); all
  are zero on the Windows build's fixtures.
  Each update rebuilds the Gaussian
  working model from the graph (block precisions from the children's own
  parameters, the latent column probed), runs the component list forward or
  fully reversed (collapsed flips and refreshes along sheared lines, swaps of
  correlated gates, a joint refresh and exact latent draws in Gaussian labels;
  Newton-Laplace tilted flips, refreshes and swaps with a row-marginal screen
  and slice updates of the latents in selection labels, all on the corrected
  sampler's surrogate) and accepts the sweep with one correction against the
  graph's own exact densities of the changed members and children, in every
  label, so the inner representation never changes the target. Fallback to the
  exact likelihood at the first `adaptOff()` with at least 50 updates that are
  evidence about the surrogate and an acceptance below the threshold; it
  registers with the corrected-sampler report and its counters reach the
  joint-sampler record. Test hook `RoBMA_location_block_test` (mutants,
  frozen-state probe, factory decision); exactness fixtures in
  `test-00-location-block-native.R` and `test-00-location-block-exactness.R`.
  Its latents are also the cluster effects `gamma[k]` of clustered `RoBMA()`
  and `BMA()` fits, and "auto" gives it independent fits with a gated
  moderator that no other joint sampler takes (without latents). The Newton
  fits only build proposals: they evaluate each cluster or factor child at
  rule 1 of its schedule (`SelNormCoarseSettings::rule`; back to the sampler's
  rule at the first `adaptOff()` when the stage-2 acceptance of the tilted
  flips and swaps was below 0.9 over at least 50 of them while adapting), and
  a line whose selection children are all step-switch rows uses the rows'
  closed-form normalizer derivatives (`StepRowNormalizer`) instead of central
  differences; Metropolis-Hastings ratios and the exact correction keep the
  sampler's rule, and the step-switch rows' densities are the distribution's
  own evaluations with the parameters read once per update. An elliptical
  slice refresh of all active coefficients (component E, latents sheared) is
  exact but stays behind the test switch `RoBMA_location_block_test_config`:
  measured against the per-group refreshes it lost effective samples of
  individual coefficients on the 18-moderator tutorial model. The fit's
  sampler record carries an evidence record per gate and chain (the updates
  and the gate's transitions by direction between each update's entry and the
  state it retains, so a rejected sweep counts nothing and a swap counts for
  both gates; and, named `inner_*`, the inner stage's flips by direction,
  screened proposals and log odds of the other state) and, for the
  screened indicators, per indicator and chain (`sampler_declarations$robma$
  evidence`), for consumers that tell decided from stuck indicators.
  With a publication-bias mixture of none, step-selection, PET and PEESE
  branches the block also owns the bias indicator and the PET and PEESE slope
  sources (`src/samplers/LocationBlockBias.h` for the declaration and
  verification, `LocationBlockBiasMoves.h` for the moves). R declares them with
  one tag (`robma_location_bias`: indicator, family codes, sources, intercept
  slot) unless `sampler.probabilities` names the indicator; the factory
  verifies the indicator's fixed categorical law, each source's scalar prior
  with fixed parameters and support exactly `[0, Inf)` entering only through
  `source * equals(indicator, label)`, the probed slope columns of the means
  (affine, additive, without effect at other labels), the kernel mode (NORMAL
  exactly on the labels none, PET and PEESE), and that the label reaches the
  likelihood otherwise only through parameters the Gaussian working model does
  not read; gate markers and the nested weight node (`BayesTools_nested_a`,
  whose density enters every label score) are the only other stochastic
  descendants. A failed check declines the bias members and keeps the
  location members (the factory decision names the reason); the selection
  weights keep `GatedCorrectedBlock` or `GatedSlice`. A label is evaluated by
  writing the indicator and its label nodes (weights, kernel modes, rules,
  the nested node's parameters) without the means: every child is scored at
  the block's working means. Sweep conditions: the slot list is fixed for the
  fit (label-dependent components decide when they run whether they act), the
  surrogate inner target is a property of the fit, composite components
  reverse their order in the reversed sweep, the correction writes the
  indicator back to its entry value before reading the entry's exact
  densities and counts every child and the nested node as changed when the
  label changed, and a sweep counts as evidence about the surrogate only when
  an endpoint is a selection label. Components, first in the forward sweep and
  mirrored in the reversed one: D draws dormant slopes from their priors; L0
  is Liu's Metropolized Gibbs step on the label's surrogate full conditional,
  or after a state-independent coin with probability 1/8 (`label_coin`) the
  plain Gibbs draw over the label scores (the current label included), with
  the other labels' masses summed directly at both endpoints. No move of the
  block is a periodic two-state move: a Metropolis step between two states
  with ratio exactly one alternates deterministically, so an evenly thinned
  chain could keep one state. The first stage of the gate flips and the
  location swaps therefore uses Barker's acceptance `r / (1 + r)` instead of
  `min(1, r)` after a coin with probability 1/4 per attempt
  (`two_state_coin`); the later stages are unchanged, since the ratio of the
  two directions' first-stage acceptances is `r` under both (test switch
  `coin` of `RoBMA_location_block_test_config` forces either way, and L0's
  coin too). L2 uses the plain Metropolis first stage: it is attempted on a
  state-independent quarter of the sweeps, so it cannot alternate
  deterministically, and its recorded pass probabilities are those of the
  kernel that runs;
  L1 (Gaussian labels) proposes label, intercept gate and slab, both slopes
  and the latents from the collapsed Gaussian cells (an independence
  Metropolis-Hastings step computed from the inner target and both proposal
  densities; prior cells propose the slope from its prior); S (PET and PEESE)
  slice-samples the active slope on its collapsed conditional, draws the
  intercept and latents from their conditionals and accepts with the ratio of
  the inner target to that law; L2 swaps between the families with label,
  gate, intercept and slope moved together along the shears `K_t` (the
  intercept's line) and `K_b` (the slope's conditional-mean shear), holding
  `E_0 = E - K_t t - K_b theta`: label proposals by the selection labels'
  row-marginal screens or the Gaussian labels' evidences (5% uniform),
  absolute evidences at the base point (`t = 0`, no slope, `E_0`: the
  Gaussian plane in closed form, a selection label's value there times its
  Laplace odds along the line), stage 1 on the evidences and the label
  Hastings factor, stage 2 on the inner target with the complete proposal
  densities. L2 is the identity when its base point is not representable:
  each label's evaluation of the children there (G0's and a fitted selection
  label's absolute evidence) sums the absolute values of their terms, and
  when the largest sum `M_0` exceeds `likelihood_representable_limit` (1e9)
  the swap makes no proposal in that sweep (nothing further is drawn, no
  proposal, attempt or inner-target change is recorded, and
  `swap_unrepresentable` counts the sweep in the block's counters and, after
  adaptation, in `evidence$bias`), while adapting too, without an error.
  `M_0` is a function of what the move holds fixed (`E_0` and the unchanged
  state, evaluated equal up to rounding from both ends; the probe records
  both ends' `M_0`), so the pair is skipped from either end and the kernel
  stays exact; entering a selection label, its terms join `M_0` after the
  target label's draw, and the reverse move checks them before its own. The
  two directions evaluate their evidences and laws at base points rebuilt
  from their own ends, equal in exact arithmetic; the damped Newton fits
  reject a step only when it lowers their objective by more than
  `1e-10 (1 + |f|)` (`fit_tolerance`), so on ordinary data both ends' fits
  take the same steps and the reverse ratios are the inverses to about
  1e-13 (accepting on `f(next) >= f0` let rounding decide the last step and
  left residuals up to 8.1e-8). At extreme rows the fits amplify the
  rounding of the working means instead: with one row at `1000 / 1e-6`
  (`M_0` 1.2e7, the Laplace means up to 7e-4 SD apart) the two ends' log
  ratios differed by up to 6.5e-4. So after stage 1 and the move, L2
  rebuilds its proposal end as an attempt that starts there computes it
  (`swap_back`, the same `swap_base` with the pair's selection fit, label
  proposals, stage 1 and proposal densities, from the working state alone:
  `_h`, the line, the fits and every cached child value rebuilt), and the
  move is the identity (`swap_irreproducible`, in the block's counters and,
  after adaptation, `evidence$bias`; the attempt's stage-1 records stay and
  nothing else is recorded) when the geometry or a fit fails there, the base
  point there is not representable, `|stage 1 + reverse stage 1|` exceeds
  `swap_reciprocity_tolerance` (1e-9, about 1e4 times the residual of
  ordinary data; a missing Hastings term is at least 1e-3), or so does the
  sum of the two ends' log ratios of the proposal densities,
  `log q(x|x') - log q(x'|x)` at each end's base point (a guard, not a
  requirement of exactness: with stage 2's exact ratio below, only stage 1's
  reciprocity enters detailed balance). Given the working means, both ends
  evaluate the same two base points, so the rule decides the same from
  either end of the evaluated pair; after adaptation the counters keep the largest residuals of the attempts
  that went on to stage 2 (`max_swap_stage1_residual`,
  `max_swap_full_residual`).
  The rebuilt end shares the forward end's working means (moved by the
  swap's shift), while the chain's own reverse attempt, starting at the
  proposal in a later sweep, reads the graph's means there and moves them by
  that sweep's shifts; the two differ by rounding, which at extreme rows the
  fits amplify (at `1000 / 1e-6` rule residuals up to 4.7e-4 with
  bitwise-equal base points). So the base point is also not representable
  for the swap when the working means' representability criterion
  `Q = sum_i delta_i (1 + |z_i|)` exceeds `swap_working_tolerance` (1e-8),
  over every child's rows, with `delta_i = 4 eps size_i / sd_i` (`size_i` the
  sum of the absolute terms of the row's mean at the end's state before its
  move to the base point: the working mean, the design terms with the
  intercept's, the latent's term and the active slope's term;
  `sd_i = W_ii^-1/2`; `eps = 2^-52`; 4 eps the rounding the tracked
  comparison allows a probed slope) and `z_i = sd_i (W r)_i` the
  standardized working residual at the base point
  (`swap_working_criterion`). The working means were read and carried where
  the intercept and the slopes hold their values; the move to the base point
  removes those terms from the means but not from their rounding, so `size_i`
  is taken before it (taken at the base point, it gave data shifted far
  from zero the `Q` of ordinary data, 1.5e-12 at `|y| / SE` of 1,400, while
  re-read means moved the rule's inputs by up to 3.3e-6 and 85 to 97% of the
  attempts that passed stage 1 were the identity after it). `Q` is thus a
  function of each end, not of the base point alone: every attempt
  evaluates it at its own end wherever `M_0` is (`swap_unrepresentable`) and
  at its rebuilt end (`swap_irreproducible`), and the reverse attempt at the
  same two ends; after adaptation the counters keep the largest `Q` over the
  attempts' own base points that pass the check (also in attempts that end
  before a proposal) and over the proposal ends that were rebuilt
  successfully (`max_swap_working`). `Q` is a criterion, not a bound. It
  adds one term per row, so it grows linearly with the number of rows; a
  row's term grows with the row's size in units of its conditional SD,
  `size_i / sd_i` (about `|y_i| / sd_i` for an outcome whose mean the
  intercept carries), and with its standardized residual at the base point,
  `1 + |z_i|` (the base point holds the intercept at zero, so for such an
  outcome `|z_i|` is of the same order); for rows alike `Q` scales as
  `n (|y| / SE)^2`. Measured as the largest `Q` after adaptation of one
  chain: 2.3e-11 on the real cases (`gender_3lvl_dec_def_long`,
  `assink_V_mods`, `che_int`); 9.1e-11 on a 10,000-row clustered SMD fixture
  (`|y| / SE` 0 to 12, median 2); on clustered SMD fits with a two-level
  moderator, 2.1e-10 with 1,000 rows at a median `|y| / SE` of 10, 6.9e-9
  with 1,000 rows at 61 and 8.0e-9 with 10,000 rows at 20; 4.8e-10 on 1,000
  clustered Fisher-z rows at a median of 43 and 3.9e-9 on 10,000 such rows;
  on the 10-row bias fixture with one row at `2 / 1e-3` up to 4.7e-11 and
  at `20 / 1e-4` up to 4.1e-9 (chain seeds 1-40), and at least 4e-6 at
  `1000 / 1e-6` and `8000 / 1e-6`. Outcomes located far from zero relative
  to their standard errors are skipped by `Q` (`swap_unrepresentable`), and
  the chain keeps its other label moves; the skip starts at these medians of
  `|y| / SE`: on the 10-row bias fixture with every estimate and the
  intercept's prior mean shifted and the SEs scaled by 10 or 100, every
  attempt from 980 (none at 930), and with the fixture's SEs of about 0.2,
  where its cluster latents take up part of the location at the base point,
  some attempts from 1,210 and every one from 1,860 (none at 1,170); on the
  clustered SMD fits (L2 in every sweep), 196 to 198 of 200 sweeps,
  adaptation included, with 1,000 rows from 81 (none at 61; at 121 the
  skipped `Q` was 2.4e-8 to 2.7e-8) and with 10,000 rows at 30 (none at
  20). Just below, the rule's residuals reach its tolerance instead on the
  10-row fixture: at a median of 700 to 930 (SEs
  times 10) 64 to 73% of the attempts that pass stage 1 are the identity
  (`swap_irreproducible`), at 470 to 560 about 30%. `Q`'s two passes over
  the rows (the sizes and the criterion) cost less than the run-to-run
  variation of the timing: L2 took 0.116 to 0.117 ms per attempt with and
  without them on the bias fixture (four alternating runs of 10,000
  sweeps), and 10.5 to 10.7 ms on 1,000 Fisher-z rows (two alternating runs
  each). The probe rebuilds each
  dry run's proposal end from working means read from the graph at the
  proposal, with the forward sweep's latent slopes (kind 7; a sweep that
  starts at the proposal probes the slopes again, which the probe does not
  reproduce), and records the change of the rule's inputs against the
  rebuild from the carried means. Over chain seeds 1-40 (L2 in every sweep,
  about 32,000 dry runs per fixture), in the dry runs the whole rule admits,
  that change is up to 8e-12, a floor of rounding, and 58 `Q` on the default
  fixture (99.9% within 5.8 `Q`), up to 3.4e-10 and 11 `Q` at `2 / 1e-3`, and
  up to 3.3e-8 and 19 `Q` at `20 / 1e-4`; on relocated fixtures below the
  skip (seed 4127) up to 2.5e-8 and 4.5 `Q`. The selection labels' Laplace
  fits can move the rule's inputs by up to about 1e-5 in rare admitted
  attempts, at initial and at ordinary states, whatever `Q` is: the fits run
  two damped Newton iterations whose step halvings are decided by
  comparisons of the objective, and the fresh rebuild's fit, with `G0`
  equal, can stop at another mean. With one row at `1 / 1e-3` (chain seeds
  1-40, 100 draws, about 12,000 rule-admitted dry runs), 3 changed by more
  than 1e-7: by 6.9e-6 at seed 17's first dry run, at the chain's initial
  state (`M_0` 2.5e6, `Q` 7.5e-9, the Laplace means 8.8e-6 law SDs apart),
  by 3.4e-6 at seed 11 (`Q` 7.6e-9), and by 2.2e-7 at an ordinary state
  (seed 1, the 75th dry run: `M_0` 7.4, `Q` 7e-12, the means 4.6e-7 law SDs
  apart); at `2 / 1.6e-3` by up to 9.5e-7 (seed 12, `Q` 6.1e-9). At
  `2 / 1e-3` `Q` skips the initial states' cases: with `size_i` taken at
  the base point, 7 admitted dry runs there changed by up to 1.6e-5, each
  among its chain's first three dry runs, where `Q` from the end's sizes is
  2.4e-8 to 3.1e-8.
  Stage 2 takes each proposal density from the base point of the end that
  draws it, `r = delta + log q'(x|x') - log q(x'|x)`, the exact
  Metropolis-Hastings ratio of the two evaluated ends and the exact negative
  of the reverse move's evaluated at the same two ends
  (a symmetrized ratio of the two ends' own-base ratios, `(full - full') / 2`,
  is off by half the sum of the two ends' discrepancies of the proposal
  densities, up to 3.4e-4 on that fixture, which the rule bounds only when
  one end's discrepancy dominates). The swap is exact for the inner target
  given the working means, up to `swap_reciprocity_tolerance` in stage 1,
  for the two ends it evaluates; the working means are carried through an
  update, and the swap runs only where `Q <= 1e-8` at both ends and elsewhere
  is the identity. The chain's own reverse attempt starts from working means
  re-read at the proposal and from latent slopes probed again there; what is
  measured is a rebuild with re-read means and the forward sweep's latent
  slopes (kind 7), which agrees with the rebuilt end to within 58 `Q` and
  3.3e-8 in a log ratio on the default, `2 / 1e-3` and `20 / 1e-4` fixtures
  and by up to about 1e-5 in the rare admitted attempts above. That
  agreement is measured, not guaranteed, and the effect of the re-probed
  slopes is not measured.
  The rebuild adds one base-point evaluation with one fit to every
  attempt that passes stage 1. On ordinary fixtures nothing is skipped; at
  `1000 / 1e-6` every attempt after adaptation is skipped at its base point
  by `Q` (without the `Q` rule most are the identity after stage 1, and the
  swaps that reach stage 2 have residuals below 2e-13), at `20 / 1e-4` a few
  are the identity after stage 1. The probe records the rule's residuals of
  each dry run, `M_0` and `Q` of each dry run skipped at its base point, and
  the change of the rule's inputs when the proposal end is rebuilt from
  working means read from the graph, after each accepted
  live swap an attempt started at its state, forced back, whose end must
  equal the rebuilt one bitwise (`Q` included: both builds of that end take
  their sizes from the same working state), and for each live attempt that reaches
  stage 2 the acceptance's input with stage 1, the exact ratio and the
  forward end's own ratio from the attempt's quantities (kind 10) (mutant 19
  scores the reverse with the forward end's base point, as the swap did
  before, and skips nothing after stage 1; mutant 20 runs Barker's first
  stage; mutant 24 leaves `Q` out of the rule; mutant 25 takes an entered
  selection label's base point as not representable after the target label's
  draw; mutant 26 takes stage 2's ratio from the forward end's base point
  alone; mutants 27 and 28 run only the dormant draws and L2, rejecting every
  proposal at stage 2, 28 with the snapshot leaving the evidence flag as the
  swap's evaluations set it; mutant 29 takes the rebuilt end's `Q` as 10
  times its tolerance, so that on ordinary fixtures, where the own end, `M_0`
  and both residuals pass, only that check refuses, and mutant 30 is mutant 29
  without it; mutant 31 takes the rebuilt end's sizes from the attempt's own
  end). The gate of the target cell is drawn and
  scored with the same two log probabilities, `-log(1 + e^-d)` and
  `-log(1 + e^d)` of the difference `d` of the cell's two evidences
  (`swap_gate_log_probabilities`): normalized by their total, evidences of
  size 1e9 rounded the smaller cell's draw probability to zero while its
  score kept 4.1e-8. With JAGS's Mersenne-Twister generator, whose uniforms
  lie on the grid `k 2^-32` (JAGS's other generators, Super-Duper and
  RngStream, which the fits' initial values select, are not on that grid),
  the draw `u < p` realizes any probability to within an absolute `2^-32`
  (2.3e-10), as any such draw does; the relative error is large only for
  very small cells: at that pair (`d` = 17) the off cell's 4.14e-8 is drawn
  with probability 4.12e-8 (relative log error -4.6e-3), at `d` = 20 the
  relative log error is -0.10, and a cell below `2^-32` (`d` above 22.2) is
  never drawn while its score stays `e^-d`; large relative errors concern
  only events below about 1e-7 per attempt.
  The swap's evaluations at its base point, at its rebuilt end and at a
  proposal it discards leave the sweep's evidence flag about the surrogate
  (`_not_exact`) as it was (the working state's snapshot restores it; the
  exact correction's ratio depends on the sweep's endpoints only); an
  accepted proposal's evaluations keep it. Only the swap resets the flag at
  a discarded proposal: the other moves' discarded proposals still mark the
  sweep. The flag enters no acceptance ratio: it decides which sweeps count
  as evidence about the surrogate, in the counters and, while adapting, in
  the acceptance rate that can switch the sampler to its exact kernel at the
  end of adaptation (`selection.coarse_min_acceptance`). The tilted fits of the location components (flips, refreshes and
  swaps, whose reverse fit is evaluated at inputs that differ from the
  forward ones by rounding) take their steps by the same rule, so they too
  agree between the two ends up to rounding. L2 runs with probability 1/4 per sweep (`default_swap_probability`;
  with its Laplace fits it took 15.7% of the sampling sweep on the Assink
  full-V moderator fit and 10.6% on the 18-moderator tutorial fit at one, 8.5%
  and 5.6% at 1/2, above the 5% of the rule, and 4.4% and 2.9% at 1/4), decided
  before any component and independent of the state. The test hook's
  probe records, at each component's slot, L1's analytic against computed
  ratio, S's acceptance ratio, L2's reciprocity, rebuilt base quantities and
  probability flows, the fitted plane against the children's densities, and
  L0's buffered scores against the graph's; mutants 6-11 remove parts of L2
  or take the joint slot from the entry label (a target term left out of
  both a component's ratio and its recorded inner-target change is restored
  by the exact correction, so the mutants remove terms from the ratio only);
  the trace switch records every
  L2 attempt (`RoBMA_location_block_test_trace`). The evidence record adds
  `evidence$bias` (per chain the proposals and acceptances of L0, L1, S and
  L2 and the sweeps L2 skipped, `swap_unrepresentable` and
  `swap_irreproducible`; by source family,
  `leave_selection_*` and `leave_gaussian_*`, for L2 the count, sum and maximum of the stage-1 log ratio and the sum of
  `min(1, exp(stage 1))` of all attempts and the count, sum and maximum of
  the stage-1 log ratio of the swaps the exact correction kept, and for L0,
  from the label scores it computes and without another likelihood
  evaluation, its updates, the sum of their probabilities of moving to the
  other family on the inner target of the kernel that ran (after the coin
  Liu's `p_j / max(1 - p_i, 1 - p_j)` or the Gibbs draw's `p_j`; the screened
  indicator's record has the same convention), the
  moves it made there and the maximum of `log P(other family) - log P(own
  family)`), `evidence$bias_labels` (sweeps ending in each label) and
  `evidence$bias_sojourns` (per chain, from all retained draws of a
  monitored indicator: the first draw's family, the initial run's length and
  whether it is censored, the crossings in each direction and the longest
  run in each family). These records are descriptive (no status, warning or
  default depends on them): a family that is never left although swaps out
  of it pass stage 1 points to a stuck chain; very negative ratios and
  near-zero L0 probability sums throughout are consistent with, but do not
  prove, decisive evidence; the sums are conditional on each chain's state,
  so a censored initial run in some chains only, while others visit the
  other family, points to chains that do not mix. Exactness fixtures (labels none, selection, PET,
  PEESE; nested; a second selection label) in
  `test-00-location-block-exactness.R`; frozen-state identities and the
  sweep's discrete-chain check in `test-00-location-bias-native.R`.
- Heterogeneity scale of the location block (`src/samplers/LocationScale.cc`,
  declared by `R/sampler-joint-location-scale.R` with the location block when
  it has latents): the block also owns the heterogeneity sources of the
  sampled term (`.mv` allocation layout: the allocation SD and the gamma
  sources of the component's Dirichlet split, or the SD alone for one term;
  clustered fits: the stochastic sources of `tau` and `rho` behind their
  gates). With `v_j = a^2 w_j`, the scale group `T_c` multiplies the sampled
  component's variance by `c^2`, keeps the others and divides the latents by
  `c` (gamma sources keep their total), so no child density changes; one
  slice update of `u = log c` from 0 (width 0.5, a total budget of 40
  expansions split at random between the sides, Neal 2003; a fixed cap per
  side is not a valid kernel) of `log f(u) = log p_a(a') + log p_w(w') -
  (2L - 1) log a' - sum(Z^2) / (2 c^2) + (2 - K) u` from the sources' own
  `PDF_PRIOR` densities (Liu and Sabatti 2000) needs no likelihood
  evaluation. One update: the location sweep, the sources' ordinary samplers
  (built by asking the factories after the location block's on exactly these
  nodes and embedded unchanged, so their adaptation, fallback and
  corrected-sampler record are those of an ordinary fit), then the move, which
  is the identity unless every declared gate indicator is at its declared
  level. The move acts on a set defined by the state alone: finite sources
  strictly inside their supports and the share (clustered) or every
  allocation weight (gamma) in `[2^-20, 1 - 2^-20]`, within which the
  variances the graph rebuilds from the rounded weights keep a relative
  accuracy of about `2^-33`; from a state outside it the move is the identity
  (counted as `outside`) and a transformed point outside it has zero density
  on the orbit, so the move is a kernel for the posterior restricted to the
  set. `d = (1 - w_s) + w_s c^2` is formed from the sources as a sum of
  positive parts. After each move's write every parent value of the
  stochastic children of the sources and latents (markers excluded) that the
  write can change is compared with its value before (`1e-6 (1 + |value|)`);
  a change restores the entry and raises the error `RoBMA_location_scale_
  invariance` naming the child (measured cost: 0.12 ms per sweep, 0.3%, on
  `gender`; 0.004 ms, 0.07%, on `assink_V_mods`). A relative tolerance on
  parent values has no unit-free meaning (in units of 1e-5 a mean term
  `0.5 sigma` changed the rows' log density by -0.63 within it), so after
  adaptation every move also evaluates the children's log densities before
  and after the write, each paired with its own value (a long double sum).
  A change of their sum above `max(1e-6, 1e-9 M)` (`M` the sum of their
  absolute values before the write, which carries the precision guard below)
  restores the entry and raises the same error naming the child whose log
  density changed most: a design diagnostic, since the orbit leaves the
  likelihood unchanged in exact arithmetic (checked on every move: when only
  the first and every 64th were, a hand-tagged graph passed a small first
  move and committed +1.40 on the second). Below it the move is a proposal
  that the children's likelihood corrects: the write is kept with
  probability `min(1, exp(change))`, without a draw when the change is not
  negative and with one uniform otherwise (the entry restored on rejection).
  The slice update is reversible for `f` and the paired changes of a move and
  of its reverse are exact negatives, so the corrected kernel is exact for
  the model as the graph computes it: a valid clustered graph inside the chart
  at `M = 1.94e8` changes by rounding up to 0.006 per move (within its
  tolerance 0.19), and a hand-tagged graph whose change stays below the
  tolerance (a child of `sigma` alone spanning 0.85 in log density, a
  constant child far from its mean setting `M` near 9e8) is sampled as the
  stock samplers sample it, while without the correction (scale mutant 5)
  it is not. Counters: `likelihood_checks` (the moves evaluated),
  `likelihood_changed` (a nonzero change), `likelihood_rejected`,
  `likelihood_change_max` (the largest absolute change, log units) and, with
  the timing switch, `seconds_likelihood`. While adapting only the parent
  values are compared (no evaluation, no correction). The factory declines
  (everything as before; draws identical) unless
  the sources are free, unobserved, fixed-parameter scalars with the declared
  supports (a gamma group recognized by `find_gamma_groups` behind the
  declared ratios), the block owns all `K` latents, no source has a gate
  marker or shares a block marker with another node, the sources' samplers
  claim exactly them, and the chart check passes: once at construction
  (chain 0, state restored) at values the declaration alone fixes (the
  sources at interior points of their supports, the latents at nonzero
  values, the indicators at their declared levels, every bias label), the
  charts of `a`, `w` and `sigma_s` match and `sigma_s' = c sigma_s` for `c`
  in {0.7, 1.3} (moved halfway towards one until the point stays in the
  move's set), so that ownership does not depend on initial values or differ
  between chains and workers; the likelihood's invariance is the per-move
  check's (a dependence it finds is an error, not a decline). The tag passes the
  gamma sources through their ratio nodes so that it adds no child to the
  gamma group. The record `sampler_declarations$robma$joint$scale` holds the
  factory decisions (layout, L, K, status, reason, invariance check, embedded
  samplers) and the move's counters. Test hooks
  `RoBMA_location_scale_test_config` (mutants, decline, timing) and
  `RoBMA_location_scale_test`; tests `test-00-location-scale-native.R`, the
  quadrature fixtures in `test-00-location-block-exactness.R`.
- Fallback and record: each sampler object (node, block or indicator x chain,
  per model build) keeps its own adaptation flag (JAGS has no accessor) and
  counts its adaptive evidence proposals and their accepted corrections
  (`CorrectedTally`). Evidence is a proposal during which at least one
  surrogate used a normalizer that is not the exact one (the trailing
  `surrogate_not_exact` flag of the three `surrogateLogDensity`: a cluster or
  factor child's coarse rung, never its `exact_fallback`; a dense child's
  anchor integral or unit normalizer; no child under the NORMAL kernel, whose
  anchor normalizer is the exact one); for an indicator, a change of level
  whose surrogate at either level was not exact; for cached step-switch rows,
  a mismatch beyond rounding (the row cache, above). A Gaussian label's proposals
  are therefore not evidence (their correction is one). At its first
  `adaptOff()` (a second call changes nothing) and with
  `selection.coarse_min_acceptance > 0` (default 0.90, frozen with the
  runtime settings; settings without it, from 4.1.1 fits, mean 0) it switches
  for good to the exact kernel, an ordinary slice update of the full
  conditional through the same `logDensity()`, only when it made at least 50
  adaptive evidence proposals and their acceptance (rate^(1/d) for a block) is
  below the threshold; with fewer it stays corrected, and the per-fit warning
  over evidence proposals after adaptation is the safety net. (Earlier
  development versions also fell back after fewer than 50 adaptive proposals or
  without an inexact surrogate, which disabled every corrected sampler of a
  chain that adapted in a Gaussian label; those two reasons are retired and
  their record elements dropped.) A sampler whose children can never use an
  inexact surrogate is not created (`child_kind`, fixed parents), except for
  cached step-switch rows, whose gain is the cheaper evaluation. Both kernels
  are exact and the choice uses adaptation draws only, so no switch happens
  after adaptation. Adaptation ends where BayesTools
  calls `rjags::adapt(end.adaptation = TRUE)` during fitting or reconstruction,
  or where JAGS ends it on setting a monitor (`coda.samples()`); rjags'
  `update()` alone does not, so a hand model updated without either keeps its
  samplers adapting (and counting as adaptive). Per-fit counters
  (`RoBMA_selnorm_sampler_report`, read and cleared) count objects, fallbacks
  (`fallback_low_acceptance`), post-adaptation evidence proposals and
  acceptances of samplers still corrected (burn-in and all sampling
  iterations, including thinned iterations) and the lowest rate
  of one sampler over its evidence proposals since the last clearing read (a registry
  of live samplers under a mutex; a sampler destroyed earlier leaves its rate
  in a retired minimum; samplers with fewer than 1,000 such proposals do not
  count); the three
  process-wide counters of `selection_sampler_info()` count every corrected
  proposal (evidence or not, adaptation included; an indicator's changes of
  level), and the test hook `RoBMA_corrected_sampler_test` counts updates,
  proposals, evidence, fallbacks and exact scores per sampler kind. `.fit()`
  passes a report callback (`BayesTools::JAGS_fit(runtime_report = )`, also on
  extensions) for selection fits, sums the retained runs over processes, adds
  the input fit's record once and stores `attr(fit, "selection_sampler")`
  (`complete = FALSE` when part of the chain was drawn without a record: an
  extension of a fit saved before the record existed, a process whose report
  failed, or draws added without any report); the counters
  are discarded at the `start` phase in `.selection_runtime_configure()`, never
  in the stored setup closure. A record with 1,000 or more counted proposals
  and a pooled or worst-sampler rate below the threshold warns
  (`RoBMA_selection_sampler_low_acceptance`, family
  `RoBMA_selection_sampler_warning`) and keeps one message in
  `attr(fit, "warnings")`. Development versions also recorded
  `allocation_block_withheld` and warned on fallbacks of ensembles without
  their allocation block; the accessor drops that element and a refresh of the
  warnings drops that stored message (`.selection_sampler_fallback_pattern`).
- Precision guard (`CorrectedSampler.{h,cc}`, one constant
  `likelihood_representable_limit = 1e9`): data that the input checks admit
  can hold an estimate that no state of the model fits within tens of
  thousands of standard deviations (`yi = 1000`, `sei = 1e-6` with the effect
  fixed at zero); its log density (1e10 to 1e26) leaves double precision unable
  to resolve the differences the samplers need, and the reviewed fixtures drew
  the prior from these and from the stock JAGS and BayesTools samplers alike.
  Wherever one of this package's samplers evaluates the exact log density of
  the state its chain is in (the location block's exact correction over the
  entry's changed children, the block kernels' entry density in both kernels,
  the scalar corrected slice's entry density and the exact slicer's first
  evaluation, the screened indicator's exact score of the current level and
  its enumeration's, the scale move's evaluation of the children) it sums
  the absolute values of the
  likelihood terms of that evaluation (`M`, the children's log densities,
  priors excluded): the terms the move differences, at the chain's state,
  not every value the model can give (a chain can sit in such a state while
  other chains sample the same data). After adaptation `M > 1e9` raises
  `[RoBMA_likelihood_not_representable]` naming the node with the largest term
  and its value, and `.sampler_error_hint()` replaces the sampler hint with one
  about the data and the chain's state (the ordinary samplers return draws for
  such data without a warning, so no sampler option is a remedy; other initial
  values or a longer adaptation can avoid a state that other chains avoid).
  While adapting nothing is raised
  (a chain may start far from the data), and a proposed state is never checked
  (a huge term rejects it as before). Below the limit a sum of `n` terms
  accumulated in double precision loses up to about `(n - 1) u M` (`u =
  2^-53`; measured up to 9e-6 at `M = 1e9`, `n = 10,000`): the error of one
  total. A difference of two totals, as the bias-side samplers keep, carries
  the sum of both totals' allowances plus the rounding of the subtraction
  (the two totals can round in opposite directions: with two terms per total
  near `M = 8e8`, 1.17e-7 against the one-total 8.9e-8). The guard evaluates
  the terms of the current state, not those of the proposed total. At the R
  boundary `.stop_fit_errors()` signals a fitting error whose
  message carries a bracketed name of RoBMA's samplers with that name as its
  condition class and the parent class `RoBMA_sampler_error`
  (`RoBMA_likelihood_not_representable`, `RoBMA_location_untracked_mean`,
  `RoBMA_location_scale_invariance`; the message is kept and the classes are
  listed in `?selection_sampler_info`). The guard sees only the state the chain is in; L2's evidences are
  evaluated at its base point, so L2 applies the same limit there and makes
  no proposal above it (`swap_unrepresentable`, the L2 paragraph above), a
  rule of what the move holds fixed that keeps the kernel exact. L2's
  frozen-state fixture with one row at `y = 1000`, `se = 1e-6`, which its
  cluster latent absorbs, has base-point sums of about 5e17 in the first
  adapting sweeps (the relative reciprocity residual 1.17 before the rule)
  and is skipped there; after adaptation the state (`M` = 1.3e7) and the base
  point (`M_0` = 1.2e7) are representable, but its working means are not for
  the swap (`Q` at least 4e-6, the tolerance 1e-8), and every attempt is
  skipped at its base point (`swap_unrepresentable`, the L2 paragraph above).
  Without the `Q` rule (mutant 24) swaps run there, and the attempts whose
  proposal end, rebuilt with Laplace laws up to 7e-4 SD from the forward
  end's, does not reproduce them to `swap_reciprocity_tolerance` are the
  identity (`swap_irreproducible`).
  The test hook `RoBMA_corrected_sampler_test` reports the largest `M` per site
  after adaptation (`likelihood_magnitude`); tests in
  `test-00-selection-corrected-block.R` (block, gated, scalar, both kernels,
  the limit's boundary pair), `test-00-selection-corrected-weights.R`,
  `test-00-location-block-native.R` and `test-00-location-scale-native.R`.
- BayesTools samplers (`R/sampler-declarations.R`): `.fit()` passes
  `sampler_gates` and `sampler_blocks` to `BayesTools::JAGS_fit()`, which writes
  marker nodes that the `BayesTools::GatedSlice` and `BayesTools::FactorSlice`
  samplers claim; the model, priors and posterior are unchanged. Gates
  (`sampler.gates`, default on): every fit whose bias prior is a bias mixture.
  Blocks (`sampler.blocks`, default on): fits with a selection likelihood (not
  all sources conditioned) and an estimated heterogeneity parameter, as a
  function of the combined prior list chosen by prior-object predicates (no
  node names): a single model blocks every eligible prior (simple continuous
  scalars and treatment/independent factors without random-effect metadata, the
  single weight function, the plain allocation SD and Dirichlet allocation
  weights, and directly sampled scalar random-effect SDs (a random intercept's,
  or the one SD shared by the levels of `id()`, `diag(hom = TRUE)`, `cs()`,
  `ar1()`, `car()`) and the scalar correlations of `cs()`, `hcs()`, `ar1()`,
  `har()`, `car()`, with simple continuous scalar priors outside an
  allocation) except the unit-level latents `theta`, `gamma`, `sampling_z`,
  `pi`, `phi`; an ensemble blocks the allocation SD and weights only. Mixtures,
  spike-and-slab, ordered, point and discrete priors, inclusion priors and other
  random-effect priors (factor or vector SDs, summaries) are never declared;
  LKJ correlations are not prior-list entries. BayesTools accepts the same
  random-effect leaves (direct SD or correlation metadata, no summary, no
  allocation) and refuses the others with an error.
  The values of `selection.sampler` are tiers: `"slice.coarse"` (corrected
  sampler and declarations), `"slice"` (declarations only) and `"jags.default"`
  (no declaration and no corrected sampler in any model, as `"default"` in
  4.1.1). Under `"slice.coarse"` gated selection weights go to the corrected
  gated block and other gated nodes stay with the gated sampler, a
  declared block goes to the corrected block sampler when it has a
  coarse-capable selection child, to the factor slice sampler otherwise, and
  eligible indicators go to the screened indicator sampler;
  `sampler.gates = FALSE` removes the gates the gated block reads.
  The released `"coarse_corrected"` and `"default"` are deprecated aliases of
  `"slice.coarse"` and `"jags.default"`; settings frozen in fits keep their
  string and `.selection_sampler_corrected()` reads both spellings.

Joint JAGS updates use the bound prior sources and the selection execution
plan that generated the observations. A joint factory must claim every
stochastic coordinate it writes. Check its actual assignment on the compiled
fit and on a saved-fit extension. The scientific likelihood includes every
selection normalizer, the original covariance, and the declared prior mass of
each label. Use proposal approximations only with an exact target correction.
The Gaussian moderator update requires the exact conjugate observation graph.

For mixing comparisons, keep the previously enabled gated and corrected
samplers in the control. Charge model compilation, adaptation, burn-in and
monitored sampling. Count physical crossings separately from the Monte Carlo
error of a posterior model probability. A recorded conditional probability
belongs to the update at which its indicator was drawn. Sum disjoint labels
from that one update to estimate a family probability. Use the monitored
joint indicators for events spanning different updates.

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

The opt-in bounded post-fit CDF primitive may return endpoints for finite
scores with absolute value at least eight: the Gaussian Mills bound is below
the existing uniform CDF error. Preserve positive-blend/product error
propagation, the normal-minimum floor guard and the ordinary exact-mode retry.
JAGS and sampling copies do not opt into this approximation. Raw-rule cache
identity includes the primitive policy version.

Singleton mixture components may use the existing Gaussian recurrence only
with normal intermediate values and a transactional output buffer; otherwise
retain the original scalar expression. Compression, omission and tail bounds
remain charged, and partial failed recurrences must never enter the row sum.

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
