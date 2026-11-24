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
  children's log density is `GatedSlice`'s declaration error. Gate open: one
  update of the shared block kernel, whose axes learn from these updates only;
  until they are first refreshed (100 open adaptive updates) the coordinate
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
  weights included), proposes by Liu's Metropolized Gibbs sampler on `exp(s)`
  (a level other than `k` with probability `p_k' / (1 - p_k)`, accepted on the
  surrogate with `min(1, (1 - p_k) / (1 - p_k'))`; the identity when no other
  level has mass), a move reversible for the surrogate law, and for `k' != k`
  evaluates the exact scores of `k` and then `k'` (so the last exact
  evaluation is the proposed state) and accepts with
  `min(1, exp((e_k' - s_k') - (e_k - s_k)))`; a self-proposal costs no exact
  evaluation. The plain surrogate draw from `exp(s)` (current level included)
  with the same correction stays behind the test hook (`finite_proposal = 0`).
  Liu's proposal is the default: on the tutorial and Assink fits it raised the
  bias label's switches by 37% to 2x for 2-10% more time. The exact
  kernel is base::Finite's enumeration, arithmetic and single uniform. The
  `finite_proposal = 1` mutant (no correction) fails the exactness fixtures.
  Batch scoring of the six selection labels was not built.
- Row cache of step-switch rows (`src/samplers/StepRowCache.{h,cc}`): the
  surrogate term of the singleton `dselnorm_step_switch` children of the
  corrected gated block and of a screened indicator that reaches their weights
  or kernel mode (`CorrectedBlockOptions::row_cache`, the kernel's
  `row_cache` argument; the declared-block sampler and the scalar corrected
  slice keep these rows as exact terms). In step mode with observed bin
  `o > 0` a row's density is `pw [log phi(sign y; sign m, sd) + log omega_o -
  log A]`, `A = sum_b omega_b p_b`; while only the weights or only the label
  move, `p_b`, both Gaussian terms (normal and step mode) and the fixed
  parameters' checks are constants, so `refresh()` at the start of each
  corrected update (and of the screened test scores) reads each row once and
  caches them, and `value()` reads per row the kernel mode, the weights (their
  validity once per weights vector) and the mean and SD, in the order of
  `cpp_selnorm_kernel_lpdf`. A row goes to its node's own `logDensity` when a
  parameter other than the mean, SD, weights and mode is not fixed or fails
  the distribution's check, its SD or SE is infinite, no bin probability
  exceeds 1e-280 (at the refresh), its mean or SD differs bitwise from the
  refresh (PET and PEESE labels), or its cached `A` is not above 1e-280
  (StepRowNormalizer's linear threshold; the log-scale rows). Telescoping
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
  affine guards) and declines otherwise. Each update rebuilds the Gaussian
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
  sampler record carries an evidence record per gate and chain (flips by
  direction, screened proposals, log odds of the other state) and, for the
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
  is Liu's Metropolized Gibbs step on the label's surrogate full conditional;
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
  densities. The full ratio of the reverse move equals the inverse of the
  forward ratio up to the rounding of the refitted Laplace laws (at most
  8.1e-8 relative over 13 chain seeds of each fixture variant: the base
  point rebuilt from either end agrees to rounding, but the damped Newton
  fit can take a different last step). L2 runs with probability 1/4 per sweep (`default_swap_probability`;
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
  by the exact correction, apart from the identical-state shortcut, so the
  mutants remove terms from the ratio only); the trace switch records every
  L2 attempt (`RoBMA_location_block_test_trace`). The evidence record adds
  `evidence$bias` (per chain the proposals and acceptances of L0, L1, S and
  L2, and for L2 by direction the count, sum and maximum of the stage-1 log
  ratio of all attempts and of the swaps the exact correction kept) and
  `evidence$bias_labels` (sweeps ending in each label): a family that is
  never left although swaps out of it pass stage 1 points to a stuck chain,
  while very negative ratios throughout are consistent with, but do not
  prove, decisive evidence. Exactness fixtures (labels none, selection, PET,
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
  level. The factory declines (everything as before; draws identical) unless
  the sources are free, unobserved, fixed-parameter scalars with the declared
  supports (a gamma group recognized by `find_gamma_groups` behind the
  declared ratios), the block owns all `K` latents, no source has a gate
  marker or shares a block marker with another node, the sources' samplers
  claim exactly them, and the invariance check passes: once at construction
  (chain 0, state restored) with the indicators forced to their levels, the
  latents at nonzero values and every bias label, for `c` in {0.7, 1.3}
  every parent value of every stochastic child of the sources and latents
  (markers excluded) is unchanged within `1e-9 (1 + |value|)`, the charts of
  `a`, `w` and `sigma_s` match and `sigma_s' = c sigma_s`. The tag passes the
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
