# Maintainer decisions pending

## Numerical-faithfulness visual baselines

No committed visual baseline was regenerated automatically. The standard test
profile retains representative, human-reviewed visual tests; the cases below
remain in the certification gallery until reviewed.

### Selected-normal funnel and Q-Q plots

- Issue: exact selected-normal boundary evaluation changes the top funnel
  endpoints by less than one SVG pixel. Exact paired PIT-tail evaluation moves
  one Q-Q point by 0.01 SVG pixel.
- Impact: certification snapshots
  `funnel-selmodel-pos-brma-ggplot.svg` and `qqnorm-selmodel-base.svg` fail;
  selection-model behavioral tests and the base funnel comparison remain in the
  standard profile.
- Suggested change: visually confirm the `.new.svg` files, then accept these two
  baselines. The differences follow the intended removal of endpoint and PIT
  clamps.

### Radial mathematical labels

- Issue: the current R graphics toolchain renders plotmath fractions a fraction
  of a pixel differently. No radial source changed, and the same displacement is
  present in both metafor and RoBMA panels. Ten default-label radial snapshots
  differ; the two plain-label snapshots are unchanged.
- Impact: default-label variants fail in certification. The standard profile now
  runs the unchanged plain-label base and ggplot snapshots plus all radial data,
  validation, and dispatch assertions.
- Suggested change: visually confirm and accept the ten `.new.svg` files as a
  graphics-toolchain refresh. Do not add coordinate-tolerance canonicalization;
  it could hide real plot regressions.

### Metafor forest gallery

- Issue: current metafor emits 9,999 prediction-shade rectangles instead of 99.
  A separate prediction-interval snapshot changed from `[0.14, 1.75]` to
  `[0.15, 1.72]` because its posterior prediction is stochastic.
- Impact: `forest-simple-prediction-shade.svg` and
  `forest-simple-ilab-columns.svg` fail only in certification. Stable forest
  comparison, prediction-bar, and fitted-family snapshots remain standard.
- Suggested change: keep the dense shade variant release-only. Make the ilab
  prediction deterministic before accepting a new human-reviewed baseline;
  do not accept a stochastic snapshot as a correctness oracle.

## Cached vignette fits

- Issue: the 16 vignette cache directories contain about 220 MB of fitted
  objects created before BayesTools required its canonical parameter registry.
  Current summaries and plots reject those objects and request a refit.
- Impact: `R CMD build` cannot rebuild the cached vignette output, so the
  default `devtools::check()` stops before package checking begins. Standard
  package tests use current fits and are unaffected.
- Suggested change: regenerate the committed vignette fits with the current
  RoBMA and BayesTools versions, one vignette at a time, then review the rebuilt
  output. Do not add an in-memory cache migration or compatibility layer: it
  would infer metadata for old fitted objects and could conceal stale model
  semantics. Until regeneration, use `devtools::check(vignettes = FALSE)`
  to verify the remainder of the package.

## Covariance validation before analysis subsetting

- Issue: `V` is validated after explicit subsetting and missing-row removal.
  An invalid block that is entirely excluded from the analysis is therefore not
  inspected.
- Impact: the fitted covariance is always validated, but the package does not
  certify every unused entry supplied in the original matrix.
- Suggested change: decide whether `V` describes the submitted dataset or only
  retained observations. If the former, validate the complete matrix before any
  row removal; if the latter, document the current retained-submatrix contract.

## IWMDE nonregular-prior warning scope

- Issue: stable log-density classification covers primitive supported priors,
  but not every transformed, induced-linear, or non-simple prior density. The
  warning is also emitted only after a successful estimate.
- Impact: some nonregular ordinates can fail before receiving the more specific
  zero/infinite/undefined prior-density explanation.
- Suggested change: add an explicit induced log-density/provenance interface and
  run its classification as a preflight step. Do not reconstruct induced
  densities from posterior samples.

## Selected-normal extreme-scale contract

- Issue: CDF/PIT and normalizers now use paired score-space tails, but moment,
  RNG, threshold, and active p-hacking paths still contain some division-first
  coordinates at scales near binary64 overflow.
- Impact: ordinary statistical ranges are covered; a package-wide guarantee for
  every representable extreme input is not yet defined.
- Suggested change: decide which public selected-normal operations promise
  extreme-scale support. Extend only those kernels with focused reference tests;
  do not restore a general exact-arithmetic subsystem.

## STEP fast-route selection

- Issue: `.selection_telescope_probabilities()` routes by fixed `z = 8` and
  `omega = 100` thresholds.
- Impact: the thresholds choose between equivalent numerical algorithms, not a
  statistical model, but their accuracy/performance contract is implicit.
- Suggested change: replace them with either a documented operation-specific
  error certificate or one canonical log-partition route. Do not tune the
  constants without such a contract.

## Cluster selected-normal quadrature policy

- Issue: clustered selected-normal AGHQ uses fixed finite-difference and
  curvature floors (`1e-3` and `1e-8`) without a public accuracy or
  nonconvergence policy.
- Impact: these choices can affect difficult boundary/singular cluster
  likelihoods and cannot be classified as performance-only fallbacks.
- Suggested change: specify quadrature order, convergence checks, error
  reporting, and boundary behavior first; then replace fixed floors with those
  rules and certification cases.

## Exact transformed-coordinate qCMDE/IWMDE

- Issue: sample-based affine detection was removed because it could silently
  change density ordinates and Bayes factors. qCMDE/IWMDE now require the fitted
  coefficient coordinate unless a dedicated exact transformation exists.
- Impact: automatically unscaled formula coefficients can require
  `standardized_coefficients = TRUE`; KDE remains available on display scales.
- Suggested change: expose structural coefficient maps from BayesTools and
  compose them into IWMDE primitive/linear parameter specifications. Restore
  transformed-coordinate support only from that provenance, including the exact
  Jacobian.
