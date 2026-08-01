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
