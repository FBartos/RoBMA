# BayesTools prior-density ordinate contract

## Purpose

RoBMA needs to classify a requested qCMDE/IWMDE prior ordinate before any
posterior-density computation. The classification explains whether the prior
density at the exact requested value is regular, zero, infinite, undefined, or
contains a point mass.

BayesTools already owns primitive priors, induced linear densities,
transformations, mixture construction, and their provenance. The classifier
therefore belongs in BayesTools. RoBMA must not reconstruct an induced prior
from posterior samples or duplicate BayesTools' density algebra.

## Required exported interface

Add an exported function with this conceptual signature:

```r
prior_density_ordinate(x, value)
```

Initial accepted `x` classes:

- scalar simple, none, point, scalar mixture, and scalar spike-and-slab
  BayesTools prior objects;
- `prior_linear_density`, including objects produced from a
  `prior_density_context` by existing deterministic BayesTools builders.

Non-scalar/vector, factor, Dirichlet, ordered, weightfunction, and other
component-valued priors require an explicit scalar marginal first. They should
arrive as a scalar `prior_linear_density` with its weights and provenance.
Valid but unsupported prior classes return `unknown`; malformed inputs error.

`value` must be one non-missing numeric scalar. Infinite values may be rejected
unless the input class has an explicitly defined ordinate there.

Return one named object with a stable schema:

```r
list(
  schema_version = "1",
  value          = value,
  behavior       = "regular",
  log_density    = -1.23,
  point_mass     = 0,
  exact          = TRUE,
  method         = "primitive",
  reason         = NULL,
  provenance     = list(...)
)
```

The result may have class `prior_density_ordinate`, but RoBMA must be able to
consume the named fields without relying on print output.

### Field semantics

`schema_version` and the field names/order are stable. `behavior` is exactly
one of:

- `regular`: the continuous density has a finite, strictly positive
  mathematical limit at `value`;
- `zero`: the continuous density is mathematically zero or tends to zero at
  `value`;
- `infinite`: the continuous density tends to positive infinity at `value`;
- `point_mass`: positive discrete probability exists exactly at `value`;
- `undefined`: the mathematical ordinate is not defined, for example because
  the requested transformation is undefined there;
- `unknown`: BayesTools cannot establish the mathematical behavior from
  deterministic provenance in the supported scope.

`log_density` is the continuous log-density ordinate when BayesTools can
evaluate it without changing the requested value. It is `-Inf` for a proven
zero, `Inf` for a proven singularity, and `NA_real_` when unavailable or not
meaningful. A representationally underflowed `-Inf` must not by itself change a
structurally established `regular` classification to `zero`.

`point_mass` is the exact discrete probability at `value`, or zero when none is
present. If it is positive, `behavior` is `point_mass`; provenance may retain
the continuous-component classification for diagnostics, but RoBMA will not
interpret a mixed-measure point as an ordinary Savage-Dickey density ordinate.

At an exact support endpoint, “limit” always means the one-sided limit from
inside the support. Thus, for example, an ordinary half-normal density is
regular at zero. Values outside the exact support are zero.

`exact` states whether `behavior` follows from the prior definition and exact
transformation/mixture provenance. Probabilities and parameters are still
represented and normalized with ordinary binary64 arithmetic; `exact` does not
promise rational or arbitrary-precision evaluation. `unknown` always has
`exact = FALSE`.

`method` is one enumerated machine-readable label: `primitive`, `point`,
`finite_mixture`, `scalar_affine`, `linear_normal`, `named_transform`, or
`unsupported_provenance`.

`reason` is `NULL` for an ordinary regular result and otherwise a concise
diagnostic suitable for inclusion in a RoBMA warning.

`provenance` must be deterministic and compact. It should identify the source
prior family, transformations, linear weights, and mixture components needed to
audit the classification. It must not contain draws, large numerical grids, or
closures whose environment captures fitted objects.

## Required mathematical rules

Classification is structural. Never infer mathematical behavior from:

- a sampled prior or posterior distribution;
- KDE output;
- a finite numerical grid or interpolation;
- FFT clipping or a zero outside the current approximation range;
- probing `value +/- epsilon`;
- binary64 underflow or overflow alone.

The first implementation should classify only cases whose behavior can be
established cleanly:

1. Supported primitive continuous priors, including truncation boundaries,
   support exclusion, and known interior zeros/singularities.
2. Point priors and finite scalar mixtures, with exact component weights and exact
   point locations.
3. Nonzero scalar affine transformations, including the absolute Jacobian.
4. Linear combinations that reduce analytically to a normal or multivariate
   normal marginal.
5. Existing named monotone transformations (`lin`, `exp`, `exp_lin`, `tanh`)
   only when their inverse, domain, and Jacobian limit are structurally known
   for the supplied arguments.

For finite mixtures:

- positive point mass at the exact requested value yields `point_mass`;
- zero-weight components are ignored;
- without point mass, precedence is `undefined`, then `infinite`, then a proven
  finite-positive `regular` sum, then `zero` only when every positive-weight
  continuous component is proven zero, then `unknown`;
- component ordinates with very different magnitudes are combined in log
  space.

A structurally regular density remains `regular` if its finite mathematical
log-density lies outside binary64 range. In that case return the representable
result (`-Inf` where appropriate), keep `exact = TRUE`, and set `reason` to
explain that the behavior—not the floating-point ordinate—is authoritative.

General convolutions, products, arbitrary user transformations, and transformed
or induced densities whose boundary limits are not analytically encoded must
return `unknown`. A useful numerical density estimate may be retained separately
in `log_density`, but it must not promote `unknown` to a structural class.

## Existing BayesTools machinery to reuse

Build on the deterministic objects and provenance already produced in:

- `R/priors-density-context.R`;
- `R/priors-linear-density.R`;
- `R/priors-linear-density-combinations.R`;
- `R/priors-density.R`.

In particular, reuse exact point locations, transformation definitions, and
compact adaptive-evaluation arguments/context. Extend this metadata where
necessary instead of adding a parallel density system.

Current `singular_density_points` markers are not authoritative: some are
created after inspecting grid-interpolated heights, and transformation code does
not always transform the marked location. Regenerate singular markers from
structural component provenance and propagate their exact transforms, or ignore
legacy markers for classification. Test both a product singularity and its
transformed location.

Current `density_evaluator` closures call density-scale `mpdf()` and may
underflow. They are insufficient classification metadata and must not be
returned in provenance. Prefer a compact recognized analytic specification or
the existing adaptive arguments/context. Mixture builders currently discard
some component metadata; mixture classification must inspect deterministic
component provenance rather than the collapsed grid.

For row-varying contexts, provenance may record dimensions, unique-row count, a
deterministic hash, and normalized component classifications. It must not retain
the raw weight matrix.

The existing `.prior_linear_density_height()` remains useful for numerical
height evaluation, but its interpolated/grid result is insufficient evidence
for `zero`, `infinite`, or `regular` classification.

## RoBMA integration contract

After this function is released, RoBMA will:

1. retain the deterministic scalar target-density object/provenance when constructing
   an IWMDE parameter specification;
2. call `BayesTools::prior_density_ordinate()` once per unique deterministic
   target density and exact requested ordinate; BayesTools recursively handles
   positive-weight mixture components and may expose branch results only in
   compact provenance;
3. run this preflight before baseline log-q, likelihood, or posterior-ordinate
   estimation;
4. issue a targeted warning/error immediately for `zero`, `infinite`,
   `point_mass`, `undefined`, and `unknown` when an ordinary point Bayes-factor
   ordinate was requested;
5. preserve the exact requested value and never substitute an epsilon;
6. keep the compact classifications in IWMDE plan/diagnostic metadata.

RoBMA will not reconstruct a missing induced density from prior or posterior
draws. Until the exported interface is available, RoBMA retains its current
primitive-prior warning and does not add a second induced-density engine.

## Tests and acceptance criteria

Use focused deterministic unit tests; no JAGS fits or large numerical grids are
needed. Cover at least:

- regular primitive interior;
- value outside support (`zero`);
- finite, zero, and infinite primitive boundary limits;
- half-normal/truncated-normal at the lower boundary (`regular`), and uniform
  endpoints under the support-relative one-sided rule;
- known interior nonlocal-prior zero;
- exact point mass;
- Bernoulli, `prior_none`, and unsupported direct vector/component priors;
- finite mixtures covering zero weights and point/infinite,
  undefined/infinite, and regular/zero precedence;
- positive, negative, and zero scalar affine weights; zero is a point mass;
- a deterministic point shift plus normal terms reducing to an analytic normal;
- analytic normal linear combination;
- supported named transform at interior, outside-support, and boundary values;
- `lin`/`exp_lin` with finite nonzero `b`; zero slope is classified as the
  corresponding point transform; `exp`/`exp_lin` range `(0, Inf)`; `tanh` range
  `(-1, 1)`; nonfinite transform arguments and source-dependent boundary limits
  return `undefined` or `unknown` as mathematically appropriate;
- arbitrary transform and general numerical convolution returning `unknown`;
- stale/manual `prior_linear_density` without adaptive provenance returning
  `unknown`, not error;
- legacy grid zeros, FFT clipping, and current approximation-range exclusion
  never promoting a structural classification;
- a mathematically regular extreme normal ordinate that underflows in ordinary
  density space but is not misclassified as `zero`;
- transformed structural singularity appearing at its transformed location;
- serialization showing fixed fields/order and compact provenance containing no
  functions, environments, grids, draws, or raw row-weight matrices.

RoBMA integration tests, after the API exists, must also verify that preflight
diagnostics occur before a mocked plan/estimator failure, are deduplicated across
adaptive repetitions, leave the requested coordinate unchanged, and do not
traverse posterior indicators to recreate prior branches.

Acceptance requires:

- no posterior/prior sampling in the implementation;
- no epsilon probes or heuristic scale thresholds;
- no arbitrary-precision, double-double, or new native arithmetic;
- no classification based solely on grid/interpolation height;
- no new large test matrix or model-fitting tests;
- documented exported API and regenerated namespace/documentation.

## Explicit non-goals

- symbolic classification of every convolution, product, or user function;
- exact arithmetic for all representable binary64 values;
- replacing BayesTools' existing plotting-density approximations;
- deciding whether RoBMA should report a Bayes factor after a nonregular result;
- backwards-compatibility adapters for undocumented internal objects.
