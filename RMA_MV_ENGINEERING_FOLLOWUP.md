# `rma.mv()` parity engineering follow-up

Date: 2026-07-13

This report follows the exhaustive
[`metafor::rma.mv()` / `RoBMA::brma.mv()` documentation parity audit](RMA_MV_DOCUMENTATION_PARITY_REPORT.md).
The original tables remain a reproducible baseline. This document records the
cross-package root-cause analysis, implemented fixes, post-fix evidence, and
remaining engineering work in RoBMA and BayesTools.

## Outcome

| Finding | Classification | Owner | Current state |
|---|---|---|---|
| Six-outcome `UN` models fail with `Invalid parent values` | BayesTools runtime defect | BayesTools | Fixed and regression-tested |
| Changes to included native `.cc` files can leave a stale DLL | BayesTools build defect | BayesTools | Fixed in Unix and Windows make rules |
| 113- and 678-level `CS` blocks generate enormous dense graphs | BayesTools architecture/resource defect | BayesTools | Fail-fast guard implemented; sparse compiler still required |
| Fixed-correlation network and `GEN` cases mix poorly | MCMC parameterization and audit-prior issue | Both | Covariance translation verified; sampler improvements recommended |
| McCurdy 1,653-row nested model | Not a structural defect | RoBMA | Deepest one-to-one block is already marginalized |
| Initial known-`R` mapping concern | Audit harness error | Audit tooling | Corrected; no package defect found |
| Orphaned known-`R` vignette fragment | RoBMA documentation defect | RoBMA | Removed |
| Random-effects vignette fails from a built source package | BayesTools packaging defect | BayesTools | Fixed and source-tar tested |
| Undeclared test dependency, unqualified utility call, and non-portable fixture path | BayesTools release hygiene | BayesTools | Fixed; clean package check |
| Random-effect labels split generated Markdown tables | Audit reporting defect | Audit tooling | Renderer fixed; reports regenerated |

The statistical parity investigation found two new package defects and one major
architectural limitation. The subsequent source-package gate exposed additional
release-hygiene defects, now fixed. The audit did **not** find evidence that the
tested `CS`, `HCS`, fixed-rho network, or `GEN` covariance algebra differs from
`metafor`.

## 1. LKJ deterministic-transform failure

### Affected examples

`help-011`, `help-012`, and `help-050` all reduce to the same six-level
`dat.craft2003` unstructured random-effect model. The original run repeatedly
failed at:

```text
mu__xREx__study_xRE_CORx_L_flat Invalid parent values
```

The exact failure reproduced with three chains, seed 1010, and 1,000 adaptation
iterations. All ten repeated restarts failed before the fix.

### Root cause

BayesTools represents an LKJ correlation matrix with stochastic canonical
partial-correlation coordinates `u` and deterministic `bt_lkj_cholesky()` /
`bt_lkj_corr()` descendants. The stochastic density correctly has open support
`0 < u < 1`. The deterministic descendants reused that strict validation.

JAGS may evaluate descendants for a finite proposal before the stochastic node
rejects an out-of-support value. A boundary or slightly outside proposal thus
raised a fatal graph error instead of receiving zero prior density. Accepting
only the exact endpoints was insufficient: adaptation can evaluate finite values
slightly beyond them.

### Implemented BayesTools fix

- Deterministic LKJ transforms now accept every finite proposal.
- Transformation-only evaluation clips each coordinate to `[0, 1]`, producing a
  finite descendant while JAGS rejects the stochastic proposal.
- LKJ density evaluation, R-side validation, bridge-sampling support, and random
  generation retain strict open support. The target distribution is unchanged.
- Regression tests cover `u = -1e-8`, `0`, `1`, `1 + 1e-8`, and a mixed
  three-dimensional boundary vector. They verify finite matrices,
  `R = L L'`, and a unit diagonal.

### Post-fix evidence

The exact three-chain/1,000-adaptation reproduction now initializes and samples.
A full audit fit (1,000 adaptation, 2,000 burn-in, 5,000 retained iterations per
chain) produced:

| Coefficient | `metafor` | `brma.mv()` posterior mean |
|---:|---:|---:|
| 1 | -0.0600 | -0.0364 |
| 2 | -0.1423 | -0.1080 |
| 3 | 0.3167 | 0.3760 |
| 4 | 0.5671 | 0.5535 |
| 5 | -0.4888 | -0.4155 |
| 6 | -0.4750 | -0.4203 |

The maximum coefficient difference was 1.44 `metafor` standard errors. This is
runtime-correctness evidence, not a final equivalence claim: coefficient R-hat
reached about 1.12, full-object R-hat 2.84, and minimum ESS 24. The crash is
fixed; this broad-prior, non-centered model still needs longer or better-tuned
sampling.

## 2. Native build invalidation

BayesTools compiles subdirectory C++ sources through an amalgamating
`native-sources.cc`. Make did not know that the resulting object depended on the
included `.cc` files, so editing the LKJ implementation could report “Nothing to
be done” and leave an old DLL loaded.

Both `src/Makevars.in` and `src/Makevars.win.common` now:

- make the shared library depend on `native-sources.o`; and
- list every included implementation file as a dependency of that object.

This was necessary for the LKJ regression to test the edited source rather than
a stale binary.

## 3. High-cardinality structured random effects

### Corrected model classification

`help-037` is a 678-level compound-symmetry model, not a 678-effect
unstructured model. Its mapping is:

| `metafor` | `brma.mv()` |
|---|---|
| `random = ~ factor(Outcome) \| factor(SampleID)` (default `CS`) | `random = ~ cs(Outcome \| SampleID)` |

`help-090` is likewise `CS`, with 113 outcome levels and 17 groups.

### Why the current graph is unsafe

The current BayesTools compiler globally one-hot expands the structured factor,
samples a dense `groups x levels` latent matrix, symbolically constructs full
`levels x levels` Cholesky/correlation matrices, performs dense transforms, and
monitors large deterministic arrays by default.

| Example | Dimensions | Symbolic Cholesky products | Dense transform products | Monitored random-block values* | Approx. retained doubles** |
|---|---:|---:|---:|---:|---:|
| `help-090` (Tanner) | 17 x 113 | 240,464 | about 1.46 million | 27,459 | about 3.1 GiB |
| `help-037` (Bakdash) | 79 x 678 | 51,944,179 | about 36.8 million | 972,930 | about 108.7 GiB |

\* Guard estimate for latent and deterministic matrix arrays; a few scalar
hyperparameters add to the actual monitored-column count.

\** Three chains x 5,000 retained draws, raw eight-byte doubles only. R/JAGS
objects, copies, indexing, and graph memory increase the real requirement.

The original 60-minute bounds were therefore symptoms of graph construction and
storage scale, not evidence against statistical equivalence.

### Implemented guard

BayesTools now estimates symbolic Cholesky work, dense transform work, and
monitored values before emitting sampled random-effect syntax. It stops with the
block name, structure, dimensions, estimates, limits, and remediation guidance
when the dense representation is unsafe.

The exact audit mappings now fail immediately and actionably:

- `help-090`: `CS; 17 groups x 113 columns`, including
  `cholesky_products=240464` and `monitored_values=27459`;
- `help-037`: `CS; 79 groups x 678 columns`, including
  `cholesky_products=51944179`, `transform_products=36774720`, and
  `monitored_values=972930`.

Experts who have independently checked resources can bypass the guard with:

```r
options(BayesTools.random_effects_complexity_multiplier = Inf)
```

This is an escape hatch, not a performance fix.

### Required compiler work

Priority order:

1. Stop retaining deterministic full `L`/`R` matrices by default; reconstruct
   summaries from `rho` or LKJ coordinates where compatible with the public raw
   posterior contract.
2. Add closed-form `CS`/`HCS` Cholesky generation and direct one-hot lookup,
   removing cubic symbolic cross-products.
3. Compile each group against only its observed factor-level subset for
   `CS`/`HCS`. The audit data reduce local squared dimensions from 12,769 to 983
   (`help-090`) and from 459,684 to 18,288 (`help-037`).
4. Extend sparse compilation carefully to `AR`/`HAR`/`CAR`; ordering and distance
   semantics must be retained.
5. Keep unrestricted `UN` blocks small or guarded. Marginal subset compilation
   under an LKJ prior is not automatically prior-equivalent and needs a separate
   design.

For `CS`, a closed-form Cholesky factor supporting the full admissible negative-
rho range was checked numerically against direct decomposition to machine
precision (maximum error `2.22e-16`). It should replace generated cubic sums in a
dedicated implementation.

## 4. Poorly mixing fixed-correlation and `GEN` cases

The fixed-rho network cases (`help-071` through `help-074`, plus `help-087`) and
the `GEN` random-slope case (`help-083`) were rechecked at the covariance level.
For every case, the BayesTools design and prior-fixed covariance reproduce the
`metafor` marginal random-effect covariance:

```text
max abs(Z G Z' - (M - V)) = 2.22e-16   # network cases
max abs(Z G Z' - (M - V)) = 3.55e-15   # GEN case
```

These are not covariance-translation defects. The practical causes are:

- an always-non-centered, dense `groups x coefficients` latent representation;
- unused group/column cells that remain in the graph;
- very broad audit location priors relative to the predictor scale; and
- prior-drawn initial values, reaching roughly +/-60 in the network cases and
  substantially larger values for the unscaled `GEN` slope.

Recommended work:

- Audit harness: choose priors on the linear-predictor scale, initialize fixed
  effects from zero or deterministic GLS estimates, and label a fit successful
  only after explicit R-hat/ESS thresholds pass.
- BayesTools: provide centered/non-centered/automatic random-block
  parameterization; prune unobserved latent cells; and support marginalizing
  Gaussian random blocks into `V + Z G Z'` where the likelihood permits it.
- RoBMA: surface failed convergence as a first-class fit/audit status rather than
  relying on the generic “fitted” state.

## 5. Nested and known-covariance cases

The 1,653-row McCurdy nested example contains five scalar-intercept blocks with
1,653, 582, 310, 126, and 804 levels. RoBMA already marginalizes the deepest
one-to-one 1,653-level block. The remaining sampled-group count is 1,822 and the
known sampling covariance stays compact and diagonal. This is linear graph
growth, not the high-cardinality structured-factor failure above.

The initially suspected known-`R` discrepancy came from an audit mapper retaining
`NULL` list entries. Correcting the mapper removed it; no RoBMA/BayesTools defect
was established.

## 6. Documentation correction

The multivariate parity vignette contained an incomplete, unfenced `predict()`
fragment referring to undefined `known_r_rows`, `fit_known_r`, and `study_R`
objects. No complete known-`R` example exists in that vignette, so the orphaned
fragment was removed. The four documented worked examples and their cached fit
objects remain unchanged.

## 7. BayesTools source-package hardening

The first clean `R CMD check` validated compilation, installation, examples, and
tests, then exposed unrelated source-package defects. They were fixed rather
than suppressed:

- The random-effects vignette sourced a helper from `tools/`, while `.Rbuildignore`
  correctly excludes that directory. The canonical helper now lives beside the
  vignette. All 60 cache-dependent non-setup chunks use `purl = FALSE`, so R's
  separate vignette-code check sources only the safe setup when the intentionally
  excluded 15 MB model cache is absent. Cache-free render, tangle, and source
  checks all pass.
- `withr`, used by fixture tests, is now declared in `Suggests`.
- `object.size()` is called as `utils::object.size()`.
- A malformed zero-byte top-level file named `2L){` was removed.
- The sole packaged path over 100 bytes was shortened by renaming the internal
  fixture ID from `fit_bias_petpeese_heterogeneous_wf` to
  `fit_bias_petpeese_hetero_wf`. The included path is now 94 bytes.

The fixture-ID rename intentionally makes an existing local fixture cache stale.
The cache was preserved; regenerate `test-00-model-fits.R` before running the
fixture-profile summary-table tests. Unit tests and source-package checks do not
depend on that local cache.

The final BayesTools check completed with 0 errors, 0 warnings, and 0 notes. On
Windows, `R CMD build` still reports that it corrected executable permissions for
`configure` and `cleanup`; verify their repository mode bits on Unix before a
release. This does not affect the Windows check result.

## 8. Verification performed

| Repository / target | Result |
|---|---|
| BayesTools complete unit profile | 6,390 pass, 0 fail, 0 warn, 7 skip |
| BayesTools compiled LKJ fit file | 778 pass, 0 fail, 0 warn |
| BayesTools structured-design file after guard override test | 1,134 pass, 0 fail, 0 warn |
| BayesTools source-package check | 0 errors, 0 warnings, 0 notes |
| RoBMA `00-input-data-mv` | 572 pass |
| RoBMA `01-brma.mv` | 128 pass |
| RoBMA `02-brma-mv-metafor` | 156 pass |
| RoBMA `02-heterogeneity-mv` | 98 pass |
| RoBMA `02-random-parameters` | 78 pass |
| Exact six-level UN reproduction | initializes and samples after fix |
| Exact `help-037` / `help-090` mappings | immediate, dimension-specific guard errors |

RoBMA tests were run against the parent-directory BayesTools source tree, so they
exercise the cross-package changes together.

## 9. Release recommendation

For the next coordinated release:

1. Ship the LKJ deterministic-transform and native dependency fixes in
   BayesTools, with the new boundary regression.
2. Ship the dense-complexity guard in BayesTools and document its expert option.
3. Ship the RoBMA vignette and parity-artifact corrections.
4. Do not claim support for the 113/678-level examples as practically fitted
   until sparse `CS`/`HCS` compilation and posterior-monitor reductions land.
5. Treat sampler parameterization and convergence status as the next correctness-
   adjacent milestone; the algebra is correct, but the current runtime behavior
   can produce nominal fit objects with unusable diagnostics.
6. Regenerate the renamed fixture cache and verify `configure`/`cleanup` mode bits
   from a Unix checkout before release packaging.

Spatial structures, arbitrary `W` weights, and remaining partial fixed-parameter
API mappings listed in the baseline syntax crosswalk remain feature gaps, not
newly discovered regressions.
