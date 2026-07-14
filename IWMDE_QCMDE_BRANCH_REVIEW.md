# IWMDE/qCMDE Branch Review

**Branch:** `BF-coef-test`
**Comparison base:** `master` at `09d07b90`
**Review date:** 2026-07-13
**Scope:** current branch implementation of qCMDE/IWMDE posterior densities,
point ordinates, Savage-Dickey Bayes factors, marginal means, plotting,
known-`V` integration, diagnostics, tests, documentation, and maintainability.

The findings through "Proposed decision order" are the pre-implementation
review and are retained as the decision record. The post-implementation review
at the end of this document supersedes their implementation-status and test-
coverage statements.
File/line references in that retained section describe the pre-implementation
snapshot; the module split intentionally superseded several paths.

## Executive conclusion

No critical mathematical correctness defect was found in the qCMDE or IWMDE
estimators. The implemented estimators have the right high-level form:

- qCMDE numerically normalizes each evaluated joint ordinate over the focal
  parameter and averages the resulting conditional densities over posterior
  nuisance draws.
- IWMDE evaluates the Chen importance ratio
  `w(theta_i | xi_i) q(theta_star, xi_i) / q(theta_i, xi_i)` and averages it.
- Active mixture mass, fixed denominators, one-sided support transformations,
  requested-versus-evaluated boundary values, and known-`V` joint likelihoods
  are handled deliberately rather than by post-hoc density rescaling.

The strongest independent validations were positive:

- In an exactly solvable known-`V` Gaussian model, qCMDE had density-curve
  `L1` error `2.13e-5`; IWMDE had `L1` error `0.00400`.
- In the ordinary random-effects example, qCMDE at all 15,000 draws produced
  `BF10 = 60.897`, versus `60.901` from an independent closed-form
  conditional-normal/Rao-Blackwell calculation.
- Interior multilevel regression, a known-`V` heterogeneity boundary, and a
  derived PET moderator all approached bridge-sampling results as the row
  budget increased.

The main release risks are not a wrong estimator formula. They are:

1. The default point-BF row budget can produce 10-17% error in a difficult
   tail while the current warning policy emits no warning.
2. Full reliability diagnostics are not available through the public result,
   although they decide whether an ordinate is accepted as BF-grade.
3. Existing bridge-oracle tests use log-BF tolerances of `0.25` and `0.35`,
   which are too loose to detect several meaningful regressions.
4. The subsystem is large and structurally repetitive: 14 implementation
   files, about 13,000 lines, and 353 functions. Several core functions exceed
   150-200 lines and manipulate informal list schemas.

Recommendation: approve the estimator design, but address Findings F1-F4
before treating qCMDE/IWMDE Bayes factors as release-ready for general use.

## Severity summary

| ID | Severity | Finding | Status |
|---|---|---|---|
| F1 | High | Default BF precision and warning thresholds are too permissive | Approved |
| F2 | High | BF-grade diagnostics are not publicly inspectable | Approved |
| F3 | High test risk | Bridge-oracle tests tolerate large numerical errors | Approved |
| F4 | Medium | IWMDE MCSE conditions on weights estimated from the same draws | Approved limitation |
| F5 | Medium | Stored-result source fingerprint samples only ten posterior rows | Approved |
| F6 | Medium | User documentation does not define estimator choice or reliability policy | Approved |
| F7 | Medium | Known-`V` qCMDE can be slow even for a very small fit | Approved |
| F8 | Medium | Core implementation remains too list-heavy and repetitive | Approved |
| F9 | Low | Conditional-weight fallback discards the failure reason | Approved |
| F10 | Low | Pilot/final normalization field names are ambiguous | Approved |

## Findings

### F1. High: default BF precision and warning thresholds are too permissive

**Evidence**

- `max_samples` defaults to `500` in `R/iwmde-api.R:595-603`.
- A point ordinate is not rejected until relative MCSE reaches 100% at
  `R/iwmde-api.R:1585-1589`.
- A warning is not emitted until relative MCSE reaches 25% at
  `R/iwmde-api.R:1591-1595`.
- BF ESS warning/failure thresholds are only 20 and 4, and maximum-weight
  warning/failure thresholds are 50% and 80%, at
  `R/iwmde-api.R:1597-1625`.
- In the difficult `fit_simple` tail at `mu = 0`, the 500-row qCMDE result was
  `54.12` versus the conditional-normal oracle `60.90`, an 11.1% difference.
  IWMDE was `50.63`, a 16.9% difference. Reported relative MCSEs were 13.3%
  and 13.6%, so neither result generated a warning under the current policy.
- At 15,000 rows, both estimators converged to the independent conditional
  oracle: qCMDE `60.897`, IWMDE `61.036`, oracle `60.901`.

The current MCSE correctly signals that the 500-row estimates are noisy, but
the policy labels them BF-grade without a warning. A 10-17% BF error can change
reported evidence categories and rounded output.

**Suggested action**

Use separate controls for plotting and point-BF ordinates. For a point BF,
evaluate rows in increasing deterministic budgets until a precision target is
met or all draws/the hard cap are exhausted. A reasonable initial policy to
validate empirically is:

- target relative MCSE at or below 5%;
- warn above 5%;
- reject above 20-25%;
- also require a materially larger ESS and lower maximum weight share than the
  present 4/80% failure limits;
- return the achieved row budget and whether the hard cap was reached.

If adaptive evaluation is deferred, raise the hypothesis default above 500
without changing the cheaper plotting default.

**Decision:** Approved. Implement separate plotting and point-BF precision
policies, with adaptive point-BF evaluation and stricter reliability thresholds.

### F2. High: BF-grade diagnostics are not publicly inspectable

**Evidence**

- The internal ordinate carries MCSE, ESS, maximum weight share, row loss,
  normalization diagnostics, quadrature changes, evaluation value, and method
  at `R/iwmde-api.R:708-790`.
- `hypothesis()` exposes only `BF_error` and warning text. In direct runtime
  probes, `attr(result, "posterior_ordinate")` was `NULL`.
- `plot_iwmde_diagnostics()` and
  `plot_iwmde_marginal_means_diagnostics()` are marked `@noRd` at
  `R/iwmde-api.R:54` and `R/iwmde-api.R:166`; they are neither documented nor
  exported.
- The permissive policy in F1 means absence of a warning does not imply high
  precision.

Users therefore cannot determine why an estimate passed, how concentrated its
weights were, whether a cap was exhausted, or whether qCMDE quadrature was
stable. They also cannot perform a principled sensitivity run beyond manually
changing controls and comparing BF values.

**Suggested action**

Attach a stable, compact diagnostics object to every qCMDE/IWMDE hypothesis
result and provide a public accessor, for example
`density_diagnostics(result)`. Keep the current developer plotting functions
internal until their API is simplified. The public schema should include:

- estimator and requested/evaluation value;
- row budget, eligible/evaluated/retained rows, and active mass;
- relative MCSE, ESS, and maximum weight share;
- qCMDE quadrature/ordinate movement or IWMDE mass error;
- warning/failure thresholds and hard-cap status;
- source and estimator schema versions.

**Decision:** Approved. Add stable public access to compact BF-grade diagnostics.

### F3. High test risk: bridge-oracle tests tolerate large numerical errors

**Evidence**

- The multilevel coefficient bridge checks use log-BF tolerance `0.25` at
  `tests/testthat/test-02-hypothesis.R:403-408`. A log difference of `0.25`
  permits a ratio difference of about 28%.
- The PET factor check uses tolerance `0.35` at
  `tests/testthat/test-02-hypothesis.R:1117-1120`, permitting about 42%.
- The PET test uses only 120 estimator rows. The independent probe found 10-13%
  error at 120 rows but about 1.3-1.4% error at 2,000 rows.
- The known-`V` boundary test checks availability and diagnostics but not an
  independent bridge or analytic target.
- The exact conditional-normal oracle for the ordinary random-effects example
  is absent, although it is cheap and more direct than bridge sampling.

The internal unit tests are extensive and valuable, especially the analytic
normal-normal and correlated-linear-contrast tests in
`tests/testthat/test-02-iwmde-fast-paths.R:256-481`. The missing layer is a
tighter set of end-to-end statistical oracles.

**Suggested action**

Add a small numerical certification matrix:

1. Exact known-`V` conjugate normal density and marginal likelihood.
2. Closed-form conditional-normal qCMDE ordinate for ordinary random effects.
3. Bridge BF for an interior primitive coefficient.
4. Bridge BF for a derived linear/factor coefficient.
5. Bridge BF for a one-sided heterogeneity boundary.
6. At least one normal GLMM and one publication-bias active-branch case, using
   a separately fitted nested model where the Savage-Dickey prior condition is
   satisfied.

Use a combined uncertainty rule rather than only a fixed tolerance: compare
log-BF differences against qCMDE/IWMDE MCSE plus bridge MCSE, with a modest hard
maximum such as `0.05-0.10`. Increase row budgets in certification tests instead
of allowing 28-42% discrepancies.

**Decision:** Approved. Strengthen end-to-end numerical certification and tighten
the bridge-oracle acceptance criteria.

### F4. Medium: IWMDE MCSE conditions on weights estimated from the same draws

**Evidence**

- IWMDE estimates conditional-normal, Gamma, beta, or marginal-normal weight
  parameters from posterior values at `R/iwmde-density.R:807-910` and
  `R/iwmde-density.R:968-1340`.
- The same posterior rows are then used to evaluate the importance
  contributions.
- `.iwmde_batch_mcse()` computes batch variability of the final contribution
  matrix at `R/iwmde-density.R:1489-1537`. It treats the fitted weight function
  as fixed and has no term for weight-estimation uncertainty.
- Rows are deterministically thinned, then chain identifiers are discarded at
  `R/iwmde-api.R:2165-2184` and `R/utilities.R:538-602`.

No material undercoverage was demonstrated in the reviewed examples. At larger
budgets, reported IWMDE errors were compatible with observed oracle errors.
This remains a methodological uncertainty in the reported `BF_error`, most
relevant for high-dimensional conditioning, small branches, or fallback
weights.

**Suggested action**

Retain the current feasible BF error computation. Document `BF_error` as a
conditional Monte Carlo error estimate, not a complete standard error, and make
clear that uncertainty from estimating the IWMDE weight function is excluded.
Do not add repeated-chain bootstrap, block bootstrap, cross-fitting, or other
expanded BF error computations.

**Decision:** Approved with scope limitation. Documentation only; no more
extensive BF error computation will be implemented.

### F5. Medium: stored-result source fingerprint samples only ten posterior rows

**Evidence**

- `.iwmde_source_fingerprint()` hashes posterior dimensions, column names, and
  only the first and last five posterior rows at `R/iwmde-plan.R:499-526`.
- A runtime mutation of row 7,500 in a 15,000-row posterior matrix produced an
  identical source fingerprint.
- Stored marginal-mean densities/ordinates are reused by exact request-key
  equality at `R/iwmde-result.R:243-250` and `R/iwmde-result.R:280-340`.
- RoBMA and BayesTools versions are attached after request-key construction at
  `R/iwmde-result.R:194-210`, so version changes do not themselves invalidate a
  result unless the manually maintained schema version changes.

Normal package update operations that change draw count will invalidate the
fingerprint. The risk is same-shaped replacement or mutation of posterior draws,
or failure to bump the estimator schema after an algorithm change.

**Suggested action**

Hash all posterior columns and rows used by the estimator, or attach an
immutable fit-draw fingerprint when the fit is created. Include an explicit
estimator algorithm version in the request key and add a test requiring that it
changes when result semantics change. A package version alone is too broad;
the current ten-row content hash is too weak.

**Decision:** Approved. Strengthen draw fingerprints and estimator schema
versioning.

### F6. Medium: user documentation does not define estimator choice or reliability policy

**Evidence**

- `vignettes/v22-posterior-densities-parameter-tests.Rmd:172-193` describes
  qCMDE as smoother and more stable than KDE, but does not present its defining
  conditional-density average or distinguish it from IWMDE.
- The vignette uses qCMDE throughout and does not explain when a user should
  select IWMDE.
- `density_control` is listed without defaults or tuning guidance in
  `R/hypothesis.R:97-108` and generated `man/hypothesis.Rd:64-75`.
- The vignette does not explain ESS, maximum weight share, quadrature movement,
  active mixture mass, or what `error%(BF)` includes.
- The normal approximation is offered alongside qCMDE/IWMDE, but its failure
  modes at boundaries and tail ordinates are not illustrated.

**Suggested action**

Document:

- the qCMDE and IWMDE estimating equations;
- the Savage-Dickey prior compatibility requirement;
- qCMDE as the preferred robust default where supported, with IWMDE as a
  potentially faster but weight-sensitive alternative;
- all control defaults, diagnostics, warning/failure thresholds, and the
  interpretation of `BF_error`;
- support exclusions: non-known-`V` `brma.mv()` random formulas, semantic
  random-effect quantities, and focal weightfunction coordinates;
- a sensitivity workflow that increases `max_samples` and normalization budget;
- why a global normal approximation is useful for near-normal interiors but is
  not a validation oracle for skewed tails or one-sided boundaries.

Relevant methodological references include [Chen (1994)](https://www.tandfonline.com/doi/abs/10.1080/01621459.1994.10476815)
and the open-access [marginal-density estimator review and application](https://pmc.ncbi.nlm.nih.gov/articles/PMC6602590/).

**Decision:** Approved. Expand the user-facing estimator and reliability
documentation before release.

### F7. Medium: known-`V` qCMDE can be slow even for a very small fit

**Evidence**

- On the six-row known-`V` boundary fit with only 240 posterior draws, qCMDE
  calls took about 23-48 seconds. IWMDE calls took about 2-4 seconds.
- The qCMDE path evaluates multiple pilot/final/validation grids and retains a
  full grid-by-row matrix in `R/iwmde-density.R:245-489`.
- Repeated separate public calls rebuild the context and estimator cache; cache
  reuse is scoped to one orchestration call.

The extra qCMDE work buys robustness and excellent accuracy, but the cost is
surprising for a tiny model and may be prohibitive for larger known-`V` blocks,
GLMM quadrature, several marginal means, or interactive plotting.

**Suggested action**

Add performance benchmarks by likelihood family, parameter type, row count,
and grid size. Profile the known-`V` predictor and normalization paths. Consider
bounded grid chunks, reusable per-fit likelihood state, and public progress
messages for expensive computations. Do not trade away the independent final
validation grid without an equivalent numerical guard.

**Decision:** Approved. Benchmark and profile known-`V` qCMDE without weakening
its numerical validation safeguards.

### F8. Medium: core implementation remains too list-heavy and repetitive

**Evidence**

- The subsystem contains about 13,000 lines and 353 functions across 14
  `R/iwmde-*.R` files.
- `R/iwmde-api.R` alone contains 88 functions and 2,248 lines, mixing public
  validation, diagnostic policy, provenance adapters, row selection, and
  marginal-mean orchestration.
- Long functions include `.iwmde_density_grid()` (247 lines),
  `.iwmde_plan_prepare_contract()` (209),
  `.iwmde_density_iwmde()` (177), and
  `.iwmde_plan_diagnostic_result()` (169).
- The unconstrained and logit conditional-normal weight implementations repeat
  centering, scaling, covariance regularization, conditioning, and variance
  checks at `R/iwmde-density.R:1020-1188` and
  `R/iwmde-density.R:1200-1290`.
- Plans, row states, density results, diagnostics, attributes, and provenance
  are informal lists. Missing fields often become `NULL` and are normalized
  later, making schema drift difficult to detect statically.

The current file split is substantially better than a monolith, and tests cover
many helpers. The remaining complexity makes future likelihood-family changes
or BayesTools coordination unnecessarily risky.

**Suggested action**

Refactor incrementally after behavior is frozen:

1. Split public orchestration, control policy, diagnostics, and provenance out
   of `iwmde-api.R`.
2. Split qCMDE normalization, IWMDE weights, aggregation, and MCSE out of
   `iwmde-density.R`.
3. Extract one conditional-Gaussian regression kernel used by identity and
   transformed targets.
4. Add constructors and validators for plan, row-state, density-result, and
   diagnostic schemas, with required fields and schema versions.
5. Replace repeated list-field copying with narrow conversion functions at
   module boundaries.
6. Keep fast-path/scalar parity tests as characterization tests during the
   refactor.

**Decision:** Approved. Refactor incrementally after behavior and certification
tests are frozen.

### F9. Low: conditional-weight fallback discards the failure reason

**Evidence**

- `.iwmde_chen_try_weight()` catches every conditional-weight error and returns
  the fallback without recording the condition at `R/iwmde-density.R:912-918`.
- The final method label records that a fallback was used, but not why.

The fallback remains mathematically valid if it is a normalized weight, so this
is not a correctness defect. It makes performance and reliability regressions
hard to diagnose, especially when a new posterior column makes the conditional
covariance singular.

**Suggested action**

Return a structured weight result containing `method`, `fallback_from`, and
`fallback_reason`. Surface aggregate fallback counts/reasons in diagnostics and
test expected fallback explicitly.

**Decision:** Approved. Preserve valid fallbacks and expose structured fallback
reasons in diagnostics.

### F10. Low: pilot/final normalization field names are ambiguous

**Evidence**

- qCMDE returns `normalization_integral` for the pilot grid and
  `normalization_final_integral` for the final grid at
  `R/iwmde-density.R:461-462`.
- In the exact Gaussian probe, these were `0.996238` and effectively `1.0`.
  Reading only `normalization_integral` suggests a residual mass error that the
  final estimator does not have.
- IWMDE uses `normalization_integral` for its actual support-grid diagnostic,
  so the same field name has estimator-dependent semantics.

**Suggested action**

Before release, rename qCMDE fields to `pilot_normalization_integral` and
`final_normalization_integral`, and IWMDE to
`support_grid_normalization_integral`. Keep one derived, consistently named
`normalization_relative_error` for policy checks. Bump the result schema.

**Decision:** Approved. Rename the normalization fields and bump the result
schema before release.

## Numerical validation

### Oracle validity

Bridge sampling and Savage-Dickey agree only when the constrained model uses
the full model's conditional prior for all shared parameters. The reviewed
model pairs were inspected for that nesting condition. Bridge MCSEs were
combined on the log marginal-likelihood scale.

The global normal approximation was treated as an independent diagnostic, not
as ground truth. It is informative for an interior, near-normal coefficient and
deliberately expected to fail for a skewed tail or one-sided boundary.

### 1. Exactly solvable known-`V` Gaussian model

Model: six observations, known block `V`, fixed `tau = 0`, and
`mu ~ Normal(0, 0.5)`.

| Quantity | Analytic | Computed | Difference |
|---|---:|---:|---:|
| Posterior mean | 0.1007682 | 0.1058362 sampled | 0.0050680 |
| Posterior SD | 0.0996723 | 0.0969018 sampled | -0.0027705 |
| Log marginal likelihood | 2.1863401 | 2.1923661 bridge | 0.0060260 |
| Density at 0 | 2.4009713 | 2.4010223 qCMDE | 0.0000511 |
| Density at posterior mean | 4.0025406 | 4.0026257 qCMDE | 0.0000851 |

Bridge log-marginal-likelihood MCSE was `0.005797`; its analytic error was about
`1.04` MCSE. qCMDE curve `L1` error was `2.13e-5`; IWMDE curve `L1` error was
`0.00400` with maximum relative MCSE `0.00629`.

### 2. Ordinary random-effects tail ordinate

Target: `mu = 0` in the vignette's `fit_simple` BCG model.

| Method/budget | BF10 | Relative to conditional oracle | Reported relative MCSE |
|---|---:|---:|---:|
| qCMDE, 500 | 54.117 | -11.1% | 13.3% |
| IWMDE, 500 | 50.630 | -16.9% | 13.6% |
| qCMDE, 1,000 | 67.189 | +10.3% | 6.61% |
| IWMDE, 1,000 | 63.529 | +4.32% | 8.68% |
| qCMDE, 5,000 | 60.778 | -0.20% | 3.70% |
| IWMDE, 5,000 | 60.666 | -0.39% | 3.99% |
| qCMDE, 15,000 | 60.897 | -0.006% | 3.43% |
| IWMDE, 15,000 | 61.036 | +0.22% | 3.48% |
| Closed-form conditional-normal average | 60.901 | reference | 3.43% |
| Bridge sampling | 65.646 | +7.79% | log MCSE 0.00219 |
| Global normal approximation | 218.702 | +259% | not applicable |

The closed-form calculation conditions on each sampled `tau` and analytically
computes `p(mu = 0 | tau, y)`, then averages over posterior `tau` draws. Its
agreement with qCMDE directly validates the qCMDE likelihood/replacement path.
The bridge difference is about 2.2 conditional-ordinate MCSEs and is not
evidence of a qCMDE formula error.

### 3. Interior multilevel regression coefficient

Target: `mu_vi = 0` in `konstantopoulos2011_3lvl2`, compared with the nested
model without the moderator.

| Method | Small budget | Large budget | Bridge BF10 | Large-budget error |
|---|---:|---:|---:|---:|
| qCMDE | 0.379237 at 120 | 0.381513 at 5,000 | 0.382425 | -0.24% |
| IWMDE | 0.376709 at 120 | 0.381968 at 5,000 | 0.382425 | -0.12% |
| Global normal | 0.369638 | n/a | 0.382425 | -3.34% |

Combined bridge log-MCSE was `0.01959`. Both likelihood-aware estimators
converged cleanly and materially outperformed the existing `0.25` log-tolerance.

### 4. Known-`V` one-sided heterogeneity boundary

Target: `tau = 0` in `brma.mv_block_mvn`, compared with
`brma.mv_block_mvn_fixed_random_null`.

| Method | BF10 | Relative to bridge |
|---|---:|---:|
| Bridge sampling | 0.367052 | reference |
| qCMDE | 0.371894 | +1.32% |
| IWMDE | 0.362858 | -1.14% |
| Global normal approximation | 1.002432 | +173.1% |

Combined bridge log-MCSE was `0.03416`. The fit had only 240 posterior draws,
so budgets above 500 did not add rows. The result validates one-sided support,
boundary evaluation provenance, known-`V` joint likelihood, and bridge paths.

### 5. Derived PET factor coefficient

Target: `Preregistered[Pre-Registered] = 0` in
`dat.lehmann2018-PETreg`, compared with `dat.lehmann2018-PET`.

| Method/budget | BF10 | Relative to bridge |
|---|---:|---:|
| Bridge sampling | 2.731818 | reference |
| qCMDE, 120 | 3.088477 | +13.1% |
| IWMDE, 120 | 3.018492 | +10.5% |
| qCMDE, 2,000 | 2.696388 | -1.30% |
| IWMDE, 2,000 | 2.694257 | -1.37% |
| Global normal approximation | 2.574806 | -5.75% |

Combined bridge log-MCSE was `0.01249`. This validates linear replacement,
factor aliasing, PET likelihood evaluation, and both estimators. It also shows
why the current 120-row test and `0.35` tolerance should be replaced by a larger
budget and tighter criterion.

### 6. Prior-density stability

For the multilevel coefficient qCMDE check, changing deterministic prior
integration from 100 to 100,000 points changed BF10 from `0.380431` to
`0.380710`; seeds 1-3 produced identical values. The displayed `BF_error` is
therefore dominated by posterior-ordinate error in this case. Documentation
should still define whether prior integration error is included generally.

## Coverage assessment

| Area | Implementation/test status | Independent numerical oracle in this review |
|---|---|---|
| Ordinary normal meta-analysis | Supported; actual-fit tests | Exact conditional-normal and bridge |
| Normal meta-regression | Supported | Bridge and global normal |
| Multilevel `brma()` | Supported | Bridge coefficient |
| Known-`V` `brma.mv()` | Supported | Exact conjugate and boundary bridge |
| Known-`V` random-formula paths | Supported selectively | No separate nested oracle |
| Non-known-`V` `brma.mv()` random formulas | Explicitly unsupported | Guard verified by tests |
| GLMM | Marginal/conditional likelihood paths and actual-fit tests | No independent nested bridge oracle |
| PET/PEESE | Supported | PET factor bridge |
| Selection/weightfunction active branches | Supported for non-omega focal targets | Internal and actual-fit tests; no bridge oracle here |
| Product-space mixtures | Active-mass/branch logic tested | Bridge unavailable by design |
| Marginal means/linear targets | Supported | PET factor plus analytic linear tests |
| One-sided parameters | Gamma-weight/boundary path supported | Known-`V` tau bridge |
| Bounded parameters | Logit-normal/beta path tested internally | No independent end-to-end oracle |
| Focal weightfunction coordinates | Explicitly unsupported | Guarded |
| Semantic random-effect quantities | Explicitly unsupported | Guarded |

The largest remaining numerical coverage gaps are GLMM, bounded focal
parameters, known-`R`/random-formula targets, and a genuine selection-model
active branch compared with a separately fitted nested model.

## Confirmed strengths

These areas should be preserved:

- qCMDE display, requested ordinate, and normalization grids are separated.
- qCMDE uses independent pilot/final/validation normalization grids.
- Failed row normalizers or weights do not silently redistribute their mass;
  the original denominator is retained and row loss is diagnosed.
- Finite baseline eligibility is frozen in the plan and included in provenance.
- Boundary requested and numerical evaluation values are separate.
- Known-`V` uses the joint likelihood target rather than multiplying correlated
  estimate-level likelihoods.
- Product-space active mass and point masses are represented separately.
- Primitive, indexed, affine display, scale-intercept, and linear target
  replacement paths are explicit.
- Predictor, batch, and scalar implementations have parity tests.
- Unsupported random-formula and focal weightfunction-coordinate cases fail
  explicitly rather than silently using an unrelated KDE ordinate.

## Test evidence

Commands used during this review:

```r
devtools::test(filter = "01-", reporter = "llm")
devtools::test(
  filter = "iwmde|hypothesis|bridgesampling",
  reporter = "llm"
)
```

Results:

- Cache-generation suite: 346 passed, 0 failed, 0 warnings, 1 intentional skip.
- IWMDE/hypothesis/bridge suite after cache generation: 2,523 passed, 0 failed,
  0 warnings, 0 skipped.
- Before cache regeneration, the same targeted suite had 426 passes and 50
  cache-related skips; this is why the official `test-01-*` cache workflow was
  run before accepting the fit-backed result.

Independent runtime probes additionally covered all five numerical cases and
prior-density stability described above. Probe scripts were disposable and
were removed after recording results.

## Proposed decision order

1. **Precision policy:** adaptive BF budgets and stricter thresholds, or an
   explicitly approved static alternative.
2. **Diagnostic surface:** stable result attribute plus public accessor.
3. **Certification tests:** exact conditional oracle and tighter bridge matrix.
4. **BF error scope:** retain the current computation and document it as a
   conditional Monte Carlo error; do not add expanded uncertainty estimation.
5. **Documentation:** equations, method choice, controls, diagnostics, and
   limitations.
6. **Performance:** benchmark and profile known-`V` qCMDE.
7. **Provenance:** strengthen draw fingerprint and estimator schema versioning.
8. **Maintainability:** split modules and introduce validated result schemas
   without changing estimator behavior.

No estimator rewrite is recommended. The numerical evidence supports the
current qCMDE/IWMDE mathematics; work should focus on reliability policy,
observability, test power, and reducing change risk.

## Post-implementation review

**Review date:** 2026-07-13
**Scope:** implementation of the approved F1-F10 decisions plus independent
core, downstream, and tests/documentation reviews.

### Outcome summary

All approved changes were implemented. The estimator definitions were retained;
the work changed reliability policy, observability, numerical certification,
GLMM likelihood evaluation, provenance, and module boundaries. Independent
review found no qCMDE/IWMDE formula, Jacobian, support-transformation, or known-
`V` likelihood defect.

One support boundary was narrowed instead of certified with weaker tolerances:
ordinary-GLMM IWMDE density estimation is rejected explicitly. Point-null error
remained above the approved `0.10` hard oracle bound even when conditional MCSE
diagnostics passed, while descriptive curves failed the density-policy MCSE,
ESS, and normalization gates. qCMDE remains the supported ordinary-GLMM route.

### Finding disposition

| ID | Status | Implementation |
|---|---|---|
| F1 | Implemented | Density plots retain a bounded row policy. Point ordinates use deterministic adaptive budgets and stop only when both the requested relative-MCSE target and all BF-grade gates pass, or the hard cap/all eligible rows are exhausted. Warning/rejection gates now cover MCSE, finite terms, ESS, maximum contribution share, row loss, estimator stability, and quadrature stability. |
| F2 | Implemented | Exported `density_diagnostics()` returns a versioned compact schema for accepted results and structured rejected-result conditions, including budgets, reliability metrics, thresholds, provenance, convergence state, and fallback reasons. |
| F3 | Implemented with explicit exclusion | Added exact/analytic normal oracles and nested bridge-sampling oracles for ordinary normal, multilevel, PET, known-`V` heterogeneity, ordinary GLMM qCMDE, and an active selection-model branch. Tests combine estimator and bridge uncertainty and retain a `0.05-0.10` hard bound. Ordinary-GLMM IWMDE density estimation now fails explicitly rather than weakening point-BF or curve-reliability tolerances. |
| F4 | Implemented as approved limitation | Documentation defines `BF_error` as conditional Monte Carlo error. It excludes uncertainty from learning IWMDE weights and prior-ordinate integration. No bootstrap, cross-fitting, or other expanded BF error computation was added. |
| F5 | Implemented | Source provenance hashes the complete estimator-relevant posterior draws and includes separate schema and algorithm versions. Stored marginal-mean results use semantic request provenance before constructing expensive plans. |
| F6 | Implemented | Roxygen and the posterior-density vignette now define both estimators, controls, thresholds, stability metrics, Savage-Dickey requirements, support limits, sensitivity analysis, normal-approximation limits, and the conditional scope of `BF_error`. |
| F7 | Implemented | `tools/profile-iwmde-qcmde.R` validates caches and profiles a row/grid matrix across known-`V` qCMDE, known-`V` IWMDE, ordinary normal, Poisson GLMM, and selection likelihoods. Independent final/validation normalization grids remain intact. |
| F8 | Implemented | Added validated, versioned plan, row-state, density-result, and diagnostic schemas with estimator-specific semantic validation. Conditional-Gaussian weighting uses a shared kernel. Developer diagnostics, context, method/attribute adapters, parameter orchestration, caches, IWMDE weighting, aggregation/MCSE, qCMDE normalization, known-`V` prediction, and marginal-means integration now live in separate modules. Characterization tests protect batch/scalar and fast/general path parity. |
| F9 | Implemented | Conditional-weight failures retain `fallback_from` and a normalized structured reason. Public diagnostics aggregate fallback counts and reasons. |
| F10 | Implemented | qCMDE exposes `pilot_normalization_integral` and `final_normalization_integral`; IWMDE exposes `support_grid_normalization_integral`. Policy reports a separate estimator-specific stability metric, and the schema version was advanced. |

### Additional correctness changes

- Adaptive point-ordinate evaluation cannot terminate merely because the MCSE
  target passes while an improvable ESS, contribution-share, finite-term, or
  stability gate still fails.
- Rejected ordinate diagnostics survive in a structured `iwmde_ordinate_error`,
  so users can inspect why certification failed.
- Explicit qCMDE/IWMDE plot requests fail when the estimator is unsupported or
  rejected; they cannot silently substitute KDE.
- `density_diagnostics()` reports row-loss fractions and their active warning
  and rejection thresholds directly.
- Expected absence of row-level marginal-likelihood auxiliaries is classified
  narrowly; unrelated BayesTools and prior-construction errors propagate.
- Ordinary binomial and Poisson GLMM likelihoods use native two-dimensional
  adaptive Gauss-Hermite quadrature with explicit convergence failures.
- The zero-effect/zero-heterogeneity binomial case uses its exact powered beta-
  binomial integral, including boundary-cell studies.
- Default ordinary-GLMM evaluation requires an untruncated beta nuisance prior
  for binomial models or an untruncated normal nuisance prior for Poisson models.
  Truncated, point, and other priors fail explicitly.
- Cluster-level GLMM likelihood, LOO, and WAIC targets remain unsupported and
  fail explicitly rather than silently using an ordinary-GLMM likelihood.
- Bridge-oracle tests now verify data, likelihood, grouping/weight/selection
  metadata, and shared nuisance-prior nesting before comparing Bayes factors.

### Performance evidence

The validated profiling matrix completed successfully. Representative qCMDE
known-`V` timings were about 1.15 seconds for 40 rows/51 normalization points,
2.11 seconds for 240 rows/51 points, and 0.92 seconds for 240 rows/201 points on
the review machine. Known-`V` heterogeneity took 1.27 seconds with qCMDE and
0.83 seconds with IWMDE; ordinary normal took 1.13 seconds, Poisson GLMM qCMDE
4.53 seconds, and the active-selection case 8.16 seconds.

The Poisson-GLMM profile recorded about 1.7 GB of cumulative allocations, not
peak resident memory. A direct release probe also measured 159.5 seconds for
one targeted `alloc` marginal-mean qCMDE curve at the reliable 500-row budget.
The latter follows the general linear-combination path and is the clearest
remaining interactive-latency target. Neither measurement indicates a
correctness or convergence failure.

### Residual limitations and risks

1. IWMDE's reported error remains conditional on its learned weights, exactly
   as approved in F4.
2. Ordinary-GLMM IWMDE density estimation is intentionally unsupported;
   weakening point-BF or descriptive-curve gates would misrepresent numerical
   certification.
3. Cluster-level GLMM likelihood targets and unsupported GLMM nuisance priors
   require new methodology before support can be claimed.
4. Bounded focal parameters have internal analytic coverage but less independent
   end-to-end bridge coverage than primitive normal and one-sided targets.
5. Likelihood-aware density estimation remains unavailable for non-known-`V`
   random-formula models, semantic random-effect quantities, and focal
   selection omega/eta coordinates requiring joint replacement.
6. Known-`R` random-formula targets have fast/scalar parity coverage but no
   separately fitted nested bridge oracle.
7. GLMM and selection likelihood profiling should be retained in future release
   work because their allocation and runtime costs dominate normal models.
8. General-path qCMDE marginal-mean curves are expensive: the release probe
   required 159.5 seconds for one level at 500 rows. The vignette therefore
   shows targeted precomputation as non-executed syntax instead of adding
   several minutes to every package build. A future optimization should batch
   linear-combination likelihood replacements without changing reliability
   gates or the independent normalization grids.
9. Concurrent source worktrees must not share an explicitly configured
   `ROBMA_TEST_FILES_DIR`. Cache metadata correctly invalidates fits from a
   different RoBMA or BayesTools source fingerprint, but a check that then
   refits can replace the shared files. Final certification therefore used a
   source-specific cache directory. This is a test-workflow coordination
   constraint, not an estimator or package-runtime limitation.

### Verification record

Final verification used the portable fit-cache hashes committed concurrently
with the multivariate work. Fit regeneration and targeted certification used
an isolated cache so the simultaneous multivariate check could not replace
their metadata:

- `devtools::document()` completed and regenerated the exported documentation;
- `devtools::test(filter = "01-", reporter = "llm")` regenerated all 65
  catalogued fits: 354 passed, 0 failed, 0 warnings, 0 skipped;
- all 65 regenerated fits passed metadata validation;
- the complete qCMDE/IWMDE/hypothesis/bridge filter passed 3,030 assertions
  with 0 failures and 0 warnings; its only skip was the intentional
  `test-01-iwmde-oracle-nested.R` refit skip because the isolated cache was
  already valid;
- the final source-tree suite reported 14,358 passes, 32 expected high-Pareto-
  `k` diagnostic warnings, 55 intentional cache/visual/extended-coverage skips,
  and four visual-only snapshot mismatches from the concurrent GLMM/multivariate
  work;
- `devtools::check()` completed in 1 hour 44 minutes with 0 errors, 0 warnings,
  and 0 notes; installed-package tests took 101 minutes and passed, and the
  vignette rebuild passed;
- the posterior-density vignette also rendered directly in 28 seconds after
  explicit reliable plotting budgets replaced the rejected 500-row examples;
- the validated performance matrix completed with every retained case
  accepted;
- independent reviews of estimator/native mathematics, downstream/API
  integration, and tests/documentation found no remaining qCMDE/IWMDE release
  blocker after their findings were resolved.

The four source-tree visual mismatches are `funnel_BMA.glmm`,
`qqnorm_glmm_base`, `qqnorm_glmm_ggplot`, and `qqnorm_BMA.glmm`. Their existing
references were preserved and testthat's `.new.svg` candidates were not
approved automatically. They reflect changed GLMM fit output outside this
review's qCMDE/IWMDE scope and require maintainer visual review.
