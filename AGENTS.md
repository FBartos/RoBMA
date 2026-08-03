# AGENTS.md

Guidance for coding assistants working in RoBMA. Be extremely concise.

## Package Overview

RoBMA provides Bayesian meta-analysis, model averaging, publication-bias
adjustment, meta-regression, multilevel and multivariate models, GLMMs,
prediction, visualization, and diagnostics.

- Backend: JAGS through `runjags`/`rjags`, plus focused native kernels.
- Core dependency: BayesTools for generic Bayesian infrastructure.
- Estimation: JAGS product-space fitting.
- System requirement: JAGS 4.3.1 or newer.

## Numerical and Engineering Guardrails

### Numerical Faithfulness

- Correct inference-changing numerical errors; do not pursue arbitrary
  machine-level perfection.
- Implement the statistical model as specified. Never silently alter user
  inputs or computed values with epsilons, clamping, covariance repair, or
  similar heuristics.
- Invalid inputs must fail clearly. Mathematically valid boundary cases must
  retain their defined behavior.
- Prefer established implementations from base R, BayesTools, JAGS, or standard
  numerical libraries.
- Support statistically meaningful extreme cases. Do not promise correctness
  for every representable binary64 input unless the public API requires it.
- Base numerical warnings on structural or provenance information where
  possible. Do not infer prior densities from posterior samples.
- Use independent references, analytic identities, or human-verified results
  for numerical tests. Do not test an implementation against the same
  unverified calculation.

### Complexity Budget

- Prefer the smallest boring implementation that solves the observed problem.
- Do not add general exact-arithmetic, arbitrary-precision, or binary64
  emulation subsystems.
- Native C/C++ is justified only by JAGS integration, statistical correctness,
  or a substantial measured performance benefit.
- Before adding native code, abstractions, or source files, verify that existing
  R, BayesTools, JAGS, or package infrastructure cannot solve the problem
  cleanly.
- Avoid parallel implementations, duplicated numerical paths, speculative
  extensibility, and abstractions serving one simple call site.
- Generic Bayesian functionality belongs in BayesTools. RoBMA contains
  meta-analysis-specific integration.
- Remove superseded machinery when replacing an implementation, subject to
  task scope and maintainer approval.

### Compatibility and Release Notes

- During development of an unreleased patch version, do not preserve
  compatibility with earlier iterations of that same unreleased code. Replace
  inferior architecture completely; do not add migrations, deprecated aliases,
  schema adapters, or compatibility layers.
- For functionality present in a released version, prefer backward-compatible
  changes.
- If released architecture prevents a clean solution or appears short-sighted,
  ask the maintainer whether to preserve compatibility or make a breaking
  change. Never assume either choice.
- When feature implementation on a branch is complete, increment the package
  patch version and update `NEWS.md` in the same feature change.
- Write `NEWS.md` for the final branch state that will be squashed into the
  future release. Do not describe intermediate revisions of unreleased features.

### Testing Budget

- Add the smallest high-information regression test proving correctness. Do not
  add tests merely to increase coverage.
- The standard unit, integration, and visual suite must target at most 15
  minutes on the reference machine.
- Expensive numerical certification belongs in the certification profile and is
  reserved for feature certification and major releases.
- Limit each independently runnable certification case to one hour. There is no
  total certification limit; release certification may run all cases in
  sequence.
- Prefer focused tests during development. Do not repeatedly run the full suite
  without a concrete reason.
- Preserve human-verified `vdiffr` snapshots. Never replace meaningful visual
  comparisons with superficial render-only assertions.
- Do not weaken, skip, regenerate, or loosen an existing regression expectation
  merely to make a changed implementation pass.

### Engineering Behavior

- Before non-trivial work, state consequential assumptions and success
  criteria. Do not narrate obvious assumptions.
- Use a short visible plan for multi-step work; revise it when evidence changes.
- When requirements, code, tests, or documentation conflict, name the conflict
  and seek or record a decision instead of guessing.
- Push back with concrete correctness, complexity, maintenance, or runtime
  costs. Propose the simpler alternative and follow the maintainer's decision.
- Establish a simple correct implementation or reference before optimizing.
  Preserve behavior while optimizing and measure performance claims.
- Distinguish dead code introduced by the current change from unrelated
  existing code. Remove the former; do not remove the latter without approval.
- State uncertainty, verification performed, and verification omitted. Never
  present incomplete evidence as completion.

### Change Discipline

- Make focused, reviewable commits with informative messages.
- Commit related code, tests, and documentation together. Avoid documentation-
  only progress commits.
- Do not refactor unrelated code while fixing a numerical issue.
- For ambiguous statistical or architectural choices, stop and record the
  issue, impact, alternatives, and recommendation in
  `.agents/instructions-decisions.md`; do not guess.
- Place agent-authored investigation reports and other temporary work artifacts
  under `.agents/tmp/`, not in the package root.
- Before finishing, review the diff for unnecessary files, abstractions,
  compatibility layers, tests, and native code.

## Detailed Instructions

Read only the guide relevant to the files being changed:

- [testing.md](.agents/instructions/testing.md): test profiles, cached fits,
  metafor comparisons, and visual regression.
- [selected-normal.md](.agents/instructions/selected-normal.md): posterior
  direction and selection-kernel routing.
- [multivariate-targets.md](.agents/instructions/multivariate-targets.md):
  `brma.mv()` predictive and covariance target semantics.
- [plotting.md](.agents/instructions/plotting.md): plot data/rendering separation
  and publication-bias display semantics.
- [vignettes.md](.agents/instructions/vignettes.md): vignette structure, caching,
  style, and citations.

Do not load all guides by default. Before changing an instruction, verify every
referenced file and function against the current tree. Put unresolved maintainer
choices in `.agents/instructions-decisions.md`.

## Development Commands

```r
devtools::load_all()
devtools::document()
devtools::test(filter = "topic", reporter = "llm")
devtools::test(reporter = "llm")
devtools::check()
```

Use `Rscript tools/test-profile.R refresh-standard` only when cached fits are
missing or stale, then `Rscript tools/test-profile.R standard` for the ordinary
suite. Use `Rscript tools/test-profile.R certification --list` to select an
expensive numerical certification case.

## R Code Style

- Use `snake_case` names and `<-` assignment.
- Use two-space indentation and no tabs.
- Use `TRUE`/`FALSE`, never `T`/`F`.
- Do not use `|>` or `%>%`; name intermediate results.
- Qualify non-base calls with `::` in package source. Vignettes may attach
  packages for user-facing examples.
- Prefer clear code over forced vectorization; use `vapply()` for type-stable
  atomic iteration and loops when clearer.
- Use existing BayesTools validators and RoBMA access helpers instead of direct
  object-internal access.
- Never call `setwd()` in package code or write to user directories without
  permission.

Align assignment arrows in sequences:

```r
is_multilevel     <- .is_multilevel(x)
is_weightfunction <- .is_weightfunction(x)
```

Align named arguments in multiline calls:

```r
result <- my_function(
  first_arg  = value1,
  second_arg = value2
)
```

Leave an empty line after an opening brace in function definitions.

## Architecture

### Classes and Interfaces

- `brma`: base class.
- `brma.norm`, `brma.glmm`, `brma.mv`: likelihood/model specializations.
- `RoBMA`: product-space model averaging, extending `brma`.
- Wrapper classes prepend their class, for example
  `c("BMA.norm", "RoBMA", "brma")` and
  `c("BMA.glmm", "RoBMA", "brma.glmm", "brma")`.

Primary interfaces are `brma()`, `brma.norm()`, `brma.glmm()`, `brma.mv()`,
`bselmodel()`, `bPET()`, `bPEESE()`, `RoBMA()`, `BMA()`, `BMA.norm()`, and
`BMA.glmm()`.

For objects returned by S3 methods, follow the existing result-object family.
For method-specific results, prefer `<generic>.<input_class>` and define
matching `print.<result_class>` or `summary.<result_class>` methods. Preserve
meaningful underlying classes such as `data.frame`; do not introduce a second
naming convention for an existing family.

Use `component` in post-fit public APIs to distinguish model parts. Normalize
`component = "mods"` and `component = "location"` through shared helpers; use
`component = "bias"` for publication-bias parameters. Retain released
`parameter_mods` and `parameter_scale` plotting arguments through 4.x. Reserve
`type` for output or prediction kind.

### Selection Models

Selection post-processing is driven by `bias_indicator` and
`.selection_context()`, not by direct `omega` inspection. Use
`.extract_use_normal()` for posterior-row routing. Selected-normal PDF, CDF, and
RNG callers must pass the context. PET and PEESE coefficients are zero in
inactive product-space branches, so their offsets may be added for all rows.

Publication-bias weight functions are unsupported for binomial and Poisson
GLMMs.

### Source Ownership

- Input: `R/input-data.R`, `R/input-priors*.R`, `R/input-object.R`.
- Fitting: `R/fit.R` and model constructors.
- Posterior evaluation: `R/evaluate.R`, `R/pdf*.R`, `R/cdf.R`, `R/rng.R`.
- Outcome access: `R/outcome-helpers.R`.
- Selection routing: `R/selection-mapping.R`.
- Native R interface: `R/distributions.R`.

The focused JAGS extension lives in `src/distributions/`; shared selected-normal
kernels live in `src/selnorm/`; R-native registrations live in `src/init.c` and
`src/r-*.cc`; module registration lives in `src/RoBMA.cc`. Update matching
`Makevars*` files when adding native sources.

## Documentation and CRAN

- Use roxygen2 for exported functions and `\insertCite{key}{RoBMA}` for package
  documentation. Vignettes use Pandoc citations instead.
- Keep dependencies minimal; do not add tidyverse dependencies.
- Use `skip_on_cran()` for computationally intensive tests.
- Never regenerate reference output or visual snapshots without maintainer
  review.
