# RoBMA

RoBMA owns Bayesian meta-analysis, publication-bias adjustment, model averaging,
meta-regression, multilevel and multivariate models, GLMMs, and their post-fit
interfaces. BayesTools owns the reusable Bayesian and formula infrastructure.

When this checkout is inside BayesToolsVerse, follow the shared
[workspace guidance](../AGENTS.md), [validation policy](../.agents/instructions/validation.md),
[R environment](../.agents/instructions/r-environment.md), and
[public API contracts](../.agents/instructions/public-api.md).
They own workflow, research tools, agent scratch space, and the distinction
between tests, verification, and scenarios. A standalone checkout retains the
package contracts below; the optional parent workspace is not a dependency.
Check [pending workspace decisions](../.agents/decisions.md) when present.

## Contracts to consult

Read the guide relevant to the change and maintain it with the implementation:

- [Model interfaces](.agents/instructions/model-interface.md): BayesTools
  integration, scale completion, semantic names, classes, and public selectors.
- [Predictive targets](.agents/instructions/multivariate-targets.md): prediction,
  conditioning, sampling and random-effect covariance, LOO, and marginal likelihood.
- [Selected normal](.agents/instructions/selected-normal.md): effect direction,
  product-space routing, selection sources, and supported kernel paths.
- [Testing](.agents/instructions/testing.md): ordinary tests, verification
  profiles, fit caches, metafor comparisons, and visual regression.
- [Scenarios](.agents/instructions/scenarios.md): readable maintainer analyses,
  cached fits, output candidates, and timing records.
- [Plotting](.agents/instructions/plotting.md): plot-data/rendering separation
  and publication-bias display semantics.
- [Vignettes](.agents/instructions/vignettes.md): cached examples and citations.

Use BayesTools' compiled formula metadata, parameter maps, and public accessors.
Fix shared infrastructure in BayesTools when it owns the behavior. Keep the
model's statistical target explicit; do not infer it from posterior draws or
substitute a target with matching dimensions.

## Source ownership

- Input: `R/input-data.R`, `R/input-priors*.R`, and `R/input-object.R`.
- Fitting: `R/fit.R` and model constructors.
- Posterior evaluation: `R/evaluate.R`, `R/pdf*.R`, `R/cdf.R`, and `R/rng.R`.
- Outcome access: `R/outcome-helpers.R`.
- Selection routing: `R/selection-mapping.R`.
- Native R interface: `R/distributions.R`.

Use existing BayesTools validators and RoBMA access helpers instead of direct
object-internal access. Qualify non-base calls with `::` in package source;
vignettes may attach packages for user-facing examples. Match nearby R
formatting, including aligned related assignments/arguments and the blank line
after a function's opening brace. Prefer existing imports and do not add
tidyverse dependencies.

## Development and backend

Requires R >= 4.3.0, C++17, and JAGS >= 4.3.1 through `runjags`/`rjags`.
Fitting uses JAGS product-space models and focused native kernels. In the
workspace, use its configured R and private agent library for these commands:

```r
devtools::load_all()
devtools::document()
devtools::test(filter = "topic", reporter = "llm")
devtools::check()
```

```text
Rscript tools/test-profile.R standard
Rscript tools/test-profile.R refresh-standard
Rscript tools/test-profile.R certification --list
Rscript tools/test-scenario.R <name>
```

Start with the affected tests. `standard` uses existing valid caches;
`refresh-standard` populates missing or stale fits. `certification` is the
runner's name for deeper verification; select relevant cases separately from
routine tests. Scenario baselines require maintainer or explicitly delegated
review, as described in the scenario guide.

Native JAGS distributions live in `src/distributions/`; selected-normal kernels
are in `src/selnorm/`. R-native registrations are in `src/init.c` and
`src/r-*.cc`; module registration is in `src/RoBMA.cc`. Keep registrations and
matching `Makevars*` source lists consistent when native sources change.

## Documentation and release

- Document exports with roxygen2 and `\insertCite{key}{RoBMA}`;
  vignettes use Pandoc citations.
- For a completed feature, increment the package patch version and update
  `NEWS.md` for the final feature state.
- Preserve released interfaces, including the compatibility periods documented
  in the model-interface guide.
- Use `skip_on_cran()` for computationally intensive tests. Keep `AGENTS.md`,
  `.agents/`, and development-only material excluded through `.Rbuildignore`.
