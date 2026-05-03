# R Implementation Review

Read-only senior review of the current `R/` implementation. No code changes or tests were run during the review.

## High Priority

### Native and package loading is fragile

`R/zzz.R` calls `rjags::load.module()` without guarding errors, uses a `getwd()/src` development fallback, manually calls `library.dynam()`, and has no `useDynLib` entry in `NAMESPACE`.

Suggestion: wrap JAGS module loading in `tryCatch()`, avoid cwd-dependent paths, and either add roxygen `@useDynLib RoBMA, .registration = TRUE` or document why manual loading is required. Use a small load smoke test for representative `.Call` symbols.

### Load-time dependencies are inconsistent

`MASS::ginv()` is used in `R/hat_matrix.R`, but `MASS` is not declared. `parallel::detectCores()` is used while `parallel` is only in `Suggests`.

Suggestion: add required packages to `Imports` or remove those dependencies from load-time/runtime paths. Prefer avoiding load-time detection when a lazy default can be computed later.

### GLMM native availability check is incomplete

`.has_native_glmm()` in `R/pdf.R` checks only the binomial native symbol, but Poisson code later calls `RoBMA_glmm_pois_marginal_loglik`.

Suggestion: split native checks by outcome type or require all GLMM native symbols before enabling the native path. Add a regression test for Poisson GLMM fallback.

### Product-space constructors ignore autocompute options

`RoBMA()`, `BMA.norm()`, and `BMA.glmm()` do not honor `RoBMA.options(autocompute.loo/waic)`, while single-model constructors do.

Suggestion: centralize constructor tail logic in an `.autocompute_brma()` helper. Keep marginal likelihood guarded for product-space objects, but make LOO/WAIC behavior consistent.

### Fitted model matrix reconstruction can diverge from fitted design

`.get_model_matrix()` in `R/residuals.R` rebuilds `model.matrix()` from stored data. This may differ from the actual fitted JAGS/BayesTools design for contrasts, scaling, interactions, aliases, and transformed terms.

Suggestion: move this helper to shared fitted-design infrastructure and use the stored fitted design where possible. Add rank/QR guards for aliased or singular moderator designs.

### Outcome helpers live in the wrong module

`.outcome_data_yi()`, `.outcome_data_sei()`, `.outcome_data_vi()`, and related helpers live in `R/residuals.R`, but are used by diagnostics, plots, pdf/cdf, z-curve, and VIF code.

Suggestion: move them to a small shared file such as `R/outcome-helpers.R`. Centralize GLMM effect-size extraction there so plotting and diagnostics do not duplicate outcome policy.

### Input validation has correctness holes

Numeric `subset` accepts non-integer values; partial `ni` can propagate into `NA` prior scales for `measure = "GEN"`; all-zero Poisson counts can create invalid default log-rate priors; `prior_lograte` lacks the validation used for `prior_baserate`.

Suggestion: tighten integerish checks for subset indices, reject incomplete `ni` when needed for UISD, and add clear positive-event validation for Poisson default priors. Validate `prior_lograte` before indexing or transformation.

### Exported `as_draws*` generics shadow posterior generics

RoBMA exports `as_draws*` generics without default methods, which can break expected `posterior::as_draws*` behavior after attaching RoBMA.

Suggestion: either forward default methods to `posterior::as_draws*` or re-export posterior generics where possible. Document the returned draw format for `brma` and `brma_samples`.

## Medium Priority

### Constructor `...` flags are silently unchecked

Internal flags such as `only_data`, `only_priors`, `is_JASP`, and `is_JASP_prefix` are accepted through `...` and interpreted with permissive checks.

Suggestion: add `.validate_constructor_dots()` to type-check known internal flags and reject unsupported values early. Reuse it across all constructors.

### Bridge log posterior duplicates likelihood reconstruction

`R/marglik.R` reconstructs mu, tau, cluster, PET/PEESE, and selection likelihood pieces separately from evaluation and log-likelihood code.

Suggestion: extract shared helpers that convert parameter lists and data into likelihood matrices. This reduces drift between bridge sampling, log-likelihood, prediction, and posterior evaluation.

### Private APIs are used in several places

The code uses `BayesTools:::` internals in fitting/plotting and `bridgesampling:::` internals in bridge sampling.

Suggestion: replace private calls with local helpers, exported APIs, or namespace lookups isolated behind compatibility wrappers. Add comments where no exported API exists.

### Selection mapping has repeated native-wrapper defaults

Alpha, p-hacking, and kernel defaults are repeated across multiple native wrappers in `R/selection-mapping.R`.

Suggestion: extract one helper for selection defaults and native argument preparation. Add assertions that p-cut breaks are sorted and include required endpoints.

### Diagnostics recompute PSIS/LOO repeatedly

`influence()` calls several helpers that each recompute LOO/PSIS-derived quantities and notes.

Suggestion: create `R/diagnostics-helpers.R` with shared PSIS weights, note handling, labels, and weighted moments. Let influence-style functions pass a cached diagnostic context.

### Diagnostics output labels do not match docs

Some diagnostics documentation promises study-label row names, but `dfbetas`, `covratio`, `hatvalues`, and `influence` outputs do not consistently set them.

Suggestion: define one study-label helper and apply it to all diagnostic return objects. Add tests for label propagation.

### Conditional posterior samples lose chain metadata

Conditional predictions subset posterior rows and rebuild samples as a single chain.

Suggestion: either preserve original chain/draw IDs as attributes and expose them in `as_draws_df()`, or document clearly that conditional samples are flattened. Add tests for draw counts and metadata behavior.

### Posterior indicator handling is duplicated

Conditional prediction and model-summary code both extract indicator columns with slightly different logic and limited range checks.

Suggestion: centralize indicator extraction and validation for `<parameter>_indicator`, `bias_indicator`, and null/alternative mappings. Use it in prediction, summaries, and conditional table construction.

### Publication-bias control naming is inconsistent

Wrappers expose `bias_adjusted`, while plotting uses `sampling_bias`, with inverted semantics in places.

Suggestion: decide whether wrappers are low-level prediction aliases or user-facing helpers. If user-facing, expose `sampling_bias` consistently and keep `bias_adjusted` internal.

### GLMM CDF and residual limitations need clearer public docs

GLMM CDF/PIT paths use normal approximations on effect-size scale rather than exact count distributions.

Suggestion: document this limitation in residual, PIT, and plotting help pages. Add tests that GLMM paths return stable approximations and clear warnings where appropriate.

### Z-curve lacks GLMM guard

`as_zcurve()` does not clearly reject GLMM objects, while lower-level z-curve code reads `outcome$yi` and `outcome$sei`.

Suggestion: either reject GLMM inputs early with a clear message or route through centralized outcome helpers that support GLMM effect-size extraction.

### Plot and diagnostic arguments have unused or early-invalidated paths

Examples include unused `show_data` in `plot_weightfunction.brma()`, unused `is_multilevel` in funnel helpers, and residual-only funnel arguments validated before outcome mode is selected.

Suggestion: remove unused arguments, delay validation until mode is known, or document intentionally ignored arguments. Add `as_data` numeric tests for plot internals rather than relying only on snapshots.

## Lower Priority

### Documentation contradicts input behavior

`R/input-data.R` says `mods` is ignored when `yi` is a formula, but code errors. It also documents vector `mods`/`scale`, while parsing rejects vectors.

Suggestion: either implement the documented behavior or update documentation to match the stricter API. Prefer stricter behavior if formula handling is the intended modern interface.

### `regplot()` uses arbitrary dummy SE for bias-adjusted curves

Publication-bias prediction over moderator grids uses `mean(sei)` as a dummy standard error.

Suggestion: expose `sei`/`vi` control or document the convention. Consider plotting multiple reference SE values if bias adjustment depends materially on precision.

### `qqnorm()` envelope is nondeterministic

Envelope simulation uses unseeded RNG, making repeated plots and visual tests unstable.

Suggestion: add a `seed` or RNG-control argument, or separate deterministic data generation for tests.

### `RoBMA.options()` accepts arbitrary values

Options can be set to invalid booleans, cores, scales, or expressions and fail much later.

Suggestion: define an option schema with defaults, validators, and documentation generated from one source. Avoid evaluating arbitrary public option values.

### Summary table construction is repeated

`summary.R` repeats similar `BayesTools::JAGS_estimates_table()` blocks for common, moderator, scale, bias, and conditional summaries.

Suggestion: factor a small internal table builder. This keeps formatting and probability handling consistent.

### Dataset documentation has minor hygiene issues

Dataset docs have typos, raw URLs, and unnecessary `@return` fields.

Suggestion: clean typos, use `@source` and `\url{}` consistently, and regenerate documentation.

### Source comments contain encoding artifacts

Some comments in `R/pdf.R` contain mojibake.

Suggestion: normalize file encoding and fix the comments. This is not runtime critical but improves maintainability.

## Suggested Implementation Order

1. Fix package load, dependency, and native-symbol issues.
2. Tighten input validation for subset, `ni`, Poisson counts, and prior objects.
3. Centralize outcome extraction and fitted model-matrix helpers.
4. Centralize constructor autocompute and internal `...` validation.
5. Reduce private dependency API usage.
6. Add diagnostics/plot helper modules and PSIS caching.
7. Clean docs, argument naming, and plot/data tests.

## Overall Assessment

The core statistical implementation appears broadly coherent: product-space assumptions, `effect_direction`, selected-normal routing, PET/PEESE zeroing, and `bias_indicator` handling are internally consistent.

The main risk is infrastructure drift. Several helpers live in accidental modules, validation is duplicated, private dependency APIs are used, and diagnostics/plots reconstruct fitted state instead of reusing canonical fitted internals.
