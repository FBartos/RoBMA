# GitHub Copilot Instructions for RoBMA

You are an expert R developer specializing in Bayesian statistics, meta-analysis, and package development. You are assisting with the `RoBMA` package.

## Project Context
`RoBMA` (Robust Bayesian Meta-Analysis) is an R package that implements a framework for estimating ensembles of meta-analytic models. It uses Bayesian model averaging to combine competing models (e.g., presence/absence of effect, heterogeneity, publication bias).
- **Backend**: JAGS (via `runjags`/`rjags`) with custom C++ module for specialized distributions
- **Core Dependency**: `BayesTools` (handles priors, plotting, diagnostics, and generic Bayesian infrastructure)
- **Algorithms**: Two estimation approaches: `"bridge"` (bridge sampling, default) and `"ss"` (spike-and-slab approximation)

## Code Style & Conventions
- **Naming**: Use `snake_case` for all function names, arguments, and variables
- **Assignment**: Use `<-` for assignment (never `=`)
- **Indentation**: Use 2 spaces. No tabs
- **Boolean**: Always use `TRUE` and `FALSE`, never `T` or `F`
- **Comments**: Document complex logic with inline comments; use roxygen2 `#'` for function headers
- **Package calls**: Use explicit `::` notation (e.g., `BayesTools::check_bool()`), never `library()` or `require()`

## Documentation (roxygen2)
- All exported functions must be fully documented
- Required tags: `@title`, `@description`, `@param` (with types and defaults), `@return`, `@examples`, `@export` (if public)
- Use `\code{}` for variable names and code in docs; `\insertCite{key}{RoBMA}` for literature references
- Internal functions (starting with `.`) should have brief header comments explaining purpose and key concepts

## Testing
- **Framework**: `testthat` in `tests/testthat/`
- **Naming**: Test files are numbered `test-XX-topic.R` (e.g., `test-04-fit.R`, `test-05-methods.R`)
- **Fast testing**: When developing, run `devtools::test(filter = "topic")` to focus on specific areas, model fits are cached to speed up tests (`devtools::test(filter = "fit")` must be run first to generate cached fits)
- **CRAN tests**: Use `skip_on_cran()` for computationally intensive model fitting tests
- **Requirements**: Test edge cases, error conditions, reproducibility; ensure CRAN compliance

## Architecture & Model Fitting Pipeline

### Component Separation (BayesTools vs RoBMA)
- **BayesTools**: Generic Bayesian infrastructure (general input validation via `check_XXX()`, JAGS settings via `JAGS_check_and_list_fit_settings()`, posterior mixing, plotting)
- **RoBMA**: Meta-analysis-specific logic (ensemble construction, model averaging, publication bias adjustment)
- **Rule**: If a feature is generic (not meta-analysis specific), suggest implementing it in `BayesTools` instead

### Model Fitting Flow (`R/fit-and-marglik.R`)
1. **Different Algorithms**: `algorithm = "bridge"` fits individual models and computes marginal likelihoods via bridge sampling; `algorithm = "ss"` fits a single spike-and-slab model
2. **Model Type Detection**: Helpers `.is_model_constant()`, `.is_model_regression()`, `.is_model_multivariate()` determine model characteristics
3. **Data Ordering**: Multivariate models require special data ordering via `.order_data.mv()`
4. **JAGS Model Generation**: Model syntax created via `.generate_model_syntax()` or `.generate_model_syntax.mv()`
5. **MCMC Sampling**: Uses `runjags::run.jags()` with automatic convergence checking (if `autofit = TRUE`)
6. **Marginal Likelihood Calculation**: Bridge sampling via `bridgesampling::bridge_sampler()` (for `algorithm = "bridge"`)
7. **Convergence**: Automatic refitting until criteria satisfied (see `.balance_component_probability()`)

### Model Structure & Classes
- **Base class**: `RoBMA` (S3 object containing `models` list, `data`, `fit_control`, etc.)
- **Subclasses**: `RoBMA.reg` (meta-regression), `BiBMA` (binomial-normal for binary outcomes)
- **Class hierarchy**: `c("BiBMA", "RoBMA")` or `c("RoBMA", "RoBMA.reg")` for method dispatch
- **Key methods**: `summary.RoBMA()`, `plot.RoBMA()`, `update.RoBMA()`, `predict.RoBMA()`

### Model Averaging & Inference (`R/inference-and-model-averaging.R`)
- **Component Structure**: Each model characterized by which components are null vs alternative hypothesis
- **Averaged vs Conditional Estimates**: Separate functions for overall averaged estimates and conditional estimates given specific components
- **Regression Terms**: Handled separately via `predictors_test` and `terms_test`

## Key Files & Their Roles
- `R/RoBMA.R`: Main user interface; ensemble specification via prior lists
- `R/RoBMA-reg.R`: Meta-regression wrapper (uses formula interface)
- `R/BiBMA.R`, `R/BiBMA-reg.R`: Binomial-normal models for binary outcomes
- `R/NoBMA.R`: Publication bias unadjusted models
- `R/priors.R`: Re-exports from `BayesTools` (e.g., `prior()`, `prior_weightfunction()`, `prior_PET()`)
- `R/fit-and-marglik.R`: MCMC fitting and marginal likelihood computation (1745 lines, core logic)
- `R/inference-and-model-averaging.R`: BMA weights, Bayes factors, posterior mixing
- `R/check-input-and-settings.R`: Input validation and `check_setup()` for previewing ensembles
- `R/check-priors-and-models.R`: Prior specification validation
- `R/summary.R`, `R/summary-effect.R`, `R/summary-heterogeneity.R`: Result summaries and statistics
- `R/plots.R`: Plotting methods (forest, funnel, model weights)
- `R/zcurve.R`: Meta-analytic z-curve diagnostics for publication bias
- `R/transformations.R`: Effect size conversions (Cohen's d ↔ Fisher's z ↔ log OR ↔ r)

## JAGS Extension (src/)
The package includes a compiled JAGS module with custom distributions for meta-analysis.

### Structure
- `src/RoBMA.cc`: Module registration (registers all distributions and functions)
- `src/distributions/`: Custom JAGS distributions (e.g., `DWT1.cc` for weighted t-distribution, `DWMN1.cc` for weighted multivariate normal)
- `src/functions/`: Custom JAGS functions (e.g., `mnorm.cc`, `wmnorm.cc`)
- `src/transformations/`: Effect size transformations (d, r, z, logOR, omega)
- `src/Makevars.win`, `src/Makevars.ucrt`, `src/Makevars.in`: Platform-specific build configs
- `configure`, `configure.ac`, `configure.win`: Detect JAGS installation

### Adding New JAGS Distributions
1. Implement `.cc` and `.h` files in `src/distributions/` (follow existing patterns like `DWT1.cc`)
2. Add common functions to `src/source/` if needed
3. Register in `src/RoBMA.cc` via `insert(new YourDistribution);`
4. Add to `OBJECTS` in all `Makevars*` files

### Build Notes
- Requires JAGS ≥ 4.3.1 installed
- Windows: Uses `PKG_CXXFLAGS = -D_GLIBCXX_USE_CXX11_ABI=0` to match JAGS 4.x ABI
- On Windows, JAGS is detected via `JAGS_ROOT` (defaults to `/c/progra~1/JAGS/JAGS-4.3.1`)

## Developer Workflows

### Building & Checking
```r
# Standard R CMD check workflow
devtools::load_all()           # Load during development
devtools::document()           # Update documentation
devtools::test()               # Run tests
devtools::check()              # Full package check
```

### Common Development Tasks
- **Preview ensemble**: Use `check_setup()` to see model structure before fitting
- **Debugging fits**: Check `.is_model_constant()`, `.is_model_regression()`, `.is_model_multivariate()` helpers
- **Algorithm differences**: `algorithm = "bridge"` (individual models + bridge sampling) vs `algorithm = "ss"` (single spike-and-slab model)
- **Convergence issues**: `.balance_component_probability()` handles failed models automatically

### Vignettes (vignettes/)
- `Tutorial.Rmd`: Main introduction to RoBMA-PSMA
- `ReproducingBMA.Rmd`: Classic Bayesian model-averaged meta-analysis
- `MetaRegression.Rmd`: `RoBMA.reg()` with moderators
- `HierarchicalRoBMA.Rmd`: Multilevel RoBMA
- `HierarchicalRoBMARegression.Rmd`: Multilevel RoBMA with moderators
- `HierarchicalBMA.Rmd`: Multilevel models via `study_ids`
- `MedicineBMA.Rmd`, `MedicineBiBMA.Rmd`: Informed priors for medical meta-analysis
- `CustomEnsembles.Rmd`: Advanced ensemble customization
- `FastRoBMA.Rmd`: Spike-and-slab algorithm (`algorithm = "ss"`)
- `ZCurveDiagnostics.Rmd`: Publication bias diagnostics

## CRAN Compliance
- **Dependencies**: Minimize new dependencies; prefer base R or `BayesTools` utilities
- **No tidyverse**: Package maintains low dependency footprint
- **No side effects**: Never write to user directories without permission
- **Namespace**: Always use `::` for package functions (no `library()` in code)
- **Tests**: Wrap long-running tests with `skip_on_cran()`

## Interaction Guidelines
- When implementing features, check if they belong in `BayesTools` (generic) or `RoBMA` (meta-analysis-specific)
- For effect size code, check `R/transformations.R` first (conversions already implemented)
- When modifying JAGS distributions, update all relevant Makevars files
- Keep this instruction file concise and focused on RoBMA-specific patterns
