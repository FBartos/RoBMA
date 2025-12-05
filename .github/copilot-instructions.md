# GitHub Copilot Instructions for RoBMA

You are an expert R developer specializing in Bayesian statistics, meta-analysis, and package development. You are assisting with the `RoBMA` package.

## Project Context
`RoBMA` (Robust Bayesian Meta-Analysis) is an R package that implements a framework for estimating ensembles of meta-analytic models. It uses Bayesian model averaging to combine competing models (e.g., presence/absence of effect, heterogeneity, publication bias).
- **Backend**: JAGS (via `runjags`/`rjags`).
- **Core Dependency**: `BayesTools` (handles priors, plotting, diagnostics, and generic Bayesian infrastructure).

## Code Style & Conventions
- **Naming**: Use `snake_case` for all function names, arguments, and variables.
- **Assignment**: Use `<-` for assignment
- **Indentation**: Use 2 spaces. No tabs.
- **Boolean**: Always use `TRUE` and `FALSE`, never `T` or `F`.
- **Comments**: Comment complex logic. Use roxygen2 style `#'` for function documentation.

## Documentation (roxygen2)
- All exported functions must be fully documented.
- Required tags: `@title`, `@description`, `@param` (with types and defaults), `@return`, `@examples`, `@export` (if public).
- Use `\code{}` for variable names and code snippets in docs.
- Use `\insertCite{key}{RoBMA}` for literature references.

## Testing
- **Framework**: `testthat`.
- **Location**: `tests/testthat/`.
- **Conventions**:
  - Test files are numbered: `test-XX-topic.R` (e.g., `test-04-fit.R`).
  - Use `skip_on_cran()` for computationally intensive tests (model fitting).
- **Requirements**:
  - Write comprehensive tests for new features.
  - Test edge cases and error conditions.
  - Ensure reproducibility (set seeds if needed, though `BayesTools` often handles this).
  - Verify CRAN compliance (no writing to user directories without permission, etc.).

## Architecture & Best Practices
1.  **BayesTools First**: 
    - Use `BayesTools` for generic Bayesian functionality (priors, mixing, plotting).
    - `BayesTools` handles input validation (e.g., `check_bool`), JAGS settings (e.g., `JAGS_check_and_list_fit_settings`), and posterior mixing.
    - If a feature is generic (not specific to meta-analysis), suggest implementing it in `BayesTools` instead of `RoBMA`.
2.  **S3 Class Structure**:
    - The main object is of class `RoBMA`.
    - Implements standard S3 methods: `summary`, `plot`, `print`, `update`, `predict`.
    - `update.RoBMA` is complex; handles refitting, extending samples, and adding models.
3.  **JAGS Integration**:
    - Models are specified for JAGS.
    - Ensure efficient parameterization in JAGS code.
4.  **CRAN Compliance**:
    - Strictly follow CRAN policies.
    - Minimal dependencies (avoid adding new dependencies unless necessary).
    - No `library()` or `require()` calls in functions; use `::`.

## Key Files & Components
- `R/RoBMA.R`: Main interface function.
- `R/priors.R`: Prior distribution definitions.
- `R/fit-and-marglik.R`: Model fitting and marginal likelihood computation.
- `R/inference-and-model-averaging.R`: BMA logic.
- `R/BiBMA.R` / `R/NoBMA.R`: Specialized model wrappers.

## JAGS Extension (src/)
The package includes a compiled JAGS extension to support custom distributions and functions.
- **Location**: `src/` contains the C++ source code (`.cc`, `.h`) and subdirectories for `distributions`, `functions`, etc.
- **Adding New Distributions**: When adding a new distribution, ensure it is:
  1.  Implemented the exported functions in `src/distributions/`.
  2.  Implemented common function in `src/source/` (if applicable).
  3.  Registered in `src/RoBMA.cc`.
  4.  Added to `OBJECTS` in Makevars files.

## Interaction Guidelines
- When asked to implement a feature, first check if it fits `RoBMA` or `BayesTools`.
- If modifying the instruction file, keep it concise and high-level.
- When generating R code, prefer base R or `BayesTools` utilities. Never suggest `tidyverse` (the package aims for low dependency footprint).
