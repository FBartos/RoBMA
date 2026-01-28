# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

RoBMA (Robust Bayesian Meta-Analysis) is an R package for estimating ensembles of meta-analytic models using Bayesian model averaging. It combines models with/without effect, heterogeneity, publication bias, and moderation, using Bayes factors for component testing.

- **Backend**: JAGS (via `runjags`/`rjags`) with custom C++ module for specialized distributions
- **Core Dependency**: `BayesTools` (handles priors, plotting, diagnostics, and generic Bayesian infrastructure)
- **Algorithms**: Two estimation approaches: `"bridge"` (bridge sampling, default) and `"ss"` (spike-and-slab approximation)

**System Requirement**: JAGS 4.3.1+ must be installed (https://mcmc-jags.sourceforge.io/)

**Current Status**: Major refactor in progress from `RoBMA` class to new `brma` class. Files marked with `### obsolete - use for reference only ###` should not be used for current implementation. No backwards compatibility required (major version change).

## Supplementary Instructions

Additional detailed instructions are available in `.claude/instructions/`:

- [r-development.md](.claude/instructions/r-development.md) - R coding standards and conventions
- [testing.md](.claude/instructions/testing.md) - Unit testing workflow and cache system
- [vignettes.md](.claude/instructions/vignettes.md) - Vignette writing and pre-computation patterns
- [posterior-samples.md](.claude/instructions/posterior-samples.md) - Effect direction handling and core computation functions

Read these files when working on the corresponding areas.

## Development Commands

```r
devtools::load_all()              # Load during development
devtools::document()              # Update roxygen2 documentation
devtools::test()                  # Run all tests
devtools::test(filter = "topic")  # Run specific topic tests
devtools::check()                 # Full R CMD check
```

**Testing Note**: Long-running model fits are cached. Run `devtools::test(filter = "fit")` first to generate cached fits, then use `devtools::test(filter = "topic")` for faster iteration.

## Code Style

- **Naming**: `snake_case` for functions, arguments, variables
- **Assignment**: Always use `<-`, never `=`
- **Indentation**: 2 spaces, no tabs
- **Booleans**: Always `TRUE`/`FALSE`, never `T`/`F`
- **Package calls**: Use explicit `::` notation (e.g., `BayesTools::check_bool()`), never `library()`

### Vertical Alignment (Required)

Align assignment arrows `<-` in sequences:
```r
is_multilevel     <- .is_multilevel(x)
is_weightfunction <- .is_weightfunction(x)
```

Align function arguments when spanning multiple lines:
```r
result <- my_function(
  first_arg  = value1,
  second_arg = value2
)
```

Leave an empty line after opening brace `{` in function definitions.

## Architecture

### Component Separation
- **BayesTools**: Generic Bayesian infrastructure (input validation via `check_XXX()`, JAGS settings, plotting)
- **RoBMA**: Meta-analysis-specific logic (ensemble construction, model averaging, publication bias)
- If a feature is generic (not meta-analysis specific), it belongs in BayesTools

### Main User Interfaces
- `brma()`, `brma.glmm()`, `bselmodel()`, `bPET()`, `bPEESE()` - Primary single-model fitting functions for simple models (mostly complete by now, closely matching metafor specification, use `brma` class)
- `RoBMA(), RoBMA.glmm()` - Primary ensemble fitting functions (still to be rewritten)

### Key Current Files
- key functionality is implemented only for `brma` class so far
  - specific functionality in `R/brma.XXX.R` files
  - input handling in `R/input-data.R`,  `R/input-priors.R`,  `R/input-object.R`
  - model specification in `R/fit.R`

### Key Obsolete Files
- `R/fit-and-marglik.R` - Core MCMC fitting and marginal likelihood computation
- `R/inference-and-model-averaging.R` - BMA weights, Bayes factors, posterior mixing
- `R/check-input-and-settings.R` - Input validation, `check_setup()` for previewing ensembles
- `R/transformations.R` - Effect size conversions (Cohen's d ↔ Fisher's z ↔ log OR ↔ r)

### JAGS Extension (src/)
Custom C++ JAGS module with weighted distributions for publication bias:
- `src/distributions/` - Custom distributions (DWT1, DWMN1, etc.)
- `src/functions/` - Matrix operations, density functions
- `src/transformations/` - Effect size transformations
- `src/RoBMA.cc` - Module registration

When adding JAGS distributions: implement in `src/distributions/`, register in `src/RoBMA.cc`, update all `Makevars*` files.

## Documentation

- Use roxygen2 with `@title`, `@description`, `@param`, `@return`, `@examples`, `@export`
- Use `\insertCite{key}{RoBMA}` for literature references
- Internal functions (starting with `.`) need brief header comments

## Testing

- Framework: testthat in `tests/testthat/`
- Files numbered `test-XX-topic.R` for execution order
- Wrap computationally intensive tests with `skip_on_cran()`

## CRAN Compliance

- Minimize dependencies; prefer base R or BayesTools utilities
- No tidyverse (maintain low dependency footprint)
- Never write to user directories without permission
- Always use `::` for package functions
