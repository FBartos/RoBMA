# AGENTS.md

This file provides guidance to coding assistants when working with code in this repository.

In all interactions be extremely concise and sacrifice grammar for the sake of concision.

## Package Overview

RoBMA (Robust Bayesian Meta-Analysis) provides a framework for Bayesian meta-analysis, including model estimation, prior specification, model comparison, prediction, summaries, visualizations, and diagnostics. 
The package fits single and model-averaged meta-analytic, meta-regression, multilevel, publication bias adjusted, and generalized linear mixed models, with optional components for effect presence, heterogeneity, moderation, and publication bias.
BMA and RoBMA use Bayesian model averaging to combine competing models, weight inference by posterior model probabilities, and test model components using Bayes factors.
Users can specify prior distributions for effect sizes, heterogeneity, publication bias, including selection models and PET-PEESE, and moderators.

- **Backend**: JAGS (via `runjags`/`rjags`) with custom C++ module for specialized distributions
- **Core Dependency**: `BayesTools` (handles priors, plotting, diagnostics, and generic Bayesian infrastructure)
- **Estimation**: JAGS product-space fitting.

**System Requirement**: JAGS 4.3.1+ must be installed (https://mcmc-jags.sourceforge.io/)

## Supplementary Instructions

Additional detailed instructions are available in `.agents/instructions/`:

- [r-development.md](.agents/instructions/r-development.md) - R coding standards and conventions
- [testing.md](.agents/instructions/testing.md) - Unit testing workflow and cache system
- [vignettes.md](.agents/instructions/vignettes.md) - Vignette writing and pre-computation patterns
- [posterior-samples.md](.agents/instructions/posterior-samples.md) - Effect direction handling and core computation functions
- [helper-functions.md](.agents/instructions/helper-functions.md) - Internal helper functions for data access
- [metafor-comparison-tests.md](.agents/instructions/metafor-comparison-tests.md) - Testing against metafor reference
- [test-file-template.md](.agents/instructions/test-file-template.md) - Template for test-02-*.R files
- [s3-class-naming.md](.agents/instructions/s3-class-naming.md) - S3 class naming conventions
- [use-normal-subdispatch.md](.agents/instructions/use-normal-subdispatch.md) - Performance optimization for mixed weighted/normal samples
- [bias-indicator-extraction.md](.agents/instructions/bias-indicator-extraction.md) - RoBMA bias model identification from posterior
- [predict-newdata-workaround.md](.agents/instructions/predict-newdata-workaround.md) - Outcome requirements for predict.brma() with newdata
- [sampling-bias-parameter.md](.agents/instructions/sampling-bias-parameter.md) - sampling_bias parameter pattern for RoBMA/bPET/bPEESE/bselmodel

Read these files when working on the corresponding areas.

## Development Commands

```r
devtools::load_all()              # Load during development
devtools::document()              # Update roxygen2 documentation
devtools::test(reporter = "llm")                  # Run all tests
devtools::test(filter = "topic", reporter = "llm")  # Run specific topic tests
devtools::check()                 # Full R CMD check
```

**Testing Note**: Always use testthat LLM reporting for unit tests (`reporter = "llm"`; or set `AGENT=1` if a wrapper cannot pass `reporter`). Long-running model fits are cached. 
Run `devtools::test(filter = "fit", reporter = "llm")` first to generate cached fits, then use `devtools::test(filter = "topic", reporter = "llm")` for faster iteration.

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

### Class Hierarchy
- **brma**: Base class for single-model fits (`brma()`, `bselmodel()`, `bPET()`, `bPEESE()`)
- **RoBMA**: Extends brma for ensemble model averaging: `class(x) <- c("RoBMA", "brma")`
- S3 dispatch gives RoBMA all brma methods automatically; override only what differs

### JAGS Behavior for Publication Bias Models

**Critical insight**: JAGS sets `omega = 1` for non-weightfunction posterior samples.

This means:
- When `is_weightfunction = TRUE`, always use `.outcome_pdf.wnorm()` etc.
- Weighted normal with omega = 1 mathematically equals normal distribution
- No per-sample dispatch needed in most cases

**bias_indicator column**: For RoBMA (mixture priors), posterior samples contain `bias_indicator` column (1-indexed) tracking which bias model generated each sample. Use `.extract_use_normal()` for robust detection.

**PET/PEESE coefficients**: JAGS sets PET = 0 and PEESE = 0 for non-PET/PEESE samples. Adding `PET * sei` for all samples is correct (0 contributes nothing).

### GLMM Limitations

Weightfunction publication bias priors are **not supported** for GLMM outcomes:
- Binomial (`brma.glmm()` with OR) - no weightfunction
- Poisson (`brma.glmm()` with IRR) - no weightfunction

### Component Separation
- **BayesTools**: Generic Bayesian infrastructure (input validation via `check_XXX()`, JAGS settings, plotting)
- **RoBMA**: Meta-analysis-specific logic (ensemble construction, model averaging, publication bias)
- If a feature is generic (not meta-analysis specific), it belongs in BayesTools

### Main User Interfaces
- `brma()`, `brma.norm()`, `brma.glmm()`, `bselmodel()`, `bPET()`, `bPEESE()` - Primary single-model fitting functions for simple models (mostly complete by now, closely matching metafor specification, use `brma` class)
- `RoBMA()` - Primary ensemble fitting function with publication bias adjustment
- `BMA()`, `BMA.norm()` - Model-averaging without publication bias (wrapper, omits `model_type`)
- `BMA.glmm()` - Model-averaging for GLMM (binomial/Poisson) without publication bias

### Class Hierarchy (Extended)

Classes are layered for S3 method dispatch:
- **brma**: Base class for all Bayesian meta-analysis
- **brma.norm** / **brma.glmm**: Likelihood type (normal vs GLMM)
- **RoBMA**: Indicates mixture priors (model-averaging via product space)
- **Wrapper classes prepend**: `c("BMA.norm", "RoBMA", "brma")`, `c("BMA.glmm", "RoBMA", "brma.glmm", "brma")`

### Key Current Files
- input handling in `R/input-data.R`,  `R/input-priors.R`,  `R/input-object.R`
- model specification in `R/fit.R`
- posterior evaluation in `R/evaluate.R`, `R/pdf.R`, `R/cdf.R`, `R/rng.R`
- weighted distributions in `R/distributions.R`


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
- Reporter: always use `LlmReporter` for unit tests (`reporter = "llm"` or `AGENT=1`)
- Files numbered `test-XX-topic.R` for execution order
- Wrap computationally intensive tests with `skip_on_cran()`

## CRAN Compliance

- Minimize dependencies; prefer base R or BayesTools utilities
- No tidyverse (maintain low dependency footprint)
- Never write to user directories without permission
- Always use `::` for package functions

