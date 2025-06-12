This is an R package that must strictly follow CRAN guidelines for code, documentation, and package structure.

Always use snake_case for function and variable names throughout the codebase. Use consistent style with the existing R code.

Place all R scripts in the R/ directory, tests in tests/testthat/, and vignettes in vignettes/.

All exported functions require roxygen2 documentation with clear parameter descriptions, return values, and examples.

Use testthat for all tests and ensure tests cover edge cases, error handling, and input validation.

When deprecating functions or arguments, use lifecycle badges and provide deprecation warnings with clear migration paths.

Use message(), warning(), and stop() for user feedback instead of hardcoded messages.

Ensure all examples and vignettes are reproducible and do not rely on random seeds unless explicitly set.

Avoid using non-CRAN dependencies to maintain CRAN compliance.

Follow test-driven development practices and ensure comprehensive test coverage.

RoBMA is a meta-analysis package that implements Bayesian meta-analytic models using JAGS.
It is designed to be a user-friendly interface for conducting Bayesian meta-analyses,
while providing flexibility for advanced users to customize their analyses.

For interaction with JAGS, posterior samples, prior distributions, plotting, diagnostics, etc., 
import functions from the BayesTools package. The purpose of the BayesTools package is to contain all 
generic functionality that is not meta-analysis specific and that can be reused across different packages.

If a new generic functionality is needed, make note that such a feature should be implemented in BayesTools directly.