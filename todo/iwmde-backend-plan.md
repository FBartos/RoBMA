# IWMDE Backend Plan

## Goal

Add backend support for importance-weighted marginal density estimation (IWMDE)
following Chen (1994, JASA), collected in StatsVault as:

- `chen1994jasa_importance-weighted`
- DOI: `10.1080/01621459.1994.10476815`
- Local summary: `.StatsVault/resources/chen1994jasa_importance-weighted_summary.md`

The target is marginal posterior density evaluation for selected continuous
parameters using MCMC draws, without relying on kernel density bandwidths.

## Chen Estimator Requirements

For a parameter split `theta = (theta_j, theta_-j)`, evaluate the marginal
density at grid value `theta_star` as:

```text
mean_i w(theta_star | theta_-j_i) *
  p_unnorm(theta_star, theta_-j_i | data) /
  p_unnorm(theta_j_i, theta_-j_i | data)
```

This requires:

- saved posterior draws;
- a way to substitute one or more parameter values into each draw;
- log unnormalized posterior evaluation at original and substituted draws;
- stable log-ratio averaging;
- conditional weighting density `w(theta_star | theta_-j_i)`;
- explicit handling of point masses and model indicators.

## Current Backend Inventory

Existing pieces:

- Pointwise log-likelihood for saved posterior draws:
  - `logLik.brma()` in `R/loo.R`
  - `.log_lik.brma()` / `.log_lik_estimate.brma()` in `R/pdf.R`
- Outcome density kernels:
  - `.outcome_pdf.norm()`
  - `.outcome_pdf.wnorm()`
  - `.outcome_pdf.binom()`
  - `.outcome_pdf.binom_conditional()`
  - `.outcome_pdf.pois()`
  - `.outcome_pdf.pois_conditional()`
- Cluster-level likelihood helpers:
  - `.log_lik_cluster_norm.brma()`
  - `.log_lik_cluster_wnorm.brma()`
  - `.log_lik_cluster_glmm.brma()`
- Bridge-sampling single-draw likelihood evaluator:
  - `.log_posterior()` in `R/marglik.R`
  - helper parameter reconstruction functions in `R/marglik.R`

Conclusion: likelihood machinery exists for LOO and bridge sampling, but there
is no general API for evaluating the unnormalized posterior at substituted
draws.

## Missing Backend

### 1. Draw Extraction

Add an internal helper that converts posterior samples into a stable parameter
table/list representation:

```r
.posterior_draws_parameters.brma(object)
```

Requirements:

- preserve row order and chain metadata when possible;
- support scalar and indexed parameters;
- parse `omega`, `gamma`, `theta`, `pi`, `phi`;
- retain model indicator columns such as `bias_indicator`;
- avoid repeated `coda::as.mcmc(object[["fit"]])` calls where possible.

### 2. Substitution API

Add a helper for replacing selected parameters in one draw or many draws:

```r
.substitute_draw_parameters(parameters, values)
```

Requirements:

- support scalar parameters such as `mu`, `tau`, `rho`, `PET`, `PEESE`;
- support indexed parameters for `omega`, `gamma`, `theta`, `pi`, `phi`;
- reject unsupported substitutions early;
- preserve mixture indicators unless explicitly changed.

### 3. Pointwise Likelihood at Arbitrary Draws

Refactor `.log_posterior()` likelihood code into reusable pieces:

```r
.log_lik_parameters.brma(object, parameters, pointwise = TRUE)
.log_lik_sum_parameters.brma(object, parameters)
```

Requirements:

- reuse existing outcome PDF kernels;
- return `S x K` pointwise matrix when `pointwise = TRUE`;
- return summed log-likelihood when needed for bridge sampling;
- match `.log_posterior()` exactly for a single unchanged draw;
- handle normal, weighted normal, binomial GLMM, Poisson GLMM, weights, PET,
  PEESE, moderators, scale regression, and multilevel effects.

### 4. Log Prior at Arbitrary Draws

Add a backend for prior density evaluation:

```r
.log_prior_parameters.brma(object, parameters)
```

Requirements:

- evaluate priors from `.create_fit_priors()`;
- evaluate formula priors for moderator and scale terms;
- handle constrained priors and transformed parameters consistently with JAGS;
- distinguish continuous prior density from discrete mixture/model weights;
- return `-Inf` for invalid support.

Risk: this is likely the hardest piece because `BayesTools::JAGS_bridgesampling`
currently owns part of prior evaluation.

### 5. Log Unnormalized Posterior

Add:

```r
.log_unnormalized_posterior.brma(object, parameters, pointwise = FALSE)
```

Definition:

```text
log p_unnorm(theta | data) = log_lik(theta) + log_prior(theta)
```

Requirements:

- use log scale only;
- expose summed likelihood for IWMDE ratios;
- optionally attach pointwise likelihood for diagnostics;
- match bridge-sampling posterior evaluation on unchanged draws.

### 6. IWMDE Core

Add an internal estimator:

```r
.iwmde_density.brma(object, parameter, grid, w = NULL, draws = NULL)
```

Core algorithm:

1. Extract posterior draws.
2. Compute denominator log unnormalized posterior for each draw.
3. For each grid value:
   - substitute `parameter = grid[g]` into every draw;
   - compute numerator log unnormalized posterior;
   - compute `log_w`;
   - average with log-sum-exp.
4. Return density values and diagnostics.

Diagnostics:

- effective sample size of importance terms;
- finite ratio rate;
- estimated integral over grid;
- warnings for unstable tails or all `-Inf` substitutions.

### 7. Default Weighting Densities

Start simple:

- unconstrained continuous parameters: normal conditional approximation from
  empirical MCMC covariance;
- positive parameters such as `tau`: log-normal or gamma approximation;
- bounded parameters such as `rho`: beta approximation;
- weightfunction parameters `omega`: defer or require user-provided `w`.

For the first version, support only univariate continuous parameters.

### 8. Mixture and Spike Handling

Rules:

- Do not estimate a continuous density for a pure spike component.
- For mixture priors, condition on compatible active model indicators by default.
- Report spike probability separately when the parameter has point mass at a
  fixed value.
- Do not substitute continuous values into inactive components unless the caller
  explicitly requests product-space density behavior.

This is essential for RoBMA model-averaged posteriors.

## Public API Later

After backend validation, consider:

```r
posterior_density(object, parameter, method = "iwmde", grid = NULL, ...)
```

Return object should include:

- `parameter`;
- `grid`;
- `density`;
- `method`;
- `spike_probability`;
- diagnostics;
- citation metadata for Chen (1994).

## Implementation Order

1. Extract and normalize posterior draw parameters.
2. Refactor bridge likelihood into reusable arbitrary-draw likelihood helpers.
3. Add exact parity tests against existing `logLik.brma()` and `.log_posterior()`.
4. Implement log-prior evaluation or identify BayesTools API needed.
5. Add log unnormalized posterior helper.
6. Implement simple univariate IWMDE for non-mixture `brma.norm`.
7. Add stable log-sum-exp averaging and diagnostics.
8. Extend to PET/PEESE and weightfunction models.
9. Add mixture/spike accounting for `RoBMA`.
10. Add public API and documentation only after internal tests pass.

## Tests

Unit tests:

- unchanged draw likelihood equals existing `logLik.brma()` row sums;
- unchanged draw `.log_unnormalized_posterior.brma()` equals bridge path up to
  known constants;
- invalid substitutions return `-Inf` or clear errors;
- log-sum-exp averaging handles extreme likelihood ratios;
- IWMDE on a simple normal-normal model integrates close to 1;
- IWMDE density at posterior mean agrees with kernel density roughly in simple
  non-mixture models;
- spike models report point mass and skip continuous density.

Regression tests:

- no change in existing LOO/WAIC output;
- no change in bridge-sampling marginal likelihoods;
- weighted normal selection likelihood still matches JAGS tests.

## Open Questions

- Does BayesTools expose a supported prior log-density evaluator for full JAGS
  parameter lists, including formula and mixture priors?
- Should IWMDE operate on product-space draws or on component-conditional draws
  for model-averaged `RoBMA` objects?
- Which parameters should be supported first: `mu`, `tau`, moderator
  coefficients, `rho`, PET/PEESE, or `omega`?
- Should multilevel random effects be treated as nuisance parameters in IWMDE
  or integrated out for density plots?

