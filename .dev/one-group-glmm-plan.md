# One-Group GLMM Support Plan

## Goal

Add `metafor::rma.glmm()`-style one-group GLMM support to `brma.glmm()` and
`BMA.glmm()`:

- `measure = "PLO"` for one-group binomial proportions.
- `measure = "IRLN"` for one-group Poisson incidence rates.

The new models should work with the existing brma infrastructure: priors,
moderators, scale regression, multilevel models, prediction, residuals,
LOO/WAIC, marginal likelihood, and BMA mixture priors.

## User-Facing API

Add data arguments to `brma.glmm()` and `BMA.glmm()`:

- `xi`: event count.
- `ni`: binomial total for `PLO`.
- `mi`: optional non-event count for `PLO`; allow `ni = xi + mi`.
- `ti`: person-time/exposure for `IRLN`.

Supported calls:

```r
brma.glmm(measure = "PLO", xi = events, ni = total, data = dat)
brma.glmm(measure = "PLO", xi = events, mi = nonevents, data = dat)
brma.glmm(measure = "IRLN", xi = events, ti = time, data = dat)

BMA.glmm(measure = "PLO", xi = events, ni = total, data = dat)
BMA.glmm(measure = "IRLN", xi = events, ti = time, data = dat)
```

## Model Definitions

### PLO

One observed count per estimate:

```text
theta_i ~ Normal(0, 1)
eta_i = mu_i + theta_i * tau_i
logit(p_i) = eta_i
xi_i ~ Binomial(ni_i, p_i)
```

With multilevel models, `mu_i` already includes the cluster-level contribution.

### IRLN

One observed count per estimate with log exposure offset:

```text
theta_i ~ Normal(0, 1)
eta_i = mu_i + theta_i * tau_i
log(lambda_i) = eta_i + log(ti_i)
xi_i ~ Poisson(lambda_i)
```

Here `mu` is the log incidence rate per unit exposure.

## Code Changes

### Data Input

Files:

- `R/brma.glmm.R`
- `R/BMA.glmm.R`
- `R/input-data.R`
- `R/input-priors.R`

Tasks:

- Add `xi`, `mi`, `ti`, and reuse `ni` in GLMM signatures.
- Extend `.check_measure()` with `"PLO"` and `"IRLN"`.
- Extend `.check_and_list_data()` GLMM dispatch:
  - `"OR"` -> existing two-group binomial.
  - `"IRR"` -> existing two-group Poisson.
  - `"PLO"` -> new one-group binomial parser.
  - `"IRLN"` -> new one-group Poisson parser.
- Add outcome types:
  - `"prop"` for `PLO`.
  - `"rate"` for `IRLN`.
- Add parser `.check_and_list_data.outcome.prop()`:
  - require `xi` and either `ni` or `mi`.
  - validate integer counts, nonnegative.
  - if `mi` supplied, compute `ni <- xi + mi`.
  - check `xi <= ni`.
  - store `xi`, `ni`, optional `cluster`, `slab`, `weights`.
- Add parser `.check_and_list_data.outcome.rate()`:
  - require `xi` and `ti`.
  - validate `xi >= 0`, `ti > 0` unless `skip_validation = TRUE`.
  - store `xi`, `ti`, optional `cluster`, `slab`, `weights`.
- Update `.prepare_newdata()` to require:
  - `xi`, `ni` for `"prop"`.
  - `xi`, `ti` for `"rate"`.

### Priors

File:

- `R/input-priors.R`

Tasks:

- Keep `theta` for all GLMM outcome types.
- Add `pi` only for `"bin"` (`OR`).
- Add `phi` only for `"pois"` (`IRR`).
- Do not add nuisance baseline priors for `"prop"` or `"rate"`.
- Extend known UISD logic:
  - `PLO`: default to `sqrt(4)` unless we decide on data-based scaling.
  - `IRLN`: use data-based rate scaling `sqrt(sum(ti) / sum(xi))`, with clear handling for all-zero events.
- Review informed prior behavior:
  - no obvious medicine informed priors for `PLO`/`IRLN`.
  - either map cautiously or reject `prior_informed_field` for these measures.
- Document BMA null interpretation:
  - `PLO`: `mu = 0` means `p = 0.5`.
  - `IRLN`: `mu = 0` means rate = 1 per exposure unit.

### JAGS Syntax And Fit Data

File:

- `R/fit.R`

Tasks:

- `.create_fit_priors()`:
  - set levels for `theta` for `"prop"` and `"rate"`.
  - do not set levels for `pi`/`phi`.
- `.create_fit_data()`:
  - for `"prop"` include `xi`, `ni`.
  - for `"rate"` include `xi`, `ti`.
- `.create_model_syntax()`:
  - keep common `mu_estimate`, scale, multilevel, and `theta[i] * tau_estimate` logic.
  - branch likelihoods:
    - `"bin"`: existing two-group binomial.
    - `"pois"`: existing two-group Poisson.
    - `"prop"`: `logit(p[i]) = mu_estimate + theta[i] * tau_estimate`; `xi[i] ~ dbinom(p[i], ni[i])`.
    - `"rate"`: `log(lambda[i]) = mu_estimate + theta[i] * tau_estimate + log(ti[i])`; `xi[i] ~ dpois(lambda[i])`.
- Check current weighted GLMM TODO behavior before allowing `weights` for new outcome types.

### Evaluation, Prediction, And Likelihood

Files:

- `R/evaluate.R`
- `R/rng.R`
- `R/pdf.R`
- `R/cdf.R`
- `R/marglik.R`
- `R/predict.R`

Tasks:

- Reuse `.evaluate.brma.true_effects.glmm()` and `.evaluate.brma.theta.glmm()` for `"prop"` and `"rate"`.
- Add RNG helpers:
  - `.outcome_rng.prop(mu_samples, ni)`.
  - `.outcome_rng.rate(mu_samples, ti)`.
- In `predict.brma(type = "response")`:
  - for `"prop"`, sample one binomial count per estimate.
  - for `"rate"`, sample one Poisson count per estimate.
  - if `as_measure = TRUE`, return logit proportion for `PLO` and log rate for `IRLN`.
- Add marginal log-likelihood helpers:
  - `.outcome_pdf.prop()` integrates over `theta` only.
  - `.outcome_pdf.rate()` integrates over `theta` only.
- Add conditional log-likelihood helpers for bridge sampling:
  - `.outcome_pdf.prop_conditional()`.
  - `.outcome_pdf.rate_conditional()`.
- Extend cluster-unit likelihood dispatch for `"prop"` and `"rate"`.
- Extend CDF helpers using normal approximations on transformed scales.

### Approximate Effect Sizes And Hashing

Files:

- `R/residuals.R`
- `R/unit_level.R`
- `R/regplot.R`
- `R/wrappers.R`

Tasks:

- Extend `.outcome_data_yi()`:
  - `"prop"`: logit proportion with continuity adjustment.
  - `"rate"`: log incidence rate with continuity adjustment.
- Extend `.outcome_data_sei()`:
  - `"prop"`: approximate SE for logit proportion.
  - `"rate"`: approximate SE for log rate.
- Extend `.get_outcome_hash()`:
  - `"prop"`: hash `xi`, `ni`.
  - `"rate"`: hash `xi`, `ti`.
- Update `regplot()` dummy newdata:
  - `"prop"`: `xi = 0`, `ni = 0`.
  - `"rate"`: `xi = 0`, `ti = 0`.
- Confirm `nobs.brma()` still works through `.outcome_data_yi()`.

## Documentation

Files:

- `R/brma.glmm.R`
- `R/BMA.glmm.R`
- `R/input-data.R`
- generated `.Rd` files after `devtools::document()`

Tasks:

- Document one-group binomial and Poisson measures.
- Mention that `model` argument from `metafor::rma.glmm()` is not used.
- Explain parameter scale:
  - `PLO`: pooled logit proportion.
  - `IRLN`: pooled log incidence rate.
- Explain BMA null priors for one-group scales.

## Tests

Files:

- `tests/testthat/test-00-input-data-glmm.R`
- `tests/testthat/test-00-input-data-predict.R`
- `tests/testthat/test-01-brma.glmm.R`
- `tests/testthat/test-10-validation.R`

Tasks:

- Input tests:
  - `PLO` accepts `xi, ni`.
  - `PLO` accepts `xi, mi` and computes `ni`.
  - invalid `xi > ni` errors.
  - `IRLN` accepts `xi, ti`.
  - invalid `ti <= 0` errors unless `skip_validation = TRUE`.
  - subset, NA dropping, weights, slab, cluster.
- Fit smoke tests:
  - simple `brma.glmm(measure = "PLO")`.
  - simple `brma.glmm(measure = "IRLN")`.
  - `BMA.glmm()` for both measures with small sample settings.
- Validation tests against metafor:
  - `metafor::rma.glmm(measure = "PLO", xi = xi, ni = ni)`.
  - `metafor::rma.glmm(measure = "IRLN", xi = xi, ti = ti)`.
- Prediction tests:
  - original-data and newdata prediction.
  - `as_measure = TRUE/FALSE`.
- LOO/WAIC/marglik smoke tests if cached fits make this practical.

## Open Decisions

- Whether `PLO` prior UISD should be fixed at `sqrt(4)` or estimated from binomial information.
- How to handle all-zero events for `IRLN` prior scaling.
- Whether to reject `prior_informed_field` for `PLO`/`IRLN`.
- Whether data `weights` should remain unsupported for GLMM likelihoods until weighted binomial/Poisson JAGS distributions are implemented.
