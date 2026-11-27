# Fitting specification

The `brma` family of functions uses the following arguments to specify
the MCMC sampling and fitting settings.

## Arguments

- sample:

  numeric. Number of MCMC samples to save per chain; each chain runs
  `sample * thin` sampling iterations. Defaults to `5000`.

- burnin:

  numeric. Length of burn-in in samples per chain at the thinning
  `thin`, like `sample`: each chain runs `burnin * thin` burn-in
  iterations. Defaults to `2000`. With at least two chains, the chains
  are checked for agreement over the last `min(burnin_extend, burnin)`
  samples of burn-in, taken every `thin` iterations (at least 100
  samples), before sampling starts, and with `autofit` burn-in is
  extended until they agree (see `burnin_max_Rhat` in
  [`set_autofit_control()`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA_control.md)
  for when the check runs). With `burnin` below 100 the check does not
  run, and with `burnin` at most `burnin_extend` it checks the whole
  burn-in. The fitted JAGS object (`object$fit$burnin`) reports the
  burn-in that was run, in samples; the fitting settings
  (`object$fit_control$burnin`) keep the requested value.

- adapt:

  numeric. Length of the adaptation phase in samples per chain at the
  thinning `thin`: where a sampler adapts, each chain runs
  `adapt * thin` adaptation iterations. Defaults to `500`.

- chains:

  numeric. Number of MCMC chains. Defaults to `3`.

- thin:

  numeric. Thinning interval: one sample is retained every `thin`
  iterations. Defaults to `1`. `adapt`, `burnin` and `sample` count
  samples at this thinning, so a fit requests
  `(adapt + burnin + sample) * thin` iterations per chain before any
  automatic extensions of burn-in or sampling (without the adaptation
  iterations when no sampler adapts); see `thin` in
  [`BayesTools::JAGS_fit()`](https://fbartos.github.io/BayesTools/reference/JAGS_fit.html)
  for the largest settings.

- parallel:

  logical. Whether to run MCMC chains in parallel. Defaults to `FALSE`.

- autofit:

  logical. Whether to automatically extend the MCMC chains if
  convergence is not met. Defaults to `FALSE`.

- autofit_control:

  list of autofit control settings. See
  [`set_autofit_control()`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA_control.md)
  for details.

- convergence_checks:

  list of convergence check settings. See
  [`set_convergence_checks()`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA_control.md)
  for details.

- seed:

  numeric. Random seed for reproducibility. Defaults to `NULL`.

- silent:

  logical. Whether to suppress output. Constructors with no explicit
  default use `RoBMA.get_option("silent")` when `silent` is omitted.
  Model-averaging wrappers default to `TRUE` unless explicitly changed.

- ...:

  additional advanced arguments. Fitting functions reject unused
  arguments; currently recognized internal arguments include
  `only_data`, `only_priors`, `is_JASP`, and `is_JASP_prefix`.

## See also

[`brma`](https://fbartos.github.io/RoBMA/dev/reference/brma.md),
[`set_autofit_control`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA_control.md),
[`set_convergence_checks`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA_control.md)
