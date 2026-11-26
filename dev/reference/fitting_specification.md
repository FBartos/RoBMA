# Fitting specification

The `brma` family of functions uses the following arguments to specify
the MCMC sampling and fitting settings.

## Arguments

- sample:

  numeric. Number of MCMC samples to save. Defaults to `5000`.

- burnin:

  numeric. Number of burn-in iterations. Defaults to `2000`. With at
  least two chains, the chains are checked for agreement over the last
  `min(burnin_extend, floor(burnin / thin))` samples of burn-in, taken
  every `thin` iterations (at least 100 samples), before sampling
  starts, and with `autofit` burn-in is extended until they agree (see
  `burnin_max_Rhat` in
  [`set_autofit_control()`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA_control.md)
  for when the check runs). A larger `thin` helps when more widely
  spaced samples are less correlated; as the check needs 100 samples at
  the fit's thinning, it does not run with `burnin` below `100 * thin`,
  and when it runs with `burnin` at most `burnin_extend * thin`, the
  block holds every complete window of `thin` iterations of burn-in (the
  first `burnin %% thin` iterations run unmonitored), so raise `burnin`
  together with `thin`. The fitted JAGS object (`object$fit$burnin`)
  reports the burn-in that was run; the fitting settings
  (`object$fit_control$burnin`) keep the requested value.

- adapt:

  numeric. Number of adaptation iterations. Defaults to `500`.

- chains:

  numeric. Number of MCMC chains. Defaults to `3`.

- thin:

  numeric. Thinning interval. Defaults to `1`.

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
