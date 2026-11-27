# Control MCMC fitting process

`set_autofit_control` and `set_convergence_checks` control settings for
the autofit process of the MCMC JAGS sampler and specify termination
criteria and convergence checks.

## Usage

``` r
set_autofit_control(
  max_Rhat = 1.05,
  min_ESS = 500,
  max_error = NULL,
  max_SD_error = NULL,
  max_time = list(time = 60, unit = "mins"),
  sample_extend = "auto",
  restarts = 10,
  max_extend = 10,
  check_indicators = FALSE,
  monitor = NULL,
  burnin_max_Rhat = "auto",
  burnin_extend = 500,
  burnin_max_extend = 10
)

set_convergence_checks(
  max_Rhat = 1.05,
  min_ESS = 500,
  max_error = NULL,
  max_SD_error = NULL,
  check_indicators = FALSE,
  monitor = NULL,
  allow_not_assessable = FALSE
)
```

## Arguments

- max_Rhat:

  maximum value of the R-hat diagnostic. Defaults to `1.05`.

- min_ESS:

  minimum estimated sample size. Defaults to `500`.

- max_error:

  maximum value of the MCMC error. Defaults to `NULL`. Be aware that
  PEESE publication bias adjustment can have estimates on different
  scale than the rest of the output, resulting in relatively large max
  MCMC error.

- max_SD_error:

  maximum value of the proportion of MCMC error of the estimated SD of
  the parameter. Defaults to `NULL`.

- max_time:

  list with the time and unit specifying the maximum autofitting process
  per model. Passed to [difftime](https://rdrr.io/r/base/difftime.html)
  function (possible units are `"secs"`, `"mins"`, `"hours"`, `"days"`,
  and `"weeks"`). Defaults to `list(time = 60, unit = "mins")`.

- sample_extend:

  number of samples per chain by which each automatic extension extends
  the fit when the criteria are not satisfied, or `"auto"` (the default)
  to size each extension from the convergence check that failed:
  BayesTools projects the additional draws that the parameter furthest
  from its criteria needs (from its ESS and MCMC errors, and at least as
  many draws again for an R-hat above `max_Rhat`; diagnostics that are
  not assessable add nothing) and adds a tenth, rounded up to a multiple
  of 100; the extension is then at most twice the current draws per
  chain (rounded down to a multiple of 100), at least 1000, and at most
  the draws per chain that the time left under `max_time` allows at the
  rate observed so far, including adaptation and burn-in (rounded down
  to a multiple of 100; automatic fitting stops when that is below 100).
  Without `max_time`, each extension of a model that cannot converge
  multiplies the draws per chain by about 2.1 (a failed R-hat) and up to
  3 (large ESS or error deficits), up to `max_extend`. See
  `autofit_control` in
  [`JAGS_fit`](https://fbartos.github.io/BayesTools/reference/JAGS_fit.html).

- restarts:

  number of times new initial values should be generated in case a model
  fails to initialize. Defaults to `10`.

- max_extend:

  number of times after which the automatic fitting function is stopped.
  Defaults to `10`.

- check_indicators:

  whether model indicator variables should be included in convergence
  checks. Binary indicators are checked as state occupancies and
  categorical indicators are checked separately for every observed
  state. Eligible indicators remain included when `monitor` selects a
  narrower parameter set; auxiliary inclusion-probability coordinates
  remain excluded unless requested explicitly. Defaults to `FALSE`.

- monitor:

  optional character vector selecting the parameters used for
  convergence checks. A base name selects all indexed elements. Defaults
  to `NULL`, which checks every eligible parameter.

- burnin_max_Rhat:

  threshold of the readiness check between burn-in and sampling.
  Defaults to `"auto"`, which is `1.7 + 0.9 / chains` (`2.15` for two
  chains, `2` for three, `1.925` for four); a number fixes the threshold
  and `NULL` turns the check off. The check runs, also with
  `autofit = FALSE`, in the initial fit of a model (not when
  [`update()`](https://rdrr.io/r/stats/update.html) extends the chains,
  nor in the automatic extensions after sampling) when there are at
  least two chains, the block of the last `min(burnin_extend, burnin)`
  samples of burn-in, taken every `thin` iterations as in sampling, has
  at least 100 samples, and JAGS assigns the model a sampler. The
  iterations of the block run with the monitors on and their samples are
  checked and then dropped: the chains are ready when the split-chain
  R-hat point estimate (each chain's block split into halves) is at most
  `burnin_max_Rhat` for every parameter of the convergence check
  (`monitor`, `check_indicators`). With `autofit`, burn-in of chains
  that are not ready is extended by blocks of `burnin_extend` samples
  (`burnin_extend * thin` iterations) until they are,
  `burnin_max_extend` blocks were added, or `max_time` is exceeded;
  without it, burn-in is kept. Chains that are still not ready start
  sampling, and a warning is issued once that sampling has finished (an
  attempt that is restarted warns nothing). In a model to which JAGS
  assigns no sampler, whose chains cannot disagree, burn-in runs as with
  the check off. A fit whose chains are ready at once has the draws it
  would have with the check off, unless a monitored node has no observed
  descendants (see `autofit_control` in
  [`JAGS_fit`](https://fbartos.github.io/BayesTools/reference/JAGS_fit.html)).
  The automatic value falls with the number of chains: for chains that
  agree but are autocorrelated, the split-chain R-hat has the same
  typical value for any number of chains but varies less the more chains
  are run, so a threshold that falls with the number of chains keeps
  false alarms from becoming more frequent with them and keeps part of
  the sensitivity to one deviating chain that a fixed value loses. A
  slowly mixing model can exceed the threshold although its chains
  agree, because the statistic grows with the autocorrelation of the
  block's samples; a larger `thin` then helps when more widely spaced
  samples are less correlated, while more burn-in is the remedy for
  chains that have not yet reached the same region. With `burnin` below
  100 the check does not run, and with `burnin` at most `burnin_extend`
  it checks the whole burn-in (`burnin` counts samples at the fit's
  thinning, see
  [`fitting_specification`](https://fbartos.github.io/RoBMA/dev/reference/fitting_specification.md)).
  When the warning appears for chains that started in different
  components (for example of the publication-bias mixture), a longer
  adaptation (such as `adapt = 5000`) resolves it sooner than more
  burn-in alone.

- burnin_extend:

  number of final burn-in samples, taken every `thin` iterations, that
  the readiness check assesses, and of each extension of burn-in
  (`burnin_extend * thin` iterations); at least `100`. Defaults to
  `500`.

- burnin_max_extend:

  largest number of extensions of burn-in, with `autofit` only. Defaults
  to `10`.

- allow_not_assessable:

  whether requested sampled parameters with undefined diagnostics may be
  ignored by the convergence checks. Defaults to `FALSE`, with which the
  fitted model reports the parameters whose diagnostics are not
  assessable. Automatic fitting does not extend a model for diagnostics
  that cannot be assessed from the samples available at a check (a
  mixture component or indicator state that no chain has visited,
  samples that are constant), so a parameter that only more samples
  would reveal is not waited for; every diagnostic that can be assessed
  still extends the model.

## Value

`set_autofit_control` returns a list of autofit control settings
including `max_Rhat`, `min_ESS`, `max_error`, `max_SD_error`,
`max_time`, `sample_extend`, `restarts`, `max_extend`,
`check_indicators`, `monitor`, `burnin_max_Rhat`, `burnin_extend`, and
`burnin_max_extend`. `set_convergence_checks` returns a list with
convergence thresholds `max_Rhat`, `min_ESS`, `max_error`,
`max_SD_error`, `check_indicators`, `monitor`, and
`allow_not_assessable`.

## Details

`set_autofit_control` controls the automatic model fitting process,
determining when to stop iterating based on convergence criteria.
`set_convergence_checks` sets thresholds for determining whether a model
has converged adequately.

The autofit control manages computational resources by setting maximum
time limits and determining how many additional samples to draw if
convergence criteria are not met. Automatic fitting stops once every
diagnostic that can be assessed from the samples available at a check
meets its criterion, also when other diagnostics are not assessable
then, or at `max_time` or `max_extend`;
`set_convergence_checks(allow_not_assessable = )` decides whether the
fitted model reports the parameters whose diagnostics are not
assessable. The convergence checks determine the quality standards that
fitted models must meet. Autofit thresholds `max_Rhat`, `min_ESS`,
`max_error`, `max_SD_error`, `max_time`, `restarts`, `max_extend`, and
`burnin_max_Rhat` can be set to `NULL`; `sample_extend` must be a
positive integer or `"auto"`, and `burnin_max_Rhat` a finite number of
at least 1 or `"auto"`. The routing controls are delegated to
[`BayesTools::JAGS_check_convergence()`](https://fbartos.github.io/BayesTools/reference/JAGS_check_convergence.html).
Thresholds can be disabled with `NULL`; otherwise `max_Rhat` must be at
least 1, `min_ESS` and `max_error` must be nonnegative, and
`max_SD_error` must be between 0 and 1.

## See also

[`RoBMA()`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.md)

## Examples

``` r
# Set custom autofit control with shorter time limit
autofit_ctrl <- set_autofit_control(
  max_time = list(time = 30, unit = "mins"),
  sample_extend = 2000
)

# Set custom convergence checks with stricter criteria
conv_checks <- set_convergence_checks(
  max_Rhat = 1.01,
  min_ESS = 1000
)

if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")

  fit <- RoBMA(
    yi                 = yi,
    vi                 = vi,
    data               = dat.lehmann2018,
    measure            = "SMD",
    autofit_control    = autofit_ctrl,
    convergence_checks = conv_checks
  )
}
} # }
```
