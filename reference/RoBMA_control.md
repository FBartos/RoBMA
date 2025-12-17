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
  sample_extend = 1000,
  restarts = 10,
  max_extend = 10
)

set_convergence_checks(
  max_Rhat = 1.05,
  min_ESS = 500,
  max_error = NULL,
  max_SD_error = NULL,
  remove_failed = FALSE,
  balance_probability = TRUE
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
  function (possible units are
  `"secs", "mins", "hours", "days", "weeks", "years"`). Defaults to
  `list(time = 60, unit = "mins")`.

- sample_extend:

  number of samples to extend the fitting process if the criteria are
  not satisfied. Defaults to `1000`.

- restarts:

  number of times new initial values should be generated in case a model
  fails to initialize. Defaults to `10`.

- max_extend:

  number of times after which the automatic fitting function is stopped.
  Defaults to `10`.

- remove_failed:

  whether models not satisfying the convergence checks should be removed
  from the inference. Defaults to `FALSE` - only a warning is raised.

- balance_probability:

  whether prior model probability should be balanced across the
  combinations of models with the same H0/H1 for effect / heterogeneity
  / bias in the case of non-convergence. Defaults to `TRUE`.

## Value

`set_autofit_control` returns a list of autofit control settings and
`set_convergence_checks` returns a list of convergence checks settings.

## Details

`set_autofit_control` controls the automatic model fitting process,
determining when to stop iterating based on convergence criteria.
`set_convergence_checks` sets thresholds for determining whether a model
has converged adequately.

The autofit control manages computational resources by setting maximum
time limits and determining how many additional samples to draw if
convergence criteria are not met. The convergence checks determine the
quality standards that fitted models must meet.

## See also

[`RoBMA()`](reference/RoBMA.md),
[`update.RoBMA()`](reference/update.RoBMA.md)

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
# Use in RoBMA function
fit <- RoBMA(d = c(0.5, 0.3, 0.1), 
             se = c(0.2, 0.15, 0.1),
             autofit_control = autofit_ctrl,
             convergence_checks = conv_checks)
} # }
```
