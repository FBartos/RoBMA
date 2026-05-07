# Update a brma Fit

Extends an existing fitted `brma` object by additional MCMC samples,
updates study labels, and optionally recomputes cached fit-dependent
quantities.

## Usage

``` r
# S3 method for class 'brma'
update(
  object,
  formula. = NULL,
  ...,
  sample_extend = NULL,
  slab = NULL,
  autofit_control = NULL,
  convergence_checks = NULL,
  recompute = c("all", "drop"),
  parallel = NULL,
  cores = NULL,
  silent = NULL,
  seed = NULL,
  evaluate = TRUE
)
```

## Arguments

- object:

  a fitted `brma` object.

- formula.:

  unsupported; included for compatibility with
  [`update`](https://rdrr.io/r/stats/update.html).

- ...:

  unsupported additional arguments.

- sample_extend:

  integer. Number of additional samples per chain.

- slab:

  optional character vector of study labels. Updating labels does not
  refit or extend the model.

- autofit_control:

  list of autofit control settings. Values are merged with the existing
  settings before extending.

- convergence_checks:

  list of convergence check settings. Values are merged with the
  existing settings and used to re-check the fit.

- recompute:

  whether cached `loo`, `waic`, and `marglik` values already stored in
  `object` should be recomputed after extension (`"all"`) or dropped
  with a warning (`"drop"`).

- parallel:

  logical. Whether to extend chains in parallel.

- cores:

  integer. Number of cores to use when `parallel = TRUE`.

- silent:

  logical. Whether to suppress JAGS output during extension.

- seed:

  optional seed used before extending.

- evaluate:

  unsupported; included for compatibility with
  [`update`](https://rdrr.io/r/stats/update.html).

## Value

The updated `brma` object.

## Details

Extending a fit adds posterior samples only. It does not rerun
adaptation or burn-in. Prior, data, and model-structure updates are
intentionally not supported by this method.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- update(fit, sample_extend = 1000)
fit <- update(fit, slab = paste("Study", seq_len(nobs(fit))))
} # }
```
