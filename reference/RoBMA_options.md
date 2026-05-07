# Options for the RoBMA package

Inspect or change package-level defaults used by RoBMA.

## Usage

``` r
RoBMA.options(...)

RoBMA.get_option(name)
```

## Arguments

- ...:

  named option(s) to change. Names must be exact public option names,
  nonempty, and unique.

- name:

  a single non-missing character string matching one public option
  exactly; for available options, see details below.

## Value

`RoBMA.options()` invisibly returns a named list with all current
options after applying any changes. `RoBMA.get_option()` returns the
current value of the requested option.

## Details

The available options are:

- `max_cores`:

  number of cores to use for parallel computing (default is one fewer
  than detected logical cores, with a minimum/fallback of 1)

- `check_scaling`:

  whether to check scaling of predictors (default `TRUE`)

- `silent`:

  whether to suppress output (default `FALSE`)

- `autocompute.loo`:

  whether to automatically compute LOO (default `FALSE`)

- `autocompute.waic`:

  whether to automatically compute WAIC (default `FALSE`)

- `autocompute.marglik`:

  whether to automatically compute marginal likelihood (default `FALSE`)

- `cluster_likelihood.n_gamma`:

  number of Gauss-Hermite nodes used for cluster-unit log-likelihoods
  (default `15`)

- `default_UISD.effect`:

  default scaling of the unit information standard deviation for the
  effect size parameter (default `0.5`)

- `default_UISD.heterogeneity`:

  default scaling of the unit information standard deviation for the
  heterogeneity parameter (default `0.25`)

- `default_UISD.mods`:

  default scaling of the unit information standard deviation for the
  moderators (default `0.25`)

- `default_UISD.scale`:

  default scaling of the unit information standard deviation for the
  scale parameter (default `0.5`)

- `default_informed_priors.mods`:

  default scaling of informed priors for moderators (default `0.5`)

- `default_informed_priors.scale`:

  default scaling of informed priors for the scale parameter (default
  `0.5`)

- `default_bias_weightfunction.alpha`:

  default alpha for the weightfunction (default `1`)

- `default_bias_PET.scale`:

  default scale for the PET (default `1`)

- `default_bias_PEESE.scale`:

  default scale for the PEESE (default `5`)

Boolean options require scalar `TRUE` or `FALSE` values. `max_cores`
must be integer-like and at least 1; `cluster_likelihood.n_gamma` must
be integer-like and at least 3. Scale and alpha defaults must be finite
positive numbers.
