# PET Prior

Create PET publication-bias regression priors.

## Usage

``` r
prior_PET(
  distribution,
  parameters,
  truncation = list(lower = 0, upper = Inf),
  prior_weights = 1
)
```

## Arguments

- distribution:

  character. Prior distribution name.

- parameters:

  list. Distribution parameters.

- truncation:

  list with `lower` and `upper` truncation bounds.

- prior_weights:

  numeric prior model weight.

## Value

An object inheriting from `prior` with the `prior.PET` marker class.

## Details

This forwards to
[`BayesTools::prior_PET()`](https://fbartos.github.io/BayesTools/reference/prior_PP.html)
and uses the same distribution and parameter conventions as
[`prior()`](https://fbartos.github.io/RoBMA/reference/prior.md). By
default, PET priors are truncated to the positive half-line.
