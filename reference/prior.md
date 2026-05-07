# Prior Distribution

Create prior distribution objects used by RoBMA and brma fitting
functions.

## Usage

``` r
prior(
  distribution,
  parameters,
  truncation = list(lower = -Inf, upper = Inf),
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

An object inheriting from `prior`.

## Details

This is RoBMA's re-export of
[`BayesTools::prior()`](https://fbartos.github.io/BayesTools/reference/prior.html);
supported distribution names and parameter lists follow BayesTools.

## Examples

``` r
prior("normal", list(mean = 0, sd = 0.5))
#> Normal(0, 0.5)
prior(
  "normal",
  list(mean = 0, sd = 0.5),
  truncation = list(lower = 0, upper = Inf)
)
#> Normal(0, 0.5)[0, Inf]
```
