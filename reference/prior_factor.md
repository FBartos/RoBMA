# Factor Prior

Create priors on factor contrasts.

## Usage

``` r
prior_factor(
  distribution,
  parameters,
  truncation = list(lower = -Inf, upper = Inf),
  prior_weights = 1,
  contrast = "meandif"
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

- contrast:

  character. Contrast coding used for factor levels. Common RoBMA model
  options are `"treatment"`, `"meandif"`, and `"orthonormal"`; the
  standalone helper also accepts BayesTools contrast aliases such as
  `"independent"`.

## Value

An object inheriting from `prior`, `prior.factor`, and a
contrast-specific marker class.

## Details

Mean-difference and orthonormal contrasts require vector or multivariate
priors. Treatment/dummy and independent contrasts use univariate simple
priors per contrast coefficient.
