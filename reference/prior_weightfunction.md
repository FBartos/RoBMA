# Weightfunction Prior

Create weightfunction publication-bias priors and their weight-prior
helper objects.

## Usage

``` r
prior_weightfunction(
  side = "one-sided",
  steps = c(0.025, 0.05),
  weights = wf_cumulative(),
  reference = "most_significant",
  prior_weights = 1
)

wf_cumulative(alpha = NULL)

wf_fixed(omega)

wf_independent(prior, scale = "omega")
```

## Arguments

- side:

  character. Either `"one-sided"` or `"two-sided"`.

- steps:

  numeric vector of p-value cut points. These define `length(steps) + 1`
  p-value bins and must be ordered values in `(0, 1)`.

- weights:

  a weight-prior object created by `wf_cumulative()`, `wf_fixed()`, or
  `wf_independent()`.

- reference:

  character. Reference bin, currently `"most_significant"`.

- prior_weights:

  numeric prior model weight.

- alpha:

  optional positive cumulative-Dirichlet concentration parameters, one
  per p-value bin. If `NULL`, `prior_weightfunction()` uses
  `rep(1, length(steps) + 1)`. Cumulative weights encode monotone
  decreasing publication weights relative to the most-significant bin.

- omega:

  fixed publication weights, one per bin; values must be non-missing,
  nonnegative, and match `length(steps) + 1` when used in
  `prior_weightfunction()`.

- prior:

  continuous simple prior distribution for each non-reference weight.
  Point, discrete, mixture, and other non-simple priors are invalid.

- scale:

  latent scale for independent weights; either `"omega"`, `"log_omega"`,
  or the `"log"` alias. Direct `"omega"` priors need nonnegative
  support; `"log"` is normalized to `"log_omega"`.

## Value

`prior_weightfunction()` returns an object inheriting from `prior` and
`prior.weightfunction`; the `wf_*()` helpers return
`weightfunction_weights` helper objects with subclass markers.

## Details

Fixed weights must have one value per p-value bin (`length(steps) + 1`),
and the reference bin must have weight 1.

## See also

[`publication_bias_prior_specification`](https://fbartos.github.io/RoBMA/reference/publication_bias_prior_specification.md)

## Examples

``` r
prior_weightfunction("one-sided", steps = 0.025)
#> omega[one-sided: .025] ~ CumDirichlet(1, 1)
prior_weightfunction(
  side    = "one-sided",
  steps   = c(0.025, 0.5),
  weights = wf_fixed(c(1, 0.8, 0.6))
)
#> omega[one-sided: .025, .5] = (1, 0.8, 0.6)
```
