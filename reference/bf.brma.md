# Bayes Factor for brma Objects

Compute the Bayes factor comparing two brma models.

## Usage

``` r
# S3 method for class 'brma'
bf(x1, x2, log = FALSE, ...)

# S3 method for class 'brma'
bayes_factor(x1, x2, log = FALSE, ...)
```

## Arguments

- x1:

  a brma model object (numerator).

- x2:

  a brma model object (denominator).

- log:

  logical; if `TRUE`, the log Bayes factor is returned. Default is
  `FALSE`.

- ...:

  additional arguments (currently not used).

## Value

A list of class `"bf_default"` with components:

- `bf`: (scalar) value of the Bayes factor in favor of `x1` over `x2`.

- `log`: Boolean indicating whether `bf` corresponds to the log Bayes
  factor.

## Details

Computes the Bayes factor in favor of the model `x1` over the model
`x2`. The marginal likelihoods must first be computed using
[`add_marglik`](https://fbartos.github.io/RoBMA/reference/add_marglik.brma.md).
Both models must be fitted to the same outcome target/data, including
outcome type and, when present, weights and cluster identifiers.

## See also

[`add_marglik`](https://fbartos.github.io/RoBMA/reference/add_marglik.brma.md),
[`bridge_sampler.brma`](https://fbartos.github.io/RoBMA/reference/bridge_sampler.brma.md),
[`post_prob.brma`](https://fbartos.github.io/RoBMA/reference/post_prob.brma.md),
[`logml.brma`](https://fbartos.github.io/RoBMA/reference/logml.brma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")
  fit1 <- brma(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
  fit2 <- brma(
    yi           = yi,
    vi           = vi,
    data         = dat.lehmann2018,
    measure      = "SMD",
    prior_effect = FALSE
  )

  fit1 <- add_marglik(fit1)
  fit2 <- add_marglik(fit2)

  bf(fit1, fit2)
}
} # }
```
