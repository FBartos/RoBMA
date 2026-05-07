# Extract Log-Likelihood Matrix from brma Object

Extract the pointwise log-likelihood matrix from a brma model object.
This is an S x K or S x G matrix where S is the number of posterior
samples, K is the number of estimates, and G is the number of clusters.
This method implements the S3 \\logLik\\ generic for `brma` objects and
returns the matrix of pointwise log-likelihoods (one column per
observation, one row per sample).

## Usage

``` r
# S3 method for class 'brma'
logLik(object, unit = "estimate", ...)
```

## Arguments

- object:

  a brma model object.

- unit:

  output unit. See
  [`add_loo`](https://fbartos.github.io/RoBMA/reference/add_loo.brma.md).

- ...:

  currently unused.

## Value

An S x K or S x G matrix of log-likelihood values.

## Details

The log-likelihood is computed for each observation at each posterior
sample. For binomial and Poisson models, each observation consists of a
pair of counts (ai/ci or x1i/x2i) that together define a single effect
size estimate.

## See also

[`loo.brma`](https://fbartos.github.io/RoBMA/reference/loo.brma.md)
