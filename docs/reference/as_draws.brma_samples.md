# Convert brma_samples to posterior Draws Formats

Provides an interface to the posterior package for `brma_samples`
objects. These functions convert the posterior samples to various draws
formats supported by the posterior package.

## Usage

``` r
# S3 method for class 'brma_samples'
as_draws(x, ...)

# S3 method for class 'brma_samples'
as_draws_array(x, ...)

# S3 method for class 'brma_samples'
as_draws_df(x, ...)

# S3 method for class 'brma_samples'
as_draws_list(x, ...)

# S3 method for class 'brma_samples'
as_draws_matrix(x, ...)

# S3 method for class 'brma_samples'
as_draws_rvars(x, ...)
```

## Arguments

- x:

  a `brma_samples` object

- ...:

  additional arguments passed to the corresponding posterior function

## Value

An object of the corresponding posterior draws class.

## Details

The conversion reconstructs the MCMC chain structure from the stored
`nchains` and `niter` attributes. The samples are assumed to be ordered
with chains concatenated (i.e., all iterations from chain 1, then all
from chain 2, etc.). Conditional RoBMA samples are intentionally stored
as one flattened chain because conditioning subsets posterior rows
across chains.

## See also

[`draws`](https://mc-stan.org/posterior/reference/draws.html),
[`as_draws.brma`](https://fbartos.github.io/RoBMA/reference/as_draws.brma.md)
