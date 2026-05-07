# Convert brma Objects to posterior Draws Formats

Provides an interface to the posterior package for `brma` objects. These
functions convert the MCMC samples from a fitted `brma` model to various
draws formats supported by the posterior package, enabling the use of
posterior's rich set of diagnostics and summary functions.

## Usage

``` r
as_draws(x, ...)

as_draws_array(x, ...)

as_draws_df(x, ...)

as_draws_list(x, ...)

as_draws_matrix(x, ...)

as_draws_rvars(x, ...)

# Default S3 method
as_draws(x, ...)

# Default S3 method
as_draws_array(x, ...)

# Default S3 method
as_draws_df(x, ...)

# Default S3 method
as_draws_list(x, ...)

# Default S3 method
as_draws_matrix(x, ...)

# Default S3 method
as_draws_rvars(x, ...)

# S3 method for class 'brma'
as_draws(x, ...)

# S3 method for class 'brma'
as_draws_array(x, ...)

# S3 method for class 'brma'
as_draws_df(x, ...)

# S3 method for class 'brma'
as_draws_list(x, ...)

# S3 method for class 'brma'
as_draws_matrix(x, ...)

# S3 method for class 'brma'
as_draws_rvars(x, ...)
```

## Arguments

- x:

  an object to convert. The `brma` methods expect a fitted `brma`
  object; default methods forward non-`brma` objects to the
  corresponding posterior conversion function.

- ...:

  additional arguments passed to the corresponding posterior function.

## Value

An object of the corresponding posterior draws class.

## Details

These functions are S3 generics. Their default methods forward to the
corresponding posterior generics so attaching RoBMA preserves the usual
posterior behavior for non-`brma` objects.

The following conversion functions are available:

- `as_draws`: converts to the default draws format

- `as_draws_array`: converts to a 3-D array (iteration x chain x
  variable)

- `as_draws_df`: converts to a data frame with columns for iterations,
  chains, and variables

- `as_draws_list`: converts to a list of lists

- `as_draws_matrix`: converts to a 2-D matrix (draw x variable)

- `as_draws_rvars`: converts to random variable objects

These methods require the posterior package to be installed. For `brma`
methods, conversion is performed by first extracting the MCMC samples as
a `mcmc.list` object and then using the corresponding posterior
conversion function. `brma_samples` objects have separate methods
documented at
[`as_draws.brma_samples`](https://fbartos.github.io/RoBMA/reference/as_draws.brma_samples.md).

## See also

[`draws`](https://mc-stan.org/posterior/reference/draws.html),
[`brma`](https://fbartos.github.io/RoBMA/reference/brma.md),
[`as_draws.brma_samples`](https://fbartos.github.io/RoBMA/reference/as_draws.brma_samples.md)
