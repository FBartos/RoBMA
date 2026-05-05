# Convert brma_samples to Matrix

Converts a `brma_samples` object to a plain matrix, removing all
brma_samples-specific attributes.

## Usage

``` r
# S3 method for class 'brma_samples'
as.matrix(x, ...)
```

## Arguments

- x:

  a `brma_samples` object

- ...:

  additional arguments (ignored)

## Value

A plain matrix of posterior samples.
