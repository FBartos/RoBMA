# Check fitted RoBMA object for errors and warnings

Checks fitted RoBMA object for warnings and errors and prints them to
the console.

## Usage

``` r
check_RoBMA(fit)

check_RoBMA_convergence(fit)
```

## Arguments

- fit:

  a fitted RoBMA object.

## Value

`check_RoBMA` returns a vector of error and warning messages.
`check_RoBMA_convergence` returns a logical vector indicating whether
the models have converged.
