# Reports whether x is a RoBMA object

Functions to test whether an object is of a specific RoBMA class.

## Usage

``` r
is.RoBMA(x)

is.RoBMA.reg(x)

is.NoBMA(x)

is.NoBMA.reg(x)

is.BiBMA(x)

is.BiBMA.reg(x)
```

## Arguments

- x:

  object to be tested

## Value

returns a boolean.

`TRUE` if the object inherits from the specified class, `FALSE`
otherwise.

## Details

These functions test whether an object inherits from specific RoBMA
classes:

- `is.RoBMA`: Tests for `"RoBMA"` class (Robust Bayesian Meta-Analysis)

- `is.RoBMA.reg`: Tests for `"RoBMA.reg"` class (RoBMA with
  meta-regression)

- `is.NoBMA`: Tests for `"NoBMA"` class (Normal-normal Bayesian
  Meta-Analysis)

- `is.NoBMA.reg`: Tests for `"NoBMA.reg"` class (NoBMA with
  meta-regression)

- `is.BiBMA`: Tests for `"BiBMA"` class (Binomial-normal Bayesian
  Meta-Analysis)

- `is.BiBMA.reg`: Tests for `"BiBMA.reg"` class (BiBMA with
  meta-regression)

## See also

[`RoBMA()`](reference/RoBMA.md),
[`RoBMA.reg()`](reference/RoBMA.reg.md),
[`NoBMA()`](reference/NoBMA.md), [`BiBMA()`](reference/BiBMA.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Example with Anderson et al. 2010 data
fit <- RoBMA(r = Anderson2010$r, n = Anderson2010$n)
is.RoBMA(fit)        # TRUE
is.BiBMA(fit)        # FALSE

# Example with regression
fit_reg <- RoBMA.reg(r ~ 1, data = Anderson2010)
is.RoBMA.reg(fit_reg)  # TRUE
is.RoBMA(fit_reg)      # TRUE (inherits from RoBMA)
} # }
```
