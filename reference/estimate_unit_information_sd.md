# Estimate Unit Information Standard Deviation

Estimates the unit information standard deviation (UISD) from sample
sizes and standard errors. The UISD is used to scale weakly informative
prior distributions for meta-analytic parameters.

## Usage

``` r
estimate_unit_information_sd(sei, ni)
```

## Arguments

- sei:

  a complete numeric vector of strictly positive standard errors for
  each study.

- ni:

  a complete numeric vector of strictly positive sample sizes for each
  study, with the same length as `sei`.

## Value

Returns a single numeric value representing the estimated unit
information standard deviation.

## Details

The unit information standard deviation is computed following Equation 6
in Röver et al. (2021) : \$\$\text{UISD} = \sqrt{\frac{\sum n_i}{\sum
\text{se}\_i^{-2}}}\$\$ where \\n_i\\ is the sample size and
\\\text{se}\_i\\ is the standard error for each study.

This function is useful when you want to compute the UISD once and pass
it to multiple analyses via the `prior_unit_information_sd` argument in
[`brma()`](https://fbartos.github.io/RoBMA/reference/brma.md) or related
functions. This ensures consistent prior scaling across analyses
performed on different subsets of the same data.

## References

Röver C, Bender R, Dias S, Schmid CH, Schmidli H, Sturtz S, Weber S,
Friede T (2021). “On weakly informative prior distributions for the
heterogeneity parameter in Bayesian random-effects meta-analysis.”
*Research Synthesis Methods*, **12**(4), 448–474.
[doi:10.1002/jrsm.1475](https://doi.org/10.1002/jrsm.1475) .

## Examples

``` r
# Example with simulated data
sei <- c(0.2, 0.3, 0.25, 0.15)
ni  <- c(50, 30, 40, 80)
estimate_unit_information_sd(sei = sei, ni = ni)
#> [1] 1.439217
```
