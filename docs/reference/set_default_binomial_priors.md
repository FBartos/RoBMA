# Set default prior distributions for binomial meta-analytic models

Set default prior distributions for BiBMA models.

## Usage

``` r
set_default_binomial_priors(parameter, null = FALSE, rescale = 1)
```

## Arguments

- parameter:

  a character string specifying the parameter for which the prior
  distribution should be set. Available options are "effect",
  "heterogeneity", "baseline", "covariates", "factors".

- null:

  a logical indicating whether the prior distribution should be set for
  the null hypothesis. Defaults to `FALSE`.

- rescale:

  a numeric value specifying the re-scaling factor for the default prior
  distributions. Defaults to 1. Allows convenient re-scaling of prior
  distributions simultaneously.

## Value

A prior distribution object or a list of prior distribution objects.

## Details

The default prior are based on the binary outcome meta-analyses in the
Cochrane Database of Systematic Reviews outlined in Bartoš et al. (2023)
.

Specifically, the prior distributions are:

**For the alternative hypothesis:**

- **Effect:** T distribution with mean 0, scale 0.58, and 4 degrees of
  freedom.

- **Heterogeneity:** Inverse gamma distribution with shape 1.77 and
  scale 0.55.

- **Baseline:** No prior distribution.

- **Standardized continuous covariates:** Normal distribution with mean
  0 and standard deviation 0.29.

- **Factors (via by-level differences from the grand mean):** Normal
  distribution with mean 0 and standard deviation 0.29.

**For the null hypothesis:**

- **Effect:** Point distribution at 0.

- **Heterogeneity:** Point distribution at 0.

- **Baseline:** Independent uniform distributions.

- **Standardized continuous covariates:** Point distribution at 0.

- **Factors (via by-level differences from the grand mean):** Point
  distribution at 0.

The rescaling factor adjusts the width of the effect, heterogeneity,
covariates, factor, and PEESE-style model prior distributions. PET-style
and weight function prior distributions are scale-invariant.

## Examples

``` r

set_default_binomial_priors("effect")
#> Student-t(0, 0.58, 4)
set_default_binomial_priors("heterogeneity")
#> InvGamma(1.77, 0.55)
set_default_binomial_priors("baseline")
#> NULL
```
