# Set default prior distributions

Set default prior distributions for RoBMA models.

## Usage

``` r
set_default_priors(parameter, null = FALSE, rescale = 1)
```

## Arguments

- parameter:

  a character string specifying the parameter for which the prior
  distribution should be set. Available options are "effect",
  "heterogeneity", "bias", "hierarchical", "covariates", "factors".

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

The default prior distributions corresponds to the specification of
RoBMA-PSMA and RoBMA-regression outlined in Bartoš et al. (2023) and
Bartoš et al. (2025) .

Specifically, the prior distributions are:

**For the alternative hypothesis:**

- **Effect:** Normal distribution with mean 0 and standard deviation 1.

- **Heterogeneity:** Inverse gamma distribution with shape 1 and scale
  0.15.

- **Bias:** A list of 8 prior distributions defining the publication
  bias adjustments:

  - Two-sided: Weight function with steps 0.05.

  - Two-sided: Weight function with steps 0.05 and 0.1.

  - One-sided: Weight function with steps 0.05.

  - One-sided: Weight function with steps 0.025 and 0.05.

  - One-sided: Weight function with steps 0.05 and 0.5.

  - One-sided: Weight function with steps 0.025, 0.05, and 0.5.

  - PET-type model with regression coefficient: Cauchy distribution with
    location 0 and scale 1.

  - PEESE-type model with regression coefficient: Cauchy distribution
    with location 0 and scale 5.

  All weight functions use a unit cumulative Dirichlet prior
  distribution on relative prior probabilities.

- **Standardized continuous covariates:** Normal distribution with mean
  0 and standard deviation 0.25.

- **Factors (via by-level differences from the grand mean):** Normal
  distribution with mean 0 and standard deviation 0.25.

**For the null hypothesis:**

- **Effect:** Point distribution at 0.

- **Heterogeneity:** Point distribution at 0.

- **Bias:** No prior distribution.

- **Standardized continuous covariates:** Point distribution at 0.

- **Factors (via by-level differences from the grand mean):** Point
  distribution at 0.

The rescaling factor adjusts the width of the effect, heterogeneity,
covariates, factor, and PEESE-style model prior distributions. PET-style
and weight function prior distributions are scale-invariant.

## Examples

``` r

set_default_priors("effect")
#> Normal(0, 1)
set_default_priors("heterogeneity")
#> InvGamma(1, 0.15)
set_default_priors("bias")
#> [[1]]
#> omega[two-sided: .05] ~ CumDirichlet(1, 1)
#> [[2]]
#> omega[two-sided: .1, .05] ~ CumDirichlet(1, 1, 1)
#> [[3]]
#> omega[one-sided: .05] ~ CumDirichlet(1, 1)
#> [[4]]
#> omega[one-sided: .05, .025] ~ CumDirichlet(1, 1, 1)
#> [[5]]
#> omega[one-sided: .5, .05] ~ CumDirichlet(1, 1, 1)
#> [[6]]
#> omega[one-sided: .5, .05, .025] ~ CumDirichlet(1, 1, 1, 1)
#> [[7]]
#> PET ~ Cauchy(0, 1)[0, Inf]
#> [[8]]
#> PEESE ~ Cauchy(0, 5)[0, Inf]
```
