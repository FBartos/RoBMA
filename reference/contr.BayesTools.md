# BayesTools Contrast Matrices

BayesTools provides several contrast matrix functions for Bayesian
factor analysis. These functions create different types of contrast
matrices that can be used with factor variables in Bayesian models.

## Usage

``` r
contr.orthonormal(n, contrasts = TRUE)

contr.meandif(n, contrasts = TRUE)

contr.independent(n, contrasts = TRUE)
```

## Arguments

- n:

  a vector of levels for a factor, or the number of levels

- contrasts:

  logical indicating whether contrasts should be computed

## Details

The package includes the following contrast functions:

- `contr.orthonormal`:

  Return a matrix of orthonormal contrasts. Code is based on
  `stanova::contr.bayes` and corresponding to description by Rouder et
  al. (2012) . Returns a matrix with n rows and k columns, with k = n -
  1 if `contrasts = TRUE` and k = n if `contrasts = FALSE`.

- `contr.meandif`:

  Return a matrix of mean difference contrasts. This is an adjustment to
  the `contr.orthonormal` that ascertains that the prior distributions
  on difference between the grand mean and factor level are identical
  independent of the number of factor levels (which does not hold for
  the orthonormal contrast). Furthermore, the contrast is re-scaled so
  the specified prior distribution exactly corresponds to the prior
  distribution on difference between each factor level and the grand
  mean – this is approximately twice the scale of `contr.orthonormal`.
  Returns a matrix with n rows and k columns, with k = n - 1 if
  `contrasts = TRUE` and k = n if `contrasts = FALSE`.

- `contr.independent`:

  Return a matrix of independent contrasts – a level for each term.
  Returns a matrix with n rows and k columns, with k = n if
  `contrasts = TRUE` and k = n if `contrasts = FALSE`.

## References

Rouder JN, Morey RD, Speckman PL, Province JM (2012). “Default Bayes
factors for ANOVA designs.” *Journal of Mathematical Psychology*,
**56**(5), 356–374.
[doi:10.1016/j.jmp.2012.08.001](https://doi.org/10.1016/j.jmp.2012.08.001)
.

## Examples

``` r
# Orthonormal contrasts
contr.orthonormal(c(1, 2))
#>            [,1]
#> [1,] -0.7071068
#> [2,]  0.7071068
contr.orthonormal(c(1, 2, 3))
#>            [,1]       [,2]
#> [1,]  0.0000000  0.8164966
#> [2,] -0.7071068 -0.4082483
#> [3,]  0.7071068 -0.4082483

# Mean difference contrasts
contr.meandif(c(1, 2))
#>      [,1]
#> [1,]   -1
#> [2,]    1
contr.meandif(c(1, 2, 3))
#>            [,1] [,2]
#> [1,]  0.0000000  1.0
#> [2,] -0.8660254 -0.5
#> [3,]  0.8660254 -0.5

# Independent contrasts
contr.independent(c(1, 2))
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
contr.independent(c(1, 2, 3))
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
```
