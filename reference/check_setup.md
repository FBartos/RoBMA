# Prints summary of `"RoBMA"` ensemble implied by the specified priors

`check_setup` prints summary of `"RoBMA"` ensemble implied by the
specified prior distributions. It is useful for checking the ensemble
configuration prior to fitting all of the models.

## Usage

``` r
check_setup(
  model_type = NULL,
  priors_effect = prior(distribution = "normal", parameters = list(mean = 0, sd = 1)),
  priors_heterogeneity = prior(distribution = "invgamma", parameters = list(shape = 1,
    scale = 0.15)),
  priors_bias = list(prior_weightfunction(distribution = "two.sided", parameters =
    list(alpha = c(1, 1), steps = c(0.05)), prior_weights = 1/12),
    prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1,
    1), steps = c(0.05, 0.1)), prior_weights = 1/12), prior_weightfunction(distribution =
    "one.sided", parameters = list(alpha = c(1, 1), steps = c(0.05)), prior_weights =
    1/12), prior_weightfunction(distribution = "one.sided", parameters = list(alpha =
    c(1, 1, 1), steps = c(0.025, 0.05)), prior_weights = 1/12), 
    
    prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1,
    1), steps = c(0.05, 0.5)), prior_weights = 1/12), prior_weightfunction(distribution =
    "one.sided", parameters = list(alpha = c(1, 1, 1, 1), steps = c(0.025, 0.05, 0.5)),
    prior_weights = 1/12), prior_PET(distribution = "Cauchy", parameters = list(0, 1),
    truncation = list(0, Inf), prior_weights = 1/4), prior_PEESE(distribution = "Cauchy",
    parameters = list(0, 5), truncation = list(0, Inf), prior_weights = 1/4)),
  priors_effect_null = prior(distribution = "point", parameters = list(location = 0)),
  priors_heterogeneity_null = prior(distribution = "point", parameters = list(location =
    0)),
  priors_bias_null = prior_none(),
  priors_hierarchical = prior("beta", parameters = list(alpha = 1, beta = 1)),
  priors_hierarchical_null = NULL,
  models = FALSE,
  silent = FALSE
)

check_setup.RoBMA(
  model_type = NULL,
  priors_effect = prior(distribution = "normal", parameters = list(mean = 0, sd = 1)),
  priors_heterogeneity = prior(distribution = "invgamma", parameters = list(shape = 1,
    scale = 0.15)),
  priors_bias = list(prior_weightfunction(distribution = "two.sided", parameters =
    list(alpha = c(1, 1), steps = c(0.05)), prior_weights = 1/12),
    prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1,
    1), steps = c(0.05, 0.1)), prior_weights = 1/12), prior_weightfunction(distribution =
    "one.sided", parameters = list(alpha = c(1, 1), steps = c(0.05)), prior_weights =
    1/12), prior_weightfunction(distribution = "one.sided", parameters = list(alpha =
    c(1, 1, 1), steps = c(0.025, 0.05)), prior_weights = 1/12), 
    
    prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1,
    1), steps = c(0.05, 0.5)), prior_weights = 1/12), prior_weightfunction(distribution =
    "one.sided", parameters = list(alpha = c(1, 1, 1, 1), steps = c(0.025, 0.05, 0.5)),
    prior_weights = 1/12), prior_PET(distribution = "Cauchy", parameters = list(0, 1),
    truncation = list(0, Inf), prior_weights = 1/4), prior_PEESE(distribution = "Cauchy",
    parameters = list(0, 5), truncation = list(0, Inf), prior_weights = 1/4)),
  priors_effect_null = prior(distribution = "point", parameters = list(location = 0)),
  priors_heterogeneity_null = prior(distribution = "point", parameters = list(location =
    0)),
  priors_bias_null = prior_none(),
  priors_hierarchical = prior("beta", parameters = list(alpha = 1, beta = 1)),
  priors_hierarchical_null = NULL,
  models = FALSE,
  silent = FALSE
)
```

## Arguments

- model_type:

  string specifying the RoBMA ensemble. Defaults to `NULL`. The other
  options are `"PSMA"`, `"PP"`, and `"2w"` which override settings
  passed to the `priors_effect`, `priors_heterogeneity`,
  `priors_effect`, `priors_effect_null`, `priors_heterogeneity_null`,
  `priors_bias_null`, and `priors_effect`. See details for more
  information about the different model types.

- priors_effect:

  list of prior distributions for the effect size (`mu`) parameter that
  will be treated as belonging to the alternative hypothesis. Defaults
  to a standard normal distribution
  `prior(distribution = "normal", parameters = list(mean = 0, sd = 1))`.

- priors_heterogeneity:

  list of prior distributions for the heterogeneity `tau` parameter that
  will be treated as belonging to the alternative hypothesis. Defaults
  to
  `prior(distribution = "invgamma", parameters = list(shape = 1, scale = .15))`
  that is based on heterogeneities estimates from psychology (van Erp et
  al. 2017) .

- priors_bias:

  list of prior distributions for the publication bias adjustment
  component that will be treated as belonging to the alternative
  hypothesis. Defaults to
  `list( prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1), steps = c(0.05)), prior_weights = 1/12), prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1, 1), steps = c(0.05, 0.10)), prior_weights = 1/12), prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1), steps = c(0.05)), prior_weights = 1/12), prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1), steps = c(0.025, 0.05)), prior_weights = 1/12), prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1), steps = c(0.05, 0.5)), prior_weights = 1/12), prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1, 1), steps = c(0.025, 0.05, 0.5)), prior_weights = 1/12), prior_PET(distribution = "Cauchy", parameters = list(0,1), truncation = list(0, Inf), prior_weights = 1/4), prior_PEESE(distribution = "Cauchy", parameters = list(0,5), truncation = list(0, Inf), prior_weights = 1/4) )`,
  corresponding to the RoBMA-PSMA model introduce by Bartoš et
  al. (2023) .

- priors_effect_null:

  list of prior distributions for the effect size (`mu`) parameter that
  will be treated as belonging to the null hypothesis. Defaults to a
  point null hypotheses at zero,
  `prior(distribution = "point", parameters = list(location = 0))`.

- priors_heterogeneity_null:

  list of prior distributions for the heterogeneity `tau` parameter that
  will be treated as belonging to the null hypothesis. Defaults to a
  point null hypotheses at zero (a fixed effect meta-analytic models),
  `prior(distribution = "point", parameters = list(location = 0))`.

- priors_bias_null:

  list of prior weight functions for the `omega` parameter that will be
  treated as belonging to the null hypothesis. Defaults no publication
  bias adjustment, [`prior_none()`](reference/prior_none.md).

- priors_hierarchical:

  list of prior distributions for the correlation of random effects
  (`rho`) parameter that will be treated as belonging to the alternative
  hypothesis. This setting allows users to fit a hierarchical
  (three-level) meta-analysis when `study_ids` are supplied. Note that
  this is an experimental feature and see News for more details.
  Defaults to a beta distribution
  `prior(distribution = "beta", parameters = list(alpha = 1, beta = 1))`.

- priors_hierarchical_null:

  list of prior distributions for the correlation of random effects
  (`rho`) parameter that will be treated as belonging to the null
  hypothesis. Defaults to `NULL`.

- models:

  should the models' details be printed.

- silent:

  do not print the results.

## Value

`check_setup` invisibly returns list of summary tables.

## See also

[`check_setup.reg()`](reference/check_setup.reg.md)
[`RoBMA()`](reference/RoBMA.md)

## Examples

``` r
# check default RoBMA setup
check_setup()
#> Robust Bayesian meta-analysis (set-up)
#> Components summary:
#>               Models Prior prob.
#> Effect         18/36       0.500
#> Heterogeneity  18/36       0.500
#> Bias           32/36       0.500

# check setup with custom priors
check_setup(
  priors_effect = prior("normal", list(mean = 0, sd = 0.5)),
  priors_heterogeneity = prior("invgamma", list(shape = 2, scale = 0.3)),  
  models = TRUE
)
#> Robust Bayesian meta-analysis (set-up)
#> Components summary:
#>               Models Prior prob.
#> Effect         18/36       0.500
#> Heterogeneity  18/36       0.500
#> Bias           32/36       0.500
#> 
#> Models overview:
#>  Model  Prior Effect  Prior Heterogeneity
#>      1       Spike(0)            Spike(0)
#>      2       Spike(0)            Spike(0)
#>      3       Spike(0)            Spike(0)
#>      4       Spike(0)            Spike(0)
#>      5       Spike(0)            Spike(0)
#>      6       Spike(0)            Spike(0)
#>      7       Spike(0)            Spike(0)
#>      8       Spike(0)            Spike(0)
#>      9       Spike(0)            Spike(0)
#>     10       Spike(0)    InvGamma(2, 0.3)
#>     11       Spike(0)    InvGamma(2, 0.3)
#>     12       Spike(0)    InvGamma(2, 0.3)
#>     13       Spike(0)    InvGamma(2, 0.3)
#>     14       Spike(0)    InvGamma(2, 0.3)
#>     15       Spike(0)    InvGamma(2, 0.3)
#>     16       Spike(0)    InvGamma(2, 0.3)
#>     17       Spike(0)    InvGamma(2, 0.3)
#>     18       Spike(0)    InvGamma(2, 0.3)
#>     19 Normal(0, 0.5)            Spike(0)
#>     20 Normal(0, 0.5)            Spike(0)
#>     21 Normal(0, 0.5)            Spike(0)
#>     22 Normal(0, 0.5)            Spike(0)
#>     23 Normal(0, 0.5)            Spike(0)
#>     24 Normal(0, 0.5)            Spike(0)
#>     25 Normal(0, 0.5)            Spike(0)
#>     26 Normal(0, 0.5)            Spike(0)
#>     27 Normal(0, 0.5)            Spike(0)
#>     28 Normal(0, 0.5)    InvGamma(2, 0.3)
#>     29 Normal(0, 0.5)    InvGamma(2, 0.3)
#>     30 Normal(0, 0.5)    InvGamma(2, 0.3)
#>     31 Normal(0, 0.5)    InvGamma(2, 0.3)
#>     32 Normal(0, 0.5)    InvGamma(2, 0.3)
#>     33 Normal(0, 0.5)    InvGamma(2, 0.3)
#>     34 Normal(0, 0.5)    InvGamma(2, 0.3)
#>     35 Normal(0, 0.5)    InvGamma(2, 0.3)
#>     36 Normal(0, 0.5)    InvGamma(2, 0.3)
#>                          Prior Bias                         Prior prob.
#>                                                                   0.125
#>            omega[two-sided: .05] ~ CumDirichlet(1, 1)             0.010
#>        omega[two-sided: .1, .05] ~ CumDirichlet(1, 1, 1)          0.010
#>            omega[one-sided: .05] ~ CumDirichlet(1, 1)             0.010
#>      omega[one-sided: .05, .025] ~ CumDirichlet(1, 1, 1)          0.010
#>        omega[one-sided: .5, .05] ~ CumDirichlet(1, 1, 1)          0.010
#>  omega[one-sided: .5, .05, .025] ~ CumDirichlet(1, 1, 1, 1)       0.010
#>                              PET ~ Cauchy(0, 1)[0, Inf]           0.031
#>                            PEESE ~ Cauchy(0, 5)[0, Inf]           0.031
#>                                                                   0.125
#>            omega[two-sided: .05] ~ CumDirichlet(1, 1)             0.010
#>        omega[two-sided: .1, .05] ~ CumDirichlet(1, 1, 1)          0.010
#>            omega[one-sided: .05] ~ CumDirichlet(1, 1)             0.010
#>      omega[one-sided: .05, .025] ~ CumDirichlet(1, 1, 1)          0.010
#>        omega[one-sided: .5, .05] ~ CumDirichlet(1, 1, 1)          0.010
#>  omega[one-sided: .5, .05, .025] ~ CumDirichlet(1, 1, 1, 1)       0.010
#>                              PET ~ Cauchy(0, 1)[0, Inf]           0.031
#>                            PEESE ~ Cauchy(0, 5)[0, Inf]           0.031
#>                                                                   0.125
#>            omega[two-sided: .05] ~ CumDirichlet(1, 1)             0.010
#>        omega[two-sided: .1, .05] ~ CumDirichlet(1, 1, 1)          0.010
#>            omega[one-sided: .05] ~ CumDirichlet(1, 1)             0.010
#>      omega[one-sided: .05, .025] ~ CumDirichlet(1, 1, 1)          0.010
#>        omega[one-sided: .5, .05] ~ CumDirichlet(1, 1, 1)          0.010
#>  omega[one-sided: .5, .05, .025] ~ CumDirichlet(1, 1, 1, 1)       0.010
#>                              PET ~ Cauchy(0, 1)[0, Inf]           0.031
#>                            PEESE ~ Cauchy(0, 5)[0, Inf]           0.031
#>                                                                   0.125
#>            omega[two-sided: .05] ~ CumDirichlet(1, 1)             0.010
#>        omega[two-sided: .1, .05] ~ CumDirichlet(1, 1, 1)          0.010
#>            omega[one-sided: .05] ~ CumDirichlet(1, 1)             0.010
#>      omega[one-sided: .05, .025] ~ CumDirichlet(1, 1, 1)          0.010
#>        omega[one-sided: .5, .05] ~ CumDirichlet(1, 1, 1)          0.010
#>  omega[one-sided: .5, .05, .025] ~ CumDirichlet(1, 1, 1, 1)       0.010
#>                              PET ~ Cauchy(0, 1)[0, Inf]           0.031
#>                            PEESE ~ Cauchy(0, 5)[0, Inf]           0.031

# show detailed model overview
check_setup(models = TRUE)
#> Robust Bayesian meta-analysis (set-up)
#> Components summary:
#>               Models Prior prob.
#> Effect         18/36       0.500
#> Heterogeneity  18/36       0.500
#> Bias           32/36       0.500
#> 
#> Models overview:
#>  Model Prior Effect Prior Heterogeneity
#>      1     Spike(0)            Spike(0)
#>      2     Spike(0)            Spike(0)
#>      3     Spike(0)            Spike(0)
#>      4     Spike(0)            Spike(0)
#>      5     Spike(0)            Spike(0)
#>      6     Spike(0)            Spike(0)
#>      7     Spike(0)            Spike(0)
#>      8     Spike(0)            Spike(0)
#>      9     Spike(0)            Spike(0)
#>     10     Spike(0)   InvGamma(1, 0.15)
#>     11     Spike(0)   InvGamma(1, 0.15)
#>     12     Spike(0)   InvGamma(1, 0.15)
#>     13     Spike(0)   InvGamma(1, 0.15)
#>     14     Spike(0)   InvGamma(1, 0.15)
#>     15     Spike(0)   InvGamma(1, 0.15)
#>     16     Spike(0)   InvGamma(1, 0.15)
#>     17     Spike(0)   InvGamma(1, 0.15)
#>     18     Spike(0)   InvGamma(1, 0.15)
#>     19 Normal(0, 1)            Spike(0)
#>     20 Normal(0, 1)            Spike(0)
#>     21 Normal(0, 1)            Spike(0)
#>     22 Normal(0, 1)            Spike(0)
#>     23 Normal(0, 1)            Spike(0)
#>     24 Normal(0, 1)            Spike(0)
#>     25 Normal(0, 1)            Spike(0)
#>     26 Normal(0, 1)            Spike(0)
#>     27 Normal(0, 1)            Spike(0)
#>     28 Normal(0, 1)   InvGamma(1, 0.15)
#>     29 Normal(0, 1)   InvGamma(1, 0.15)
#>     30 Normal(0, 1)   InvGamma(1, 0.15)
#>     31 Normal(0, 1)   InvGamma(1, 0.15)
#>     32 Normal(0, 1)   InvGamma(1, 0.15)
#>     33 Normal(0, 1)   InvGamma(1, 0.15)
#>     34 Normal(0, 1)   InvGamma(1, 0.15)
#>     35 Normal(0, 1)   InvGamma(1, 0.15)
#>     36 Normal(0, 1)   InvGamma(1, 0.15)
#>                          Prior Bias                         Prior prob.
#>                                                                   0.125
#>            omega[two-sided: .05] ~ CumDirichlet(1, 1)             0.010
#>        omega[two-sided: .1, .05] ~ CumDirichlet(1, 1, 1)          0.010
#>            omega[one-sided: .05] ~ CumDirichlet(1, 1)             0.010
#>      omega[one-sided: .05, .025] ~ CumDirichlet(1, 1, 1)          0.010
#>        omega[one-sided: .5, .05] ~ CumDirichlet(1, 1, 1)          0.010
#>  omega[one-sided: .5, .05, .025] ~ CumDirichlet(1, 1, 1, 1)       0.010
#>                              PET ~ Cauchy(0, 1)[0, Inf]           0.031
#>                            PEESE ~ Cauchy(0, 5)[0, Inf]           0.031
#>                                                                   0.125
#>            omega[two-sided: .05] ~ CumDirichlet(1, 1)             0.010
#>        omega[two-sided: .1, .05] ~ CumDirichlet(1, 1, 1)          0.010
#>            omega[one-sided: .05] ~ CumDirichlet(1, 1)             0.010
#>      omega[one-sided: .05, .025] ~ CumDirichlet(1, 1, 1)          0.010
#>        omega[one-sided: .5, .05] ~ CumDirichlet(1, 1, 1)          0.010
#>  omega[one-sided: .5, .05, .025] ~ CumDirichlet(1, 1, 1, 1)       0.010
#>                              PET ~ Cauchy(0, 1)[0, Inf]           0.031
#>                            PEESE ~ Cauchy(0, 5)[0, Inf]           0.031
#>                                                                   0.125
#>            omega[two-sided: .05] ~ CumDirichlet(1, 1)             0.010
#>        omega[two-sided: .1, .05] ~ CumDirichlet(1, 1, 1)          0.010
#>            omega[one-sided: .05] ~ CumDirichlet(1, 1)             0.010
#>      omega[one-sided: .05, .025] ~ CumDirichlet(1, 1, 1)          0.010
#>        omega[one-sided: .5, .05] ~ CumDirichlet(1, 1, 1)          0.010
#>  omega[one-sided: .5, .05, .025] ~ CumDirichlet(1, 1, 1, 1)       0.010
#>                              PET ~ Cauchy(0, 1)[0, Inf]           0.031
#>                            PEESE ~ Cauchy(0, 5)[0, Inf]           0.031
#>                                                                   0.125
#>            omega[two-sided: .05] ~ CumDirichlet(1, 1)             0.010
#>        omega[two-sided: .1, .05] ~ CumDirichlet(1, 1, 1)          0.010
#>            omega[one-sided: .05] ~ CumDirichlet(1, 1)             0.010
#>      omega[one-sided: .05, .025] ~ CumDirichlet(1, 1, 1)          0.010
#>        omega[one-sided: .5, .05] ~ CumDirichlet(1, 1, 1)          0.010
#>  omega[one-sided: .5, .05, .025] ~ CumDirichlet(1, 1, 1, 1)       0.010
#>                              PET ~ Cauchy(0, 1)[0, Inf]           0.031
#>                            PEESE ~ Cauchy(0, 5)[0, Inf]           0.031
```
